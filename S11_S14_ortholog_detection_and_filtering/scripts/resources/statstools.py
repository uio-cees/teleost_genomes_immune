"""
A reusable module of statistical tools.  Originally included in Concaterpillar.

"""

import random, math, scipy.stats, scipy.special

VERYSMALL = 10e-6

def calcBLratio(seqidA, seqidB):
  """Calculate a likelihood ratio for the branch-length congruence test.

  Parameters:
    seqidA -- ID number for set A
    seqidB -- ID number for set B
    
  Return:
    Tuple containing the likelihood ratio and degrees of freedom.

  """
  
  likes = []
  nparams = []
  
  setfiles = ('set%03d' % seqidA, 'set%03d' % seqidB, 'set%03d%03d' \
              % (seqidA, seqidB))
  
  for set in setfiles:
    lnlfile = open('%s.lnl' % set)
    likes.append(float(lnlfile.readline().strip()))
    lnlfile.close()
      
    paramfile = open('%s.prm' % set)
    nparams.append(int(paramfile.readline().strip()))
    paramfile.close()
       
  # lnL^A + lnL^B - lnL^{AB}
  result = likes[0] + likes[1] - likes[2]
  df = nparams[0] + nparams[1] - nparams[2]
        
  # this should only happen in there are 2 copies of an alignment, I think
  if result <= 0:
    sys.stderr.write('dlnL is less than or equal to 0! Converting to 1E-6.\n')
    result = 1e-6
  
  return (result, df)

def correctalpha(desiredalpha, level):
  """Correct the user-defined alpha level for the number of hierarchical levels.
  
  This multiple-comparison correction serves to lower the alpha level according
  to the equation:
  
  alpha_c = 1 - (1 - alpha_u)^(1/max(levels_p, level)
  
  where alpha_c is the corrected alpha, alpha_u is the user-specified alpha,
  levels_p is the predicted number of hierarchical levels, and level is the
  current level of the hierarchy.
  
  However, since the predicted number of levels will always be at most the
  current level, the maximum max(levels_p, level) will always be equal to level.
  So this equation becomes:
  
  alpha_c = 1 - (1 - alpha_u)^level
  
  Parameters:
    desiredalpha -- user-specified alpha level
    level -- current level of the test
    
  Return:
    A float, the corrected alpha level.
    
  """
  
  correctedalpha = 1 - (1 - desiredalpha) ** (1.0 / level)
  
  return correctedalpha

def arecloseenough(x1, x2):
  """Return True iff difference is VERYSMALL or less."""

  if abs(x1 - x2) <= VERYSMALL:
    return True
  
  return False

def firstderivative(func, x, samples):
  """Estimate the first derivative of the function at x.
  
  Parameters:
    func -- function for which to calculate first derivative.  Must take 2 
            arguments, the first of which should be a float (will be about x),
            second can really be anything (samples)
    x -- point at which to estimate the first derivative (should be a float)
    samples -- other argument to pass to func... can really be anything

  Return:
    Estimated first derivative (a float, probably).
  
  """
  
  a = 0.5 * VERYSMALL 
  
  fxmina = func(x - a, samples)
  fxplusa = func(x + a, samples)
  
  return (fxplusa - fxmina) / (2 * a)

def secondderivative(func, x, samples):
  """Estimate the second derivative of the function at x.
  
  Parameters:
    func -- function for which to calculate second derivative.  Must take 2 
            arguments, the first of which should be a float (will be about x)
    x -- point at which to estimate the first derivative (should be a float)
    samples -- other argument to pass to func... can really be anything

  Return:
    Estimated second derivative (a float, probably).
  
  """
  
  a = VERYSMALL
  
  fxplus2a = func(x + 2 * a, samples)
  fxplusa = func(x + a, samples)
  fx = func(x, samples)
  
  return (fxplus2a - 2 * fxplusa + fx) / (a ** 2)  

def optimiselike(likefunc, samples, x1):
  """Find the optimum (maximum) of likefunc by the Newton-Raphson method.
  
  Parameters:
    likefunc -- function whose optimum is to be calculated. Must take 2 
                arguments, the first of which will start as x1 (must be a float)
                the second (samples) can be anything, but will probably be a
                list of numbers.
    samples -- second argument that will eventually be passed to likefunc
    x1 -- starting estimate of the optimum of likefunc
    
  Return:
    A float, the value of x at the optimum of likefunc.
  
  """

  xnplus1 = x1
  xn = x1 + 100 # makes starting condition true
  
  while not arecloseenough(xnplus1, xn):
    xn = xnplus1    
    xnplus1 = xn - (firstderivative(likefunc, xn, samples)/secondderivative(likefunc, xn, samples))
    
  return xn

def getlikeweibull(k, samples):
  """Calculate the likelihood of shape parameter k of a Weibull distribution.
  
  Other parameters of the distribution are estimated based on k and sample
  mean from samples.
  
  Parameters:
    k -- shape parameter for the distribution
    samples -- list of floating point values assumed to follow a Weibull
    
  Return:
    likelihood of shape parameter k.
  
  """
  
  N = len(samples)
  samplemean = sum(samples) / N
  
  # inverse of the scale parameter
  lambdainv = scipy.special.gamma(1 + 1/k) / samplemean 
  
  # equation for likelihood:
  # Nlog(k/lambda) + sum{(k-1)log(x_i/lambda) - (x_i/lambda)^k}
  
  sumterm = 0
  for val in samples:
    sumterm += ((k - 1) * math.log(val * lambdainv, math.e) - (val * lambdainv) ** k)
  
  # log-likelihood
  like = N * math.log(k * lambdainv, math.e) + sumterm
  
  return like
                                      
def getpval(teststat, statlist):
  """Calculate the prob. of teststat if drawn from distribution of statlist.
  
  Statlist is fit to a Weibull distribution, and the p-value for teststat
  is estimated from this distribution.
  
  Parameters:
    teststat -- some test statistic (a float)
    statlist -- a set of teststatistics generated under the null hypothesis.
    
  Return:
    The probability that teststat comes from the same distribution as statlist
    (a float, between 0 and 1)
  
  """
  
  propzero = 0
  bootvals = []
  for val in statlist:
    if val == 0:
      propzero += 1
    else:
      bootvals.append(val)
      
  propzero = float(propzero) / len(statlist)
  
  shapeinit = getstartingshape(statlist)
  
  shape = optimiselike(getlikeweibull, bootvals, shapeinit)
  scale = (sum(bootvals) / len(bootvals)) / scipy.special.gamma(1 + 1/shape)
  
  pvalue = math.exp(- (teststat/scale) ** shape)
  
  return pvalue * (1 - propzero)
  
def getstartingshape(vals):
  """Estimate a shape parameter for a Weibull distribution fitted to vals.
  
  This should be a starting value, in order to calculate the ML value of the
  shape.  Ideally, a good starting value should be calculated, but currently
  this function is a 'stub' that just returns 1.
  
  Parameters:
    vals -- a list of floats to be fit to a Weibull distribution
    
  Return:
    An estimate of the shape parameter for the Weibull distribution that fits
    vals best.
  
  """
  
  return 1

def getrawpval(teststat, statlist):
  """Calculate the prob. of teststat if drawn from distribution of statlist.
  
  If the test statistic is smaller than all bootstrap statistics, a value of 1
  will be returned.  If the statistic is larger than all bootstraps, a value of
  0 will be returned.  If the statistic is equal to 1 or more bootstrap
  statistics, a random number is chosen in the range of the possible corres-
  ponding p-values.
  
  Parameters:
    teststat -- some test statistic (a float)
    statlist -- a set of teststatistics generated under the null hypothesis.
    
  Return:
    The probability that teststat comes from the same distribution as statlist
    (a float, between 0 and 1)

  """
  
  nstats = len(statlist)
  highend = lowend = 0 # specifies the range of the p-value
  
  for i in range(nstats):
    
    if teststat == statlist[i]:
      lowend = 1 - (i + 1)/float(nstats)
      if highend < lowend:
        highend = 1 - i/float(nstats)
        
    elif teststat < statlist[i]:
      if highend == 0:
        return 1 - i/float(nstats)
      
  return random.uniform(lowend, highend)

def pchisq(x, df):
  """Calculate a prob. of drawing x from a chi-square distribution.
  
  Parameters:
    x -- a test statistic that follows a chi-square distribution
    df -- number of degrees of freedom for the distribution
    
  Return:
    The probability that x was drawn from a chi-squre distribution with df 
    degrees of freedom

  """
  
  if df % 2 == 0:
    dchi = 0.5 * math.exp(-0.5 * x)
    f = 1.0 - 2.0 * dchi
    for i in range(4, df + 1, 2):
      dchi *= x / (i - 2)
      f -= 2.0 * dchi
      
  else:
    f = 2.0 * pnorm(math.sqrt(x), 0.0, 1.0) - 1.0
    dchi = math.exp(-0.5 * x) / math.sqrt(2.0 * math.pi * x)
    for i in range(3, df + 1, 2):
      dchi *= x / (i - 2)
      f -= 2.0 * dchi
      
  return f

def pnorm(x, a, b):
  """No idea what this does, or what these args are."""
  rt2 = 1.41421356237309515
  xhi = 2.66e+1
  
  z = (x - a)/(rt2 * b)
  xt = abs(z)

  if xt >= xhi:
    if x < 0.0:
      return 0.0
    else:
      return 1.0
    
  t = 1.0 - 7.5 / (xt + 3.75)
  y = (((((((((((((((3.328130055126039e-10 * \
                     t - 5.718639670776992e-10) * t - 4.066088879757269e-9) * \
                   t + 7.532536116142436e-9) * t + 3.026547320064576e-8) * \
                 t - 7.043998994397452e-8) * t - 1.822565715362025e-7) *  \
               t + 6.575825478226343e-7) * t + 7.478317101785790e-7) * \
             t - 6.182369348098529e-6) * t + 3.584014089915968e-6) * \
           t + 4.789838226695987e-5) * t - 1.524627476123466e-4) * \
         t - 2.553523453642242e-5) * t + 1.802962431316418e-3) * \
       t - 8.220621168415435e-3) * t + 2.414322397093253e-2
  
  y = ((((( y * t - 5.480232669380236e-2) * t + 1.026043120322792e-1) * \
         t - 1.635718955239687e-1) * t + 2.260080669166197e-1) * \
       t -2.734219314954260e-1) * t + 1.455897212750385e-1
  
  if z > 0.0:
    return 1.0 - 0.5 * math.exp(-xt * xt) * y
  
  return 0.5 * math.exp(-xt * xt) * y

