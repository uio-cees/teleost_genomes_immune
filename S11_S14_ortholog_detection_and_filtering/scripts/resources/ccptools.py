"""
Functions required by Concaterpillar.

"""

from filetools import *
from seqtools import *
from statstools import *
from treetools import *
# get RAxML location and version
from concaterpillar import raxmlpath
raxmlversion = getraxmlversion(raxmlpath)

import cPickle, getopt, glob, re, sys, shutil, threading, subprocess
import Queue

mpisize = 10
mpipresent = False
try:
  import mpi
  if mpi.size > 1:
    mpipresent = True
except ImportError:
  pass


def brlensetup(model, debug=False):
  """Move to the directory BLtest, if it exists, and copy relevant files.

  If the directory BLtest does not exist, it is created.  If the file 
  catsets.pickle exists (from the congruence test), the largest concatenated
  set is assumed to be the one on which the branch length compatibility test
  should be performed.  The sequences from this set are copied into BLtest,
  along with the tree from their concatenation.  If catsets.pickle does not
  exist, all sequence files are copied over, along with the tree file global.tre
  (if it exists), which is assumed to be the tree from concatenation of all
  alignments.  If it doesn't exist, alignments are concatenated, and RAxML is
  run on the global alignment to produce the global tree.
  
  Parameters:
    model -- substitution model to use.  Must be 'DAYHOFF', 'DCMUT', 'JTT', 
             'MTREV', 'WAG', 'RTREV', 'CPREV', 'VT', 'BLOSUM62', 'MTMAM', 'GTR'
  
  Return:
    List containing names of sequence files copied to the BLtest directory
  
  """
  
  bltestdirs = glob.glob('BLtest*')
  
  try:
    catsets = recover(model=model, rebuild=False)

    setcount = 0
    for cs in catsets:
      bldir = 'BLtest-set%03d' % setcount 
      os.mkdir(bldir)
      shutil.copy('set%03d.tre' % setcount, '%s/global.tre' % bldir)
      
      for gene in cs:
        shutil.copy('%s' % gene, '%s/%s' % (bldir, gene))

      setcount += 1
      
  except IOError:
    print '\nCouldn\'t find backup from previous concaterpillar run. ' + \
          'Branch length congruence\ntest will be performed on all sequences!\n'      
    sys.stdout.flush()

    # make a "fake" catsets, with all genes in same topologically congruent set
    catsets = [getseqlist()] 
    if 'global.tre' not in os.listdir('.'):
      writephy(catseq(getaligns(catsets[0])), 'global.seq')
      buildmultitrees(('global.seq',), model, debug=debug)

  return catsets    
    
if not mpipresent:
  def brlentest(catsets, setdir, cutoff, model, nprocs = 1, debug=False):
    """Execute the branch-length congruence test.

    Parameters:
      catsets -- list of topologically congruent sets (produced by toptest)
      setdir -- "name" for branch-length test run (e.g. BLtest-set000/)
      cutoff -- uncorrected alpha-level for the test
      nprocs -- number of processes to thread off... should be number of 
                processors to use on a shared memory machine
    
    """
    
    resfile = open('%s.ccp' % setdir[:-1].lower(), 'w')
    gammacats = 4
    lastpval = 0
    q = Queue.Queue()
    lock = threading.Lock()
    threadcounter = Counter()
    
    treefile = open('global.tre')
    tree = treefile.readline().strip()
    treefile.close()
    
    nsets = 0
    for cs in catsets:
      nsets += len(cs)
      
    # current level is number of concatenations that have alredy happened (single-
    # gene sets minus remaining sets plus 1... note incrementation in loop below)
    level = nsets - len(catsets)
  
    while len(catsets) > 1:
      level += 1
      dlist = os.listdir('.')
  
      for i in range(len(catsets)):
        if 'set%03d.lnl' % i not in dlist:
          q.put({'tree': tree, 'seqfile': 'set%03d.seq' % i, \
               'gammacats': gammacats, 'model':model, 'debug':debug})
  
        for j in range(i + 1, len(catsets)):
          catsetname = 'set%03d%03d.seq' % (i, j)
  
          # if concatenated file doesn't exist, create it
          if catsetname not in dlist: 
            writephy(catseq(getaligns(('set%03d.seq' % i, 'set%03d.seq' % j))), 
                     catsetname)
          
          if 'set%03d%03d.lnl' % (i, j) not in dlist:
            q.put({'tree': tree, 'seqfile': catsetname, 'gammacats': gammacats, \
                   'model':model, 'debug':debug})
      
      threads = []
      for i in range(nprocs):
        t = PhyloThread(q, lock, target=runpuzzlebrlen, name='thread %d' % i,
                        counter=threadcounter)
        threads.append(t)
        t.start()
        
        # wait until thread really starts before starting next thread
        lock.acquire()
        lock.release()
        
      for t in threads:
        t.join()
        
      pvalmat = []
      
      for i in range(len(catsets)):
        pvalmat.append([])
        for j in range(i + 1, len(catsets)):
          result = calcBLratio(i, j)
          pvalmat[-1].append(1 - pchisq(2 * result[0], result[1]))
      
      row, col = getmax(pvalmat)
      setA = row
      setB = col + row + 1
      
      pval = pvalmat[row][col]      
      
      # Correct alpha based on current level... see argument in toptest
      alpha = correctalpha(cutoff, level)
      
      # Only tell the user we're correcting if we NEED to correct
      if pval < cutoff:
        print 'Applying alpha-level correction: new alpha = %.6f' % alpha
        resfile.write('Applying alpha-level correction: new alpha = %.6f\n' % alpha)
        sys.stdout.flush()
        resfile.flush()
      
      if pval < alpha:
        print 'Failed to concatenate %s and %s (p = %.6f).\n' % \
              (str(catsets[setA]), str(catsets[setB]), pval)
        sys.stdout.flush()
        resfile.write('\nFailed to concatenate %s and %s (p = %.6f).\n\n' % \
                      (str(catsets[setA]), str(catsets[setB]), pval))
        resfile.flush()
        break
        
      else: # p-val of chi-square test >= corrected alpha
        print 'Concatenating %s and %s (p = %.6f).\n' % \
              (str(catsets[setA]), str(catsets[setB]), pval)
        sys.stdout.flush()
        resfile.write('Concatenating %s and %s (p = %.6f).\n' % \
                      (str(catsets[setA]), str(catsets[setB]), pval))
        resfile.flush()
        
        # rename and remove necessary sets
        filedict = {}
        for i in ('.seq', '.lnl', '.prm'):
          filedict['set%03d%03d%s' % (setA, setB, i)] = 'set%03d%s' % (setA, i)
        renamefiles(filedict)
  
        # necessary only if setB is highest-numbered set
        removefiles(glob.glob('set%03d*' % setB))
        removesets(setA)
        removesets(setB)
        renamesets(setB)
        
        catsets[setA].extend(catsets[setB])
        del catsets[setB]
        
    print 'Finished concatenating sets.  The following files contain ' + \
          'the final sets:'
    resfile.write('Finished concatenating sets.  The following files contain ' + \
                  'the final sets:\n')
  
    for i in range(len(catsets)):
      print '%sset%03d.seq (%d genes): %s' % (setdir, i, len(catsets[i]), 
                                                   str(catsets[i]))
      resfile.write('%sset%03d.seq (%d genes): %s\n' % (setdir, i, 
                                                         len(catsets[i]), 
                                                         str(catsets[i])))
    resfile.close()
    
    print ''
    sys.stdout.flush()
  
def buildmatrix(likes):
  """Build a matrix of log-likelihood ratios from a matrix of log-likelihoods.

  The matrix returned is in upper-triangular form, and contains likelihood
  ratios of the form lnL^A_{T_A} + lnL^B_{T_B} - (lnL^A_{T_AB} + lnL^B_{T_AB})

  Parameters:
    likes -- a matrix of pairs of log-likelihoods [lnL^A_{T_A}, lnL^A_{T_AB}]
             with the diagonal missing (not invalid values, these cells are
             just absent)
             
  Return:
    Log-likelihood ratio matrix.
    
  """

  difmat = []
  nsets = len(likes)

  for i in range(nsets):
    row = []
    for j in range(1, nsets):
      # RAxML may estimate a better tree for alignment A from alignment AB
      # (basically, this constrains differences to be at least 0)
      row.append(max(likes[i][j - 1]) - likes[i][j - 1][1])
    difmat.append(row)
    
  for i in range(len(difmat)):
    for j in range(len(difmat[i]) + 1):
      if i == j:
        difmat[i][i:i] = [0.0]

  # add across the diagonal to get the likelihood ratio including both genes
  symmat = []
  for i in xrange(nsets):
    symmat.append([])
    for j in xrange(i):
      symmat[-1].append(difmat[i][j] + difmat[j][i])
    
  return symmat


if not mpipresent:
  def buildmultitrees(infiles, model, nprocs = 1, debug=False):
    """Run RAxML on a list of alignment filenames.
  
    Parameters:
      infiles -- list of file names on which to run RAxML
      nprocs -- number of processes to thread off... should be number of 
                processors to use on a shared memory machine
  
    """
    
    print 'I\'m about to build trees for %d data sets (using RAxML).' % (len(infiles))
    sys.stdout.flush()
    
    q = Queue.Queue()
    for i in infiles:
      cmds = {'model': model, 'seqfile': i, 'debug':debug}
      q.put(cmds)
    
    lock = threading.Lock()
    threadcounter = Counter()
    threads = []
    for i in range(nprocs):
      t = PhyloThread(name = 'thread %d' % i, jobqueue = q, 
                      lock = lock, target = buildtrees, counter = threadcounter)
      threads.append(t)
      t.start()
      
      # wait until new thread really starts before starting next thread
      lock.acquire()
      lock.release()
      
    for t in threads:
      t.join()
    
def buildtrees(cmds):
  """Run RAxML with the parameters specified.
  
  Parameters:
    cmds -- dictionary specifying the sequence file and model to use with RAxML.
  
  """
  
  raxmlpath = globals()['raxmlpath']
  raxmlversion = globals()['raxmlversion']
  infile = cmds['seqfile']
  suffix = infile[:-4]
  
  # try to run RAxML up to 5 times... this is only to deal with our sick cluster
  numtries = 5
  while numtries > 0:
    runraxml(raxmlpath, infile, suffix, cmds['model'], version=raxmlversion, debug=cmds['debug'])
    likefile = open('%s.lnl' % suffix, 'w')
    try:
      like = getlikelihood('RAxML_info.%s' % suffix)
      break
    except ConcaterpillarError:
      sys.stderr.write('Somehow, RAxML failed to complete. Trying again...\n')
    numtries -= 1
  
  if numtries == 0:
    sys.stderr.write('Failed to run RAxML. Is the raxmlpath variable correctly set (in concaterpillar.py)?\n')
    sys.exit(1)
          
  if like == None:
    sys.stderr.write('There was a problem running RAxML with file %s.\n' % 
                     infile)
    sys.stderr.write('Please check the file and try again.\n\n')
    sys.exit(1)
    
  #alphafile.write('%.6f\n' % alpha)
  #alphafile.close()
  likefile.write('%.9f\n' % like)
  likefile.close()
  
  renamefiles({'RAxML_result.%s' % suffix: suffix + '.tre'})
  
  for f in glob.glob('RAxML_*.%s' % suffix):
    os.remove(f)
  try:
    os.remove('%s.reduced' % infile)
  except OSError:
    pass

  print 'Finished analysis of %s.' % (infile)
  sys.stdout.flush()

def calcBLratio(seqA, seqB):
  """Calculate the log-likelihood ratio for the branch-length test.
  
  Parameters:
    setA: number of the first set
    setB: number of the second set
    
  Return:
    Tuple containing the likelihood ratio and number of degrees of freedom
    
  """
  
  likes = {}
  nparams = {}
  
  setnames = {'A': 'set%03d' % seqA, 'B': 'set%03d' % seqB, \
              'AB': 'set%03d%03d' % (seqA, seqB)}
  
  
  for setkey in setnames:
    lnlfile = open('%s.lnl' % setnames[setkey])
    likes[setkey] = float(lnlfile.readline().strip())
    lnlfile.close()
    
    paramfile = open('%s.prm' % setnames[setkey])
    nparams[setkey] = int(paramfile.readline().strip())
    paramfile.close()
       
  result = likes['A'] + likes['B'] - likes['AB']
  df = nparams['A'] + nparams['B'] - nparams['AB']
   
  # this will likely only happen if two alignments are identical
  if result <= 0:
    sys.stderr.write('dlnL is less than or equal to 0! Converting to 1E-6.\n')
    result = 1e-6
  
  return (result, df)


def cleanup(ask=True):
  """Check for and remove files that interfere with concaterpillar."""
  
  # remove useless files 
  rmfiles = glob.glob('boot*.*')
  rmfiles.extend(glob.glob('RAxML_*'))
  rmfiles.extend(glob.glob('*.seq.reduced'))

  if rmfiles:
    sys.stderr.write('Found the following files in the working directory:\n')
  
    for f in rmfiles:
      sys.stderr.write('%s\n' % f)
    if ask:
      while True:
        reply = raw_input('Remove files or Quit (R or Q)? ').upper()
        if reply == 'Q':
          sys.exit(1)
        elif reply == 'R':
          break
        else:
          sys.stderr.write('Invalid response!\n')
    else:
      sys.stderr.write('Deleted files without asking.\n')
    for f in rmfiles:
      os.remove(f)
          
def getalllikes(nsets):
  """Return a matrix of likelihoods from RAxML output files.

  The RAxML output files should have been produced by evaluating the
  likelihood for a data set (set A) under 2 trees: the ML tree from set A, and 
  the ML tree from a data set produced by concatenating set A with a second set
  (set B).  Files need to be named set*$.puz, where * and $ are each 3-digit 
  numbers, * representing the number assigned to set A, and $ the number
  assigned to set B.  Sets must be numbered consecutively, starting with 000.

  Parameters:
    nsets -- number of alignments for which likelihoods were reestimated

  Return:
    a matrix of pairs of likelihoods [a, b], in which a is the likelihood of 
    data set A, given the ML (estimated) tree for A, and b is the likelihood
    of data set A, given the ML (estimated) tree for data set AB.  This is a
    square matrix, in which the diagonal is missing (i.e., matrix indices
    almost correspond to data set numbers, with the row representing the first
    number in the puzzle file (* above) and the column representing the second
    ($ above), but column index is off by one (too low), where column > row).
    
  """
  likes = []
  
  for i in range(nsets):
    likes.append([])
    for j in range(nsets):
      if i != j:
        likefile = open('set%03d%03d.lnl' % (i, j))
        likes[-1].append([float(line) for line in likefile])
        likefile.close()

  return likes


def getlikelihood(ofname):
  """Extract the GAMMA-likelihood value from a RAxML output file.
  
  Parameters:
    ofname -- the name of the RAxML output file
  
  Return:
    Extracted likelihood (a float)
  
  """
  
  raxmlversion = globals()['raxmlversion']
  
  try:
    outfile = open(ofname, 'r')
  except IOError:
    raise ConcaterpillarError('Couldn\'t find RAxML file %s!' % ofname)

  like = None
  
  if raxmlversion == 7.2:
    infoline = re.compile(r'GAMMA-based Likelihood: +(-\d+\.\d+)|GAMMA  likelihood: +(-\d+\.\d+)')
  elif raxmlversion == 7.0:
    infoline = re.compile(r'GAMMA[-\s]+likelihood:* +(-\d+\.\d+)\,*')
  for line in outfile:
    gamma = infoline.search(line)
    if gamma:
      results = gamma.groups()
      like = None
      for item in results:
        if item:
          try:
            like = float(item)
            
          except TypeError:
            sys.stderr.write('Line: %s' % line)
            raise
      if not like:
        raise ConcaterpillarError('Couldn\'t read likelihood from RAxML file %s!' % ofname)
      break
  outfile.close()
  
  if like == None:
    raise ConcaterpillarError('Could not find likelihood in file %s!' % ofname)
  
  return like



def getinput():
  """Display a menu and determine at which stage Concaterpillar should begin.
     
     Deprecated in Concaterpillar 1.5
  
  TODO: add model selection option!
  
  Return:
    Dictionary summarising values selected from the menu.
  
  """

  exstat = os.system('clear')
  if exstat != 0:
    os.system('cls')
                              
  print '\n\nWelcome to Concaterpillar!'
  print '\
                     o  o\
\n __  __  __  __  __  _\\/\
\n/  \\/  \\/  \\/  \\/  \\/ .\\\
\n\\__/\\__/\\__/\\__/\\__/\\_~/\
\n     ""  ""  ""  ""\n\n' 

  options = ('t', 'i', 'b', 'c', 'p', 's')
  menuitems = {'t': 'Start NEW topological congruence test', 
               'i': 'Resume interrupted topological congruence test', 
               'b': 'Branch-length congruence test', 
               'c': 'Number of CPUs to use',
               'p': 'Alpha level (p-value cutoff) for tests', 
               's': 'Save a backup point every n rounds (0 = not at all)' }
  
  itemstates = {'t': True, 'i': False, 'b': True,
                'c': 1, 'p': 0.05, 's': 0}
  
  statemap = {True: 'Y', False: 'N'}

  while True:
    print 'OPTIONS:\n'
    for opt in options:
      if isinstance(itemstates[opt], types.BooleanType):
        print '%3s%60s%6s' % (opt.upper(), menuitems[opt], \
                              statemap[itemstates[opt]])
      else:
        print '%3s%60s%6s' % (opt.upper(), menuitems[opt], str(itemstates[opt]))
    
    while True:
      userin = raw_input('Menu item to change, Y to continue, Q to quit: ').lower()
    
      if userin == 'q':
        sys.exit(0)
      
      if userin == 'y':
        optdict = {}
        for key in itemstates:
          if itemstates[key] is False:
            continue
          optdict['-%s' % key] = itemstates[key]
            
        return optdict
    
      if userin == 't':
        if itemstates['t']:
          itemstates['t'] = False
          break
            
        else:
          if itemstates['i']:
            print 'Turn off option I first!'
          else:
            itemstates['t'] = True
            break
          
      elif userin == 'i':
        if itemstates['i']:
          itemstates['i'] = False
          break
        else:
          if itemstates['t']:
            print 'Turn off option T first!'
          else:
            itemstates['i'] = True
          
      elif userin == 'b':
        if itemstates['b']:
          itemstates['b'] = False
          break
          
        else:
          itemstates['b'] = True
          break
          
      elif userin == 'c':
        nprocs = raw_input('How many CPUs do you want to use? ')
        try:
          nprocs = int(nprocs)
          if nprocs < 1:
            raise ValueError
          itemstates['c'] = nprocs
          break
        
        except ValueError:
          print 'Must use at least 1 processor!'
        
      elif userin == 'p':
        cutoff = raw_input('What cutoff do you wish to use? ')
        try:
          cutoff = float(cutoff)
          if cutoff < 0 or cutoff > 1:
            raise ValueError
          itemstates['p'] = cutoff
          break
        
        except ValueError:
          print 'Cutoff must be a value between 0 and 1!'
          
      elif userin =='s':
        saverounds = raw_input('after how many rounds do you want to save? ')
        try:
          saverounds = int(saverounds)
          if saverounds < 0:
            raise ValueError
          itemstates['s'] = saverounds
          break
      
        except ValueError:
          print 'The number of rounds must be a non-negative integer!'
      
      else:
        print 'Invalid option!'
    
def getlikes(ofname):
  """Extract and return a list of log-likelihoods from a TREE-PUZZLE outfile.

  Parameters:
    ofname -- name of the PUZZLE output file.
    
  Return:
    Extracted log-likelihoods.

  """

  try:
    outfile = open(ofname, 'r')
  except IOError:
    raise ConcaterpillarError('Couldn\'t find PUZZLE file %s!' % (ofname))

  likes = []
  
  lnlRE = re.compile(r'log L: (-\d+.\d+)')

  for line in outfile:
    result = lnlRE.search(line)
    if result:
      likes.append(float(result.group(1)))

  outfile.close()

  return likes

def getmin(matrix, blacklistedmat):
  """Find the minimum value in a 2-D matrix.
  
  Parameters:
    matrix -- an upper-triangular matrix of likelihood ratios.
    
  Return:
    Tuple containing (row, col) of the minimum value, (-1, -1) if the matrix
    is empty or has only 1 dimension.
  
  """
  row, col = 1, 0
  while blacklistedmat[row][col]:
    col += 1
    if col >= row:
      col = 0
      row += 1
    if row >= len(matrix):
      raise ConcaterpillarError('Matrix is completely blacklisted!')

  try:
    minval = matrix[row][col]
    
  except IndexError:
    return (-1,-1)

  for i in range(len(matrix)):
    for j in range(len(matrix[i])):
      if not blacklistedmat[i][j] and (matrix[i][j] < matrix[row][col]):
        row = i
        col = j

  return row, col

def getmax(matrix):
  """Find the maximum value in a 2-D matrix.
  
  Parameters:
    matrix -- an upper-triangular matrix of likelihood ratios.
    
  Return:
    Tuple containing (row, col) of the maximum value, (-1, -1) if the matrix
    is empty or has only 1 dimension.
  
  """
  try:
    row, col = 0, 0
    max = matrix[row][col]
  except IndexError:
    return (-1,-1)

  for i in range(len(matrix)):
    for j in range(len(matrix[i])):
      if matrix[i][j] > matrix[row][col]:
        row = i
        col = j

  return row, col


def getoptdict(argv):
  """Return a dictionary of command-line options, given command line args.
  
  Terminates the program with exit code 1 and a usage message if invalid options
  are specified.

  Parameters:
    argv -- sys.argv[1:] (the script's argument vector, minus the
          name of the script itself)
          
  Return:
    a dictionary containing the values for options with short option flags as
    keys
  
  """
  
  goodopts = 'hartibc:p:s:m:dn'
  goodopttup = ('help', 'phyml', 'puzzle', 'topology', 'interrupt', 'blength',
                 'cpus=', 'pval=', 'save=', 'model=', 'debug', 'noask')
  
  try:
    opts, args = getopt.getopt(argv, goodopts, goodopttup)
  except getopt.GetoptError:
    usage(1)
    
  if len(args) != 0:
    usage(1)
    
  equivalents = []
  
  optdict = dict(opts)
  
  equivalents = [('-h', '--help'), ('-a', '--phyml'), ('-r', '--puzzle'), 
                 ('-t', '--topology'), ('-i', '--interrupt'),
                 ('-b', '--blength'), ('-c', '--cpus'), ('-p', '--pval'),
                 ('-s', '--save'), ('-m', '--model'), ('-d', '--debug'), ('-n', '--noask')]
  
  goodoptdict = {}
  
  for opt in equivalents:
    if opt[0] in optdict and opt[1] in optdict:
      usage(1)
    if opt[1] in optdict:
      goodoptdict[opt[0]] = optdict[opt[1]]
    elif opt[0] in optdict:
      goodoptdict[opt[0]] = optdict[opt[0]]
  
  if '-a' in goodoptdict or '-r' in goodoptdict:
    sys.stderr.write('WARNING: the -a and -r flags have been deprecated. Please just use -t for the\n')
    sys.stderr.write('topology test, or -i to continue the topology test.\n')
  return goodoptdict

def getseqlist():
  """Produce a list of alignment files in the working directory.

  All files with the extension .seq other than those that start with 'seq' or 
  'boot', and not including 'global.seq' are included in the list.
  
  Return:
    List of alignment file names.

  """
  
  seqlist = []
    
  for f in glob.glob('*.seq'):
    if f[:3] == 'set' or f == 'global.seq' or f[:4] == 'boot':
      continue
    seqlist.append(f)

  seqlist.sort()
  
  return seqlist

if not mpipresent:
  def mlboot(args):
    """Determine the topology test statistic by ML using a bootstrapped dataset.
  
    Parameters:
      args -- a dictionary containing the following items:
        'lock': a threading.Lock object used to prevent deadlocking
        'startcount': a Counter used to determine iteration number
        'align': alignment from which sites are sampled
        'ncharA': number of sites to sample for the first alignment
        'ncharB': number of sites to sample for the second alignment
        'endcount': a Counter used to determine how many bootstrap replicates 
                    have finished
        'niter': indicates total number of replicates to perform
        'teststats': list to which the test statistic is appended
  
    """
    
    raxmlpath = globals()['raxmlpath']
    raxmlversion = globals()['raxmlversion']
    
    args['lock'].acquire()
    iterno = int(args['startcount'])
    args['startcount'].increment()
    args['lock'].release()
    
    # names of bootstrap files
    bootseqs = ['boot%03d%s.seq' % (iterno, c) for c in ('a', 'b', 'ab')]
    
    writephy(bootstrap(args['align'], args['ncharA']), bootseqs[0])
    writephy(bootstrap(args['align'], args['ncharB']), bootseqs[1])
    
    writephy(catseq(getaligns((bootseqs[0], bootseqs[1]))), bootseqs[2])
    
    for seq in bootseqs:
      suffix = seq[:-4]
      
      runraxml(raxmlpath, seq, suffix, model=args['model'], version=raxmlversion, debug=args['debug'])
      os.rename('RAxML_result.%s' % suffix, '%s.tre' % suffix)
      
      for raxfile in glob.glob('RAxML_*.%s' % suffix):
        os.remove(raxfile)
      try:
        os.remove('%s.reduced' % seq)
      except OSError:
        pass
          
    likes = []    
    trees = []
      
    for seq in bootseqs:
      treefile = open('%s.tre' % seq[:-4], 'r')
      trees.append(string.strip(treefile.readline()))
      treefile.close()
      
    
    for i in (0, 1):
      intree = '%s.intree' % (bootseqs[i][:-4])
      treefile = open(intree, 'w')
      treefile.write('%s\n%s\n' % (trees[i], trees[2]))
      treefile.close()
      # split trees for RAxML and make 2 runs
      treefile = open(intree, 'r')
      for letter in ('A', 'B'):
        suffix = '%s%s' % (bootseqs[i][:-4], letter)
        temptreename = '%s.temptree' % suffix
        temptree = open(temptreename, 'w')
        treeline = treefile.readline().strip()
        temptree.write(treeline)
        temptree.close()
        
        if args['model'] == 'GTR':
          binstring = '%s -n %s -m GTRGAMMA -s %s -f e -t %s' \
                      % (raxmlpath, suffix, bootseqs[i],temptreename)
        else:
          binstring = '%s -n %s -m PROTGAMMA%s -s %s -f e -t %s' \
                      % (raxmlpath, suffix, args['model'], bootseqs[i],temptreename)
        try:
          runbin(binstring, (), debug=args['debug'])
        except KeyError:
          runbin(binstring, ())
        likes.append(getlikelihood('RAxML_info.%s' % suffix))
    
    treefile.close()
    for f in glob.glob('*boot%03d*' % iterno):
      os.remove(f)
    
    args['lock'].acquire()
    
    # should be lnL(A)_{T_A} + lnL(B)_{T_B} - (lnL(A)_{T_AB} + lnL(B)_{T_AB}))
    # because lnL(A)_{T_AB} might be greater than lnL(A)_{T_A}, max taken instead
    args['teststats'].append(max(likes[:2]) - likes[1] + \
                             max(likes[2:]) - likes[3])
    
  # print pretty progress counter    
    countval = int(args['endcount'])
    args['endcount'].increment()
    if countval % 10 == 0:
      sys.stdout.write('[.')
    elif countval % 10 == 9:
      sys.stdout.write('.] %d/%d\n' % (countval + 1, args['niter']))
    else:
      sys.stdout.write('.')    
    sys.stdout.flush()
    
    args['lock'].release()


if not mpipresent:
  def multinonparboot(setA, setB, niter, model, nprocs, debug=False):
    """Perform several ML-based non-parametric bootstrap replicates  
    
    Parameters:
      setA -- number of first set
      setB -- number of second set
      niter -- number of times to repartition data set and get pseudo test stat
      nprocs -- number of processes to thread off... should be number of 
                processors to use on a shared memory machine
              
    Return:
      List of pseudo test statistics

    """
      
    teststats = []
    
    alignA = getseqs('set%03d.seq' % setA)
    alignB = getseqs('set%03d.seq' % setB)
    
    # want same # of taxa in bootstraps as in real dataset
    alignA, alignB = getsequnion(alignA, alignB)
        
    ncharA = alignA.length
    ncharB = alignB.length
    
    q = Queue.Queue()
    lock = threading.Lock()
    endcount = Counter() # determines how many iterations finished
    startcount = Counter() # determines how many iterations started
    threadcounter = Counter() 
    
    for i in range(niter):
      argsdict = {}
      argsdict['lock'] = lock
      argsdict['ncharA'] = ncharA
      argsdict['ncharB'] = ncharB
      argsdict['endcount'] = endcount
      argsdict['startcount'] = startcount
      argsdict['niter'] = niter
      argsdict['teststats'] = teststats
      argsdict['model'] = model
      argsdict['debug'] = debug
  
      # comment out this block to resample from concatenated set
      if i < niter/2: 
        argsdict['align'] = alignA
      else:
        argsdict['align'] = alignB
  
      # uncomment this line to resample from concatenated set
      #argsdict['align'] = catseq((alignA, alignB))
        
      q.put(argsdict)
    
    threads = []
    for i in range(nprocs):
      t = PhyloThread(q, lock, target = mlboot, name = 'thread %d' % i,
                        counter = threadcounter)
      threads.append(t)
      t.start()
            
      # wait until new thread really starts before starting next thread
      lock.acquire()
      lock.release()
      
    for t in threads:
      t.join()
    
    return teststats

def preparesets(seqfiles):
  """Concatenate all pairs of files in seqfiles, and return the new file names.

  Copy all files in seqfiles to consecutively numbered files with names of the
  form set*.seq, where * is a 3 digit number, and prepare pairwise concetenated
  sequence files with names of the form set*$.seq, where * is the 3 digit number
  assigned to the first file and $ is the 3 digit number assigned to the second.
  Return a list of the new files.
  
  Parameters:
    seqfiles -- a list of names of files to be copied/concatenated
    
  Return:
    List of file names for new concatenated data sets.
    
  """
  sets = []
  
  for i in range(len(seqfiles)):
    print 'Reading %s...' % seqfiles[i]
    writephy(getseqs(seqfiles[i]), 'set%03d.seq' % i)
    sets.append('set%03d.seq' % i)
          
  nsets = len(seqfiles)
    
  for i in range(nsets):
    for j in range(i + 1, nsets):
      aligns = getaligns((sets[i], sets[j]))
      unionI, unionJ = getsequnion(aligns[0], aligns[1])
      if len(unionI.seqlist) < 4:
        raise ConcaterpillarError('Only %d taxa shared between alignments %s and %s (must be at least 4)' \
                                  % (len(unionI.seqlist), seqfiles[i], seqfiles[j])) 
      writephy(catseq(aligns), 'set%03d%03d.seq' 
               % (i, j))
      sets.append('set%03d%03d.seq' % (i, j))
            
  return sets



def recover(model = None, nprocs = 1, rebuild = False, debug=False):
  """Use the file results.ccp to figure out where previous run was interrupted.
  
  Parameters:
    nprocs -- number of CPUs to use in recovery process (for tree-building etc.)
    
  Return:
    Current list of concatenated sets (list of lists).
  
  """  
  # open results file
  print "Attempting to retrieve current analysis state from results.ccp..."
  results = open('results.ccp', 'r')
    
  # retrieve innitial sequences
  sequences = getseqlist()
  
  setnums = len(sequences)
  catsets = []
  inverse = {}
  
  for i in xrange(len(sequences)):
    catsets.append([sequences[i]])
    inverse[sequences[i]] = i
    
  
  # parse concatenation steps from results file  
  concatRE = re.compile(r'Concatenating[\s\S]*(\[[\s\S]*\])[\s\S]*(\[[\s\S]*\])', \
                        re.I)
  
  # redo the concatenation steps
  index1 = -1
  for line in results:
    matchline = concatRE.match(line)
    if matchline:
      setnums -= 1
      group1, group2 = [i for i in matchline.groups()]
      species1 = group1[1:-2].split(',')[0].strip().strip('\'') 
      species2 = group2[1:-2].split(',')[0].strip().strip('\'') 
      index1 = inverse[species1]
      index2 = inverse[species2]
      catsets[index1].extend(catsets[index2])
      del inverse[species2]
      del catsets[index2]
      for key in inverse.keys():
        if inverse[key] > index2:
          inverse[key] -= 1
          
  # catsets is now recovered, and setnums tells you how many set numbers there 
  # should be i.e. if setnums is 25, then set024.seq is the last index that 
  # should be there (0 - 24)
  
  # check if there was an error while renaming files
  setlist = glob.glob('set%03d*' % setnums)
  if len(setlist) > 0:
    sys.stderr.write('Sorry, concaterpillar must have crashed unexpectedly, ')
    sys.stderr.write('use your last savepoint if possible.\n')
    shutdown(0)
  
  # congruence test hasn't started yet
  if index1 < 0:
    return catsets
  
  # if you don't want to recover anything
  if not rebuild:
    return catsets
  
  # find the sets that need to be reanalyzed
  setno1 = index1
  nsets = len(catsets)
  newsets = []
  for setid in range(nsets):
    if setid < setno1:
      setname = 'set%03d%03d.seq' % (setid, setno1)
    elif setid > setno1:
      setname = 'set%03d%03d.seq' % (setno1, setid)
    else:
      continue
    
    newsets.append(setname)
  
  # analyze data sets from last round
  if newsets:
    buildmultitrees(newsets, model, nprocs, debug=debug)
    runmultipuzzle(model, setnames = newsets, nprocs = nprocs, debug=debug)
  
  return catsets



def removesets(setid):
  """Remove all concatenated set files with setid as one of the set ID #s.

  The alignment files as well as any files produced from the alignment files
  (tree files, alpha files, etc.) are deleted.

  Parameters:
    setid -- identification number (int) for concatenated sets to be deleted.

  """
  
  setlist = []
  setRE1 = re.compile(r'set\d\d\d%03d' % setid)
  setRE2 = re.compile(r'set%03d\d\d\d' % setid)

  for f in os.listdir('.'):
    if setRE1.match(f) or setRE2.match(f):
      setlist.append(f)

  setlist.sort()
  
  removefiles(setlist)

def renamesets(setid):
  """Rename sets with IDs greater than setid.

  Remove sets that contain the set with id setid (sets that SHOULD be 
  useless at this point), and rename sets with IDs greater that setid (i.e.,
  if setid = 3, set003 sets should be removed, and set004 becomes set003, set005
  becomes set004, etc.)  In other words, shift set numbers down by 1.

  Parameters:
    setid -- ID number of the set above which sets are renamed.
    
  """

  namesdict = {}
  
  # concatRE will be something like "set000123"
  concatRE = re.compile(r'set(\d\d\d)(\d\d\d)')
  simpleRE = re.compile(r'set(\d\d\d)')
  
  for f in glob.glob('set*'):
    result = concatRE.match(f)
    
    if result: # set is concatenated pair (two ID numbers)
      num1, num2 = [int(i) for i in result.groups()]
    
      # don't need to rename sets that fall below setid in both cases
      if (num1 < setid) and (num2 < setid):
        continue
    
      if num1 > setid:
        num1 -= 1
      if num2 > setid:
        num2 -= 1
        
      suffix = f.split('.')[-1]
      
      namesdict[f] = 'set%03d%03d.%s' % (num1, num2, suffix)

    # set has a single ID number
    else: 
      result = simpleRE.match(f)
      try:
        num = int(result.group(1))
      except (ValueError, IndexError):
        raise ConcaterpillarError('Filename %s is not a legitamate set name!' \
                                  % f)
      except AttributeError:
        continue
      
      if num < setid:
        continue
      suffix = f.split('.')[-1]
      namesdict[f] = 'set%03d.%s' % (num - 1, suffix)
      
  renamefiles(namesdict)
        
def runpuzzlebrlen(argsdict):
  """Run RAxML to get likelihoods required for branch-length test.
  
  The name of this function is left over from a previous incarnation of 
  Concaterpillar, where likelihoods were evaluated with TREE-PUZZLE.
  
  Parameters:
    argsdict -- dictionary containing the following items:
      'tree': Newick-format tree string (with branch lengths, whatever)
      'seqfile': name of seqfile
      'gammacats': integer, number of categories for discrete Gamma distribution
                   (should probably just be 4)
      'model': substitution model to use.  Must be 'DAYHOFF', 'DCMUT', 'JTT', 
               'MTREV', 'WAG', 'RTREV', 'CPREV', 'VT', 'BLOSUM62', 'MTMAM', 'GTR'
  
  """
  
  tree = argsdict['tree']
  seqfile = argsdict['seqfile']
  gammacats = str(argsdict['gammacats'])
  model = argsdict['model']
  suffix = seqfile[:-4]
  
  if gammacats == 4 and prefix + '.lk' in os.listdir('.'):
    print 'Already done RAxML analysis of %s.' % seqfile
    return
  taxlst = getseqs(seqfile).taxlist
  
  # number of branch lengths + 1
  np = len(taxlst) * 2 - 2 
  paramfile = open('%s.prm' % seqfile[:-4], 'w')
  paramfile.write('%d\n' % np) 
  paramfile.close()
      
  treefile = open(suffix + '.tre', 'w')
  treefile.write('%s\n' % unroot(prunetree(tree, taxlst)))
  treefile.close()
  
  raxmlpath = globals()['raxmlpath']
      
  if model == 'GTR':
    binstring = '%s -n %s -m GTRGAMMA -s %s -f e -t %s.tre'\
                % (raxmlpath, suffix, seqfile, suffix)
  else:
    binstring = '%s -n %s -m PROTGAMMA%s -s %s -f e -t %s.tre'\
                % (raxmlpath, suffix, model, seqfile, suffix)
  try:
    runbin(binstring, (), argsdict['debug'])
  except KeyError:
    runbin(binstring, ())
  print 'Finished RAxML analysis of %s.' % seqfile
  sys.stdout.flush()
  
  # extract & save likelihood from puzzle outfile
  likefile = open('%s.lnl' % suffix, 'w')
  likefile.write('%f\n' % getlikelihood('RAxML_info.%s' % suffix))
  likefile.close()
  
  # delete RAxML outfiles
  for ext in ('result', 'log', 'info'):
    os.remove('RAxML_%s.%s' % (ext, suffix))
  try:
    os.remove('%s.reduced' % seqfile)
  except OSError:
    pass
    
def createsavepoint(debug=False):
  """Create a backup archive.
  
  Backup file will contain all set files and results.ccp.  File is named 
  savepoint.tbz If the program crashes at an irrecoverable point, this file can be
  manually extracted and used as a restart point.
  
  """

  os.mkdir('savepoint')
  for setfile in glob.glob('set*'):
    shutil.copy(setfile, 'savepoint')
  shutil.copy('results.ccp', 'savepoint')
  
  # create temporary savepoint
  args = ('tar', '-cjf', 'savetemp.tbz', 'savepoint')
  if debug:
    sys.stderr.write('Running command: %s\n' % ' '.join(args))
  process = subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
  
  out, err = process.communicate()
  if err:
    if debug:
      sys.stderr.write(err)
      
  if out:
    if debug:
      sys.stderr.write(out)
  
  shutil.rmtree('savepoint')
  if os.access('savepoint.tbz', os.F_OK):
    os.remove('savepoint.tbz')
  os.rename('savetemp.tbz', 'savepoint.tgz')
  sys.stdout.write('Savepoint file savepoint.tbz written.\n')
  sys.stdout.flush()


def toptest(catsets, likes, cutoff, model, nprocs = 1, saverounds = None, debug=False):
  """Execute the topological congruence test.

  Parameters:
    catsets -- list of lists containing sets that have already been concatenated
    likes -- matrix of pairs of likelihoods (as described in the documentation
             for buildmatrix, above)
    cutoff -- uncorrected alpha level for the test 
    
  """
    
  matrix = buildmatrix(likes)
  blacklisted = []
  for row in matrix:
    blacklisted.append([])
    for val in row:
      blacklisted[-1].append(False)
  outfile = open('results.ccp', 'a')
  
  nsets = 0
  for cs in catsets:
    nsets += len(cs)
    
  # current level is number of concatenations that have alredy happened (single-
  # gene sets minus remaining sets plus 1... note incrementation in loop below)
  level = nsets - len(catsets)
  
  iterationcounter = 0
  while len(matrix) > 1 and not allblacklisted(blacklisted):#0:
    
    level += 1
    row, col = getmin(matrix, blacklisted)
    
    setno1 = col
    setno2 = row
    
    if matrix[row][col] == 0:
      sys.stdout.write("Statistic is 0, p-value will be 1, skipping bootstraps.\n")
      sys.stdout.flush()
      pval = 1
      rawp = 1
    
    else:
      
      sys.stdout.write('\nStarting bootstrap replicates:\n')
      sys.stdout.flush()
      niter = 100
      
      results = multinonparboot(setno1, setno2, niter, model, nprocs, debug=debug)  
      sys.stdout.write('\n')
      sys.stdout.flush()
      results.sort()
      
      # Uncomment to record bootstrap test stats (and real test stat) to a file
      # bootfile = open('bootvals', 'a')
      # bootfile.write('%f ' % matrix[row][col])
      # for r in results:
      #   bootfile.write('%f ' % r)
      # bootfile.write('\n')

      rawp = getrawpval(matrix[row][col], results)
      try:
        pval = getpval(matrix[row][col], results)
      except ValueError:
        print 'Poor fit of Weibull distribution, using raw p-value (%.2f)' % rawp
        pval = rawp
      else:  
        print 'Raw p-value: %.2f' % rawp
        print 'Weibull-smoothed p-value: %f' % pval
    
    # Note: the original correction method was to "predict" the number of levels
    # using the uncorrected alpha, then to correct based on this value.  But
    # since the correction is used WHENEVER the null is rejected with the
    # uncorrected alpha, it's simpler to just always correct based on the
    # minimum number of levels... i.e., the current level.  This works because
    # the corrected alpha is always less than or equal to the user-defined alpha
    
    alpha = correctalpha(cutoff, level)

    # ... but only tell the user we're correcting when we NEED to correct
    if pval < cutoff:
      print 'Applying alpha-level correction: new alpha = %.6f' % alpha
      outfile.write('Applying alpha-level correction: new alpha = %.6f\n' % alpha)
      sys.stdout.flush()
      outfile.flush()
      
    # if null is rejected
    if pval < alpha:
      print 'Failed to concatenate %s and %s (p = %.6f).\n' % \
              (str(catsets[setno1]), str(catsets[setno2]), pval)
      outfile.write('\nFailed to concatenate %s and %s (p = %.6f).\n\n' % \
                    (str(catsets[setno1]), str(catsets[setno2]), pval))
      outfile.flush()
      for i in xrange(len(blacklisted)):
        if i < setno1:
          blacklisted[setno1][i] = True
        elif i > setno1:
          blacklisted[i][setno1] = True
        if i < setno2:
          blacklisted[setno2][i] = True
        elif i > setno2:
          blacklisted[i][setno2] = True      
    
    # otherwise...
    else:
      
      for i in xrange(setno2 + 1, len(blacklisted)):
        del blacklisted[i][setno2]
      del blacklisted[setno2]
      
      print 'Concatenating %s and %s (p = %.6f).\n' % \
            (str(catsets[setno1]), str(catsets[setno2]), pval)
      outfile.write('Concatenating %s and %s (p = %.6f).\n' % \
                    (str(catsets[setno1]), str(catsets[setno2]), pval))
      outfile.flush()
          
      namesdict = {}
      removefiles(glob.glob('set%03d.*' % setno1))

      for suffix in ('seq', 'tre'):
        namesdict['set%03d%03d.%s' % (setno1, setno2, suffix)] = \
                'set%03d.%s' % (setno1, suffix)
      renamefiles(namesdict)
      
      # only necessary if setno2 is set with highest ID (should always be true)
      removefiles(glob.glob('set%03d.*' % setno2))
      removesets(setno1)
      removesets(setno2)
      renamesets(setno2)
      
      catsets[setno1].extend(catsets[setno2])
      del catsets[setno2]

      # concatenate necessary data sets
      nsets = len(catsets)
      newsets = []
      for setid in range(nsets):
        if setid < setno1:
          setname = 'set%03d%03d.seq' % (setid, setno1)
        elif setid > setno1:
          setname = 'set%03d%03d.seq' % (setno1, setid)
        else:
          continue
        writephy(catseq(getaligns(('set%03d.seq' % setid, 
                                  'set%03d.seq' % setno1))), setname)
        newsets.append(setname)
      
      # analyze new data sets
      if newsets:
        buildmultitrees(newsets, model, nprocs, debug=debug)
        runmultipuzzle(model, setnames = newsets, nprocs = nprocs, debug=debug)
      
        likes = getalllikes(len(catsets))
        
        # make a savepoint if nessecary
        if saverounds:
          if iterationcounter % saverounds == 0 and\
            not iterationcounter == 0:
            createsavepoint()
          iterationcounter += 1
        
        matrix = buildmatrix(likes)
      
      # if there are no new sets, I think we may still need to continue 
      # (if there are still sets that aren't blacklisted)
      #else:
      #  break

  
  print 'Finished concatenating sets.  The following files contain the',
  print 'final sets:'
  outfile.write('Finished concatenating sets.  ')
  outfile.write('The following files contain the final sets:\n')

  for i in range(len(catsets)):
    print 'set%03d.seq (%d genes): %s' % (i, len(catsets[i]), str(catsets[i]))
    outfile.write('set%03d.seq (%d genes): %s\n' \
                  % (i, len(catsets[i]), str(catsets[i])))

  outfile.close()
  
def allblacklisted(matrix):
  allBL = True
  
  for row in matrix:
    for item in row:
      if not item:
        allBL = False
        break
    if not allBL:
      break
  
  return allBL

def usage(exitcode):
  """Print a usage message and exit with the given exit code."""
  
  from concaterpillar import __doc__
  
  sys.stderr.write(__doc__)
  shutdown(exitcode)


class ConcaterpillarError(Exception):
  """
  Exception class designed for concaterpillar-specific errors.
  
  Extends Exception.
  
  """

  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

    
class Counter:
  """
  Implementation of an integer object that can be passed as a reference.
    
  """

  def __init__(self, initval = 0):
    self.count = initval
    
  def increment(self):
    self.count += 1
  
  def decrement(self):
    self.count -= 1
    
  def __int__(self):
    return self.count
    
  def __repr__(self):
    return str(self.count)
  
  def __str__(self):
    return repr(self)

class PhyloThread(threading.Thread):  
  """
  Thread that runs an appropriate target phylogeny programme.
  
  Extends threading.Thread
  
  """
  
  def __init__(self, jobqueue, lock, target, counter, name = None):
    """Create a PhyloThread object.
    
    Parameters:
      jobqueue -- queue from which job objects are drawn by all active threads
      lock -- shared lock object to prevent deadlocking of threads
      target -- phylogeny programme to run
      counter -- Counter object to count active threads
      name -- name assigned to thread (doesn't really DO anything)
    
    """

    threading.Thread.__init__(self, name = name)
    self.name = name
    self.jobqueue = jobqueue
    self.lock = lock
    self.target = target
    self.counter = counter # counts active threads
    
    self.lock.acquire()
    self.counter.increment()
    self.lock.release()
    
  def run(self):

    while True:
      self.lock.acquire()
      
      if self.jobqueue.empty():
        self.lock.release()
        break
      
      arg = self.jobqueue.get()
      self.lock.release()
      
      self.target(arg)
 
    self.lock.acquire() 
    self.counter.decrement()
    self.lock.release()      

if not mpipresent:
  def runmultipuzzle(model, nsets = None, setnames = None, nprocs = 1, debug=False):
    """Launch multiple threads to run RAxML for toptest.

    The name of this function is left over from a previous version of 
    Concaterpillar that used TREE-PUZZLE to calculate likelihoods
    
    This function should be called with either nsets or setnames, but not both.
    If setnames is specified, the numbers of the sets in the setnames list will
    be extracted and passed to runpuzzle.  If nsets is specified, the appropriate
    set numbers will be generated.
  
    Parameters:
      nsets -- integer specifying the number of (single) sets to be reanalysed.
      setnames -- names of sets to be reanalysed.
      nprocs -- number of processes to thread off... should be number of 
                processors to use on a shared memory machine
  
    """  
        
    print 'Starting RAxML reanalysis of trees.'
    
    q = Queue.Queue()
    lock = threading.Lock()
    threadcounter = Counter()
    
    try:
      if setnames:
        for name in setnames:
          i = int(name[3:6])
          j = int(name[6:9])
          q.put((i, j, model, debug))
          q.put((j, i, model, debug))
          
      elif nsets:
        for i in range(nsets):
          for j in range(nsets):
            if i != j:
              q.put((i, j, model, debug))
      
      else:
        sys.stderr.write('runmultipuzzle was not called properly.\n%s' % \
                         runmultipuzzle.__doc__)
        sys.exit(1)
              
    except IOError: 
      sys.stderr.write('Invalid starting conditions for tree reanalysis.\n')
      usage(1)
            
    for i in range(nprocs):
      t = PhyloThread(q, lock, runpuzzle, name = 'thread %d' % i, 
                      counter = threadcounter)
      t.start()
      
      # wait until newly-started thread really starts before starting next thread
      lock.acquire()
      lock.release()
      
    while int(threadcounter) > 0:
      time.sleep(10) 
      
    print 'Finished reanalysis.'
    sys.stdout.flush()

def runpuzzle(seqids):
  """Run RAxML on a set with trees from single and concatenated sets.
  
  The name of this function is left over from a previous version of 
  Concaterpillar that used TREE-PUZZLE to calculate likelihoods
  
  Taxa may be pruned from trees if the two data sets don't share the same taxa.
  The tree inferred from the data set with seqids[0] is only reanalysed if taxa
  were pruned, or if the unpruned tree wasn't already reanalysed.
  
  Parameters:
    seqids -- tuple containing two set ID numbers: the first is the set for
              which likelihoods will be determined, the second is the set with
              which set seqids[0] was concatenated to get the second tree to be
              evaluated
              
  Return:
    Log-likelihoods for trees evaluated with the alignment in question.
    
  """

  raxmlpath = globals()['raxmlpath']
  
  if len(seqids) == 4:
    seqid1, seqid2, model, debug = seqids
  else:
    seqid1, seqid2, model = seqids
    debug = False

  tf1 = open('set%03d.tre' % seqid1)
  tree1 = string.strip(tf1.readline())
  tf1.close()
   
  tree1 = root(cleantree(tree1))
    
  if seqid1 < seqid2:
    tf2 = open('set%03d%03d.tre' % (seqid1, seqid2))
  else:
    tf2 = open('set%03d%03d.tre' % (seqid2, seqid1))
    
  tree2 = string.strip(tf2.readline())
  tf2.close()
  tree2 = root(cleantree(tree2))
              
  align1 = getseqs('set%03d.seq' % seqid1)
  align2 = getseqs('set%03d.seq' % seqid2)
  
  taxa1 = align1.taxlist[:]
  taxa1.sort()
  taxa2 = align2.taxlist[:]
  taxa2.sort()
  
  foundlike = False

  # If the two alignments share the same taxa, and tree1 has already been 
  # reevaluated, take note.
  if taxa1 == taxa2:
    sametaxa = True
    try:
      likefile = open('set%03d.lnl' % seqid1)
      while True:
        try:
          lnL1 = float(likefile.readline())
          likefile.close()
          foundlike = True
          break
        except ValueError:
          pass
      
    except IOError:
      pass
    
  else:
    sametaxa = False
    align1, align2 = getsequnion(align1, align2)

  intreename = 'set%03d%03d.raxtre' % (seqid1, seqid2)
  treefile = open(intreename, 'w')

  if foundlike:
    treefile.write('%s\n' % tree2)

  # No likelihood yet calculated for this tree with these taxa
  else:
    treefile.write('%s\n%s\n' % (prunetree(tree1, align1.taxlist), 
                                 prunetree(tree2, align1.taxlist)))
  treefile.close()
  
  treefile = open(intreename, 'r')
  infilename = 'set%03d%03d.raxin' % (seqid1, seqid2)
  writephy(align1, infilename)
  
  runnum = 'A'
  treefile = open(intreename, 'r')
  likes = []
  for line in treefile:
    if not line.strip() == "":
      suffix = 'set%03d%03d%s' % (seqid1, seqid2, runnum)
      temptreename = 'set%03d%03d%s.raxtre' % (seqid1, seqid2, runnum)
      temptree = open(temptreename, 'w')
      temptree.write(line)
      temptree.close()
      
      if model == 'GTR':
        binstring = '%s -n %s -m GTRGAMMA -s %s -f e -t %s' \
                    % (raxmlpath, suffix, infilename, temptreename)
      else:
        binstring = '%s -n %s -m PROTGAMMA%s -s %s -f e -t %s' \
                    % (raxmlpath, suffix, model, infilename, temptreename)
      
      runbin(binstring, (), debug=debug)
      runnum = 'B'
      
      like = getlikelihood('RAxML_info.%s' % suffix)
      likes.append(like)
    
  treefile.close()
  suffix = 'set%03d%03d' % (seqid1, seqid2)
  freefiles = []
  
  freefiles.extend(glob.glob('RAxML*.%s*' % suffix))
  freefiles.extend(glob.glob('%s*.rax*' % suffix))
  
  removefiles(freefiles)
  try:
    os.remove('%s.reduced' % infilename)
  except OSError:
    pass

  if foundlike:
    likes[:0] = [lnL1]
  elif sametaxa:
    likefile = open('set%03d.lnl' % seqid1, 'w')
    likefile.write('%f\n' % likes[0])
    likefile.close()
    
  likefile = open('set%03d%03d.lnl' % (seqid1, seqid2), 'w')
  likefile.write('%f\n%f\n' % tuple(likes))
  
  return likes

  
def shutdown(exitstatus = 0):
  """ Wrapper function for sys.exit, to use instead
  
  Shuts down MPI if it is running, then exits with sys.exit.
  
  Parameters:
    exitstatus -- exit status to send to OS.
    
  """
  
  if globals()['mpipresent']:
    mpi.bcast(None)
    mpi.finalize()
    
  sys.exit(exitstatus)


###########################################################
#                                                         #
#         This whole section contains MPI code            #
#                                                         #
###########################################################

if mpipresent:
  
  def brlentest(catsets, setdir, cutoff, model, nprocs = 1, debug=False):
    """Execute the branch-length congruence test.

    Parameters:
      catsets -- list of topologically congruent sets (produced by toptest)
      setdir -- "name" for branch-length test run (e.g. BLtest-set000/)
      cutoff -- uncorrected alpha-level for the test
      nprocs -- number of processes to thread off... should be number of 
                processors to use on a shared memory machine
  
    """
        
    resfile = open('%s.ccp' % setdir[:-1].lower(), 'w')
    gammacats = 4
    lastpval = 0
    q = Queue.Queue()
    
    treefile = open('global.tre')
    tree = treefile.readline().strip()
    treefile.close()
    
    nsets = 0
    for cs in catsets:
      nsets += len(cs)
      
    # current level is number of concatenations that have alredy happened (single-
    # gene sets minus remaining sets plus 1... note incrementation in loop below)
    level = nsets - len(catsets)
  
    while len(catsets) > 1:
      
      # put each process in the correct directory
      mpi.bcast((1, os.chdir))
      for i in xrange(1, mpi.size):
        mpi.recv(i)
        mpi.send(os.getcwd(), i)
        mpi.recv(i)
        mpi.send(None, i)
    
      level += 1
      dlist = os.listdir('.')
  
      for i in range(len(catsets)):
        if 'set%03d.lnl' % i not in dlist:
          q.put({'tree': tree, 'seqfile': 'set%03d.seq' % i, \
               'gammacats': gammacats, 'model':model, 'debug':debug})
  
        for j in range(i + 1, len(catsets)):
          catsetname = 'set%03d%03d.seq' % (i, j)
  
          # if concatenated file doesn't exist, create it
          if catsetname not in dlist: 
            writephy(catseq(getaligns(('set%03d.seq' % i, 'set%03d.seq' % j))), 
                     catsetname)
          
          if 'set%03d%03d.lnl' % (i, j) not in dlist:
            q.put({'tree': tree, 'seqfile': catsetname, 'gammacats': gammacats, \
                   'model':model, 'debug':debug})
      
      
      mpi.bcast((1, runpuzzlebrlen))
      if mpi.size > globals()['mpisize']:
        coordinatewait(q)
      else:
        t = threading.Thread(target=extrathread, args=(q, runpuzzlebrlen))
        t.start()
        coordinate(q)
        t.join()
            
      # dlnL (lnL(A) + lnL(B) - lnL(AB))

      pvalmat = []
      
      for i in range(len(catsets)):
        pvalmat.append([])
        for j in range(i + 1, len(catsets)):
          result = calcBLratio(i, j)
          pvalmat[-1].append(1 - pchisq(2 * result[0], result[1]))
           
      row, col = getmax(pvalmat)
      setA = row
      setB = col + row + 1
      
      # p-val for chi-square test (note that lnL is multiplied by 2)
      #pval = 1 - pchisq(2 * matrix[row][col], dfmat[row][col])      
      pval = pvalmat[row][col]      
      
      # Correct alpha based on current level... see argument in toptest
      alpha = correctalpha(cutoff, level)
      
      # Only tell the user we're correcting if we NEED to correct
      if pval < cutoff:
        print 'Applying alpha-level correction: new alpha = %.6f' % alpha
        resfile.write('Applying alpha-level correction: new alpha = %.6f\n' % alpha)
        sys.stdout.flush()
        resfile.flush()
      
      if pval < alpha:  
        print 'Failed to concatenate %s and %s (p = %.6f).\n' % \
              (str(catsets[setA]), str(catsets[setB]), pval)
        resfile.write('\nFailed to concatenate %s and %s (p = %.6f).\n\n' % \
                      (str(catsets[setA]), str(catsets[setB]), pval))
        resfile.flush()
        break
        
      else: # p-val of chi-square test >= corrected alpha
        print 'Concatenating %s and %s (p = %.6f).\n' % \
              (str(catsets[setA]), str(catsets[setB]), pval)
        resfile.write('Concatenating %s and %s (p = %.6f).\n' % \
                      (str(catsets[setA]), str(catsets[setB]), pval))
        resfile.flush()
        
        # rename and remove necessary sets
        filedict = {}
        for i in ('.seq', '.lnl', '.prm'):
          filedict['set%03d%03d%s' % (setA, setB, i)] = 'set%03d%s' % (setA, i)
        renamefiles(filedict)
  
        # necessary only if setB is highest-numbered set
        removefiles(glob.glob('set%03d*' % setB))
        removesets(setA)
        removesets(setB)
        renamesets(setB)
        
        catsets[setA].extend(catsets[setB])
        del catsets[setB]
        
    print 'Finished concatenating sets.  The following files contain ' + \
          'the final sets:'
    resfile.write('Finished concatenating sets.  The following files contain ' + \
                  'the final sets:\n')
  
    for i in range(len(catsets)):
      print '%sset%03d.seq (%d genes): %s' % (setdir, i, len(catsets[i]), 
                                                   str(catsets[i]))
      resfile.write('%sset%03d.seq (%d genes): %s\n' % (setdir, i, 
                                                         len(catsets[i]), 
                                                         str(catsets[i])))
    resfile.close()
    
    print ''
  
  
  def buildmultitrees(infiles, model, nprocs = 1, debug=False):
    """Run RAxML on a list of alignment filenames.
  
    Parameters:
      infiles -- list of file names on which to run RAxML
      nprocs -- number of processes to thread off... should be number of 
                processors to use on a shared memory machine
                
      model -- substitution model to use.  Must be 'DAYHOFF', 'DCMUT', 'JTT', 
               'MTREV', 'WAG', 'RTREV', 'CPREV', 'VT', 'BLOSUM62', 'MTMAM', 'GTR'

  
    """
 
    print 'I\'m about to build trees for %d data sets (using RAxML).' % (len(infiles))
    sys.stdout.flush()

    q = Queue.Queue()
    for i in infiles:
      cmds = {'model': model, 'seqfile': i, 'debug':debug}
      q.put(cmds)
    
    mpi.bcast((1, buildtrees))
    if mpi.size > globals()['mpisize']:
      coordinatewait(q)
    else:
      t = threading.Thread(target=extrathread, args=(q, buildtrees))
      t.start()
      coordinate(q)
      t.join()
  
  def mlboot(args):
    """Determine the topology test statistic by ML using a bootstrapped dataset.
    
    RAxML is used to both infer trees and re-evaluate likelihoods.
  
    Parameters:
      args -- a dictionary containing the following items:
        'lock': a threading.Lock object used to prevent deadlocking
        'startcount': a Counter used to determine iteration number
        'align': alignment from which sites are sampled
        'ncharA': number of sites to sample for the first alignment
        'ncharB': number of sites to sample for the second alignment
        'endcount': a Counter used to determine how many bootstrap replicates 
                    have finished
        'niter': indicates total number of replicates to perform
        'teststats': list to which the test statistic is appended
  
    """
    
    raxmlpath = globals()['raxmlpath']
    raxmlversion = globals()['raxmlversion']
    
    iterno = int(args['startcount'])
    
    # names of bootstrap files
    bootseqs = ['boot%03d%s.seq' % (iterno, c) for c in ('a', 'b', 'ab')]
    
    writephy(bootstrap(args['align'], args['ncharA']), bootseqs[0])
    writephy(bootstrap(args['align'], args['ncharB']), bootseqs[1])
    
    writephy(catseq(getaligns((bootseqs[0], bootseqs[1]))), bootseqs[2])
    
    for seq in bootseqs:
      suffix = seq[:-4]
      
      try:
        runraxml(raxmlpath, seq, suffix, model=args['model'], version=raxmlversion, debug=args['debug'])
      except KeyError:
        runraxml(raxmlpath, seq, suffix, model=args['model'], version=raxmlversion)

      os.rename('RAxML_result.%s' % suffix, '%s.tre' % suffix)
      for raxfile in glob.glob('RAxML_*.%s' % suffix):
        os.remove(raxfile)
      try: 
        os.remove('%s.reduced' % seq)
      except OSError:
        pass
      
    # likelihood re-evaluation

    likes = []    
    trees = []
      
    for seq in bootseqs:
      treefile = open(seq[:-4] + '.tre')
      trees.append(string.strip(treefile.readline()))
      treefile.close()
  
    for i in (0, 1):
      intree = '%s.intree' % (bootseqs[i][:-4])
      treefile = open(intree, 'w')
      treefile.write('%s\n%s\n' % (trees[i], trees[2]))
      treefile.close()
      # split tree for RAxML and make 2 runs
      treefile = open(intree, 'r')
      for letter in ('A', 'B'):
        suffix = '%s%s' % (bootseqs[i][:-4], letter)
        temptreename = '%s.temptree' % suffix
        temptree = open(temptreename, 'w')
        treeline = treefile.readline().strip()
        temptree.write(treeline)
        temptree.close()
        
        if args['model'] == 'GTR':
          binstring = '%s -n %s -m GTRGAMMA -s %s -f e -t %s' \
                      % (raxmlpath, suffix, bootseqs[i], temptreename)
        else:
          binstring = '%s -n %s -m PROTGAMMA%s -s %s -f e -t %s' \
                      % (raxmlpath, suffix, args['model'], bootseqs[i], temptreename)
        try:
          runbin(binstring, (), debug=args['debug'])
        except KeyError:
          runbin(binstring, ())
        likes.append(getlikelihood('RAxML_info.%s' % suffix))
  
  
    for f in glob.glob('*boot%03d*' % iterno):
      os.remove(f)
  
    return likes

  
  def multinonparboot(setA, setB, niter, model, nprocs, debug=False):
    """Perform several ML-based non-parametric bootstrap replicates.  
  
    Parameters:
      setA -- number of first set
      setB -- number of second set
      niter -- number of times to repartition data set and get pseudo test stat
      nprocs -- number of processes to thread off... should be number of 
                processors to use on a shared memory machine
              
    Return:
      List of pseudo test statistics

    """
    
    alignA = getseqs('set%03d.seq' % setA)
    alignB = getseqs('set%03d.seq' % setB)
    
    # want same # of taxa in bootstraps as in real dataset
    alignA, alignB = getsequnion(alignA, alignB)
    
    ncharA = alignA.length
    ncharB = alignB.length
    
    q = Queue.Queue()
    endqueue = Queue.Queue()
    resultqueue = Queue.Queue()
    
    for i in range(niter):
      argsdict = {}
      argsdict['ncharA'] = ncharA
      argsdict['ncharB'] = ncharB
      argsdict['startcount'] = i
      argsdict['niter'] = niter
      argsdict['model'] = model
      argsdict['debug'] = debug
  
      # comment out this block to resample from concatenated set
      if i < niter/2: 
        argsdict['align'] = alignA
      else:
        argsdict['align'] = alignB
  
      # uncomment this line to resample from concatenated set
      #argsdict['align'] = catseq((alignA, alignB))
        
      q.put(argsdict)
      endqueue.put(i)

    mpi.bcast((2, mlboot))
    if mpi.size > globals()['mpisize']:
      coordinatelistwait(q, endqueue, resultqueue, niter)
    else:
      t = threading.Thread(target=listthread, args=(q, mlboot, endqueue, resultqueue))
      t.start()
      coordinatelist(q, endqueue, resultqueue, niter)
      t.join()
      
    teststats = []
    for i in xrange(resultqueue.qsize()):
      teststats.append(resultqueue.get())

    return teststats
  
  
  def runmultipuzzle(model, nsets = None, setnames = None, nprocs = 1, debug=False):
    """Launch multiple threads to run RAxML for toptest.
    
    The name of this function is left over from a previous version of 
    Concaterpillar that used TREE-PUZZLE to calculate likelihoods
  
    This function should be called with either nsets or setnames, but not both.
    If setnames is specified, the numbers of the sets in the setnames list will
    be extracted and passed to runpuzzle.  If nsets is specified, the appropriate
    set numbers will be generated.
  
    Parameters:
      nsets -- integer specifying the number of (single) sets to be reanalysed.
      setnames -- names of sets to be reanalysed.
      nprocs -- number of processes to thread off... should be number of 
                processors to use on a shared memory machine
  
  """  
    
    print 'Starting RAxML reanalysis of trees.'
    
    q = Queue.Queue()
    
    try:
      if setnames:
        for name in setnames:
          i = int(name[3:6])
          j = int(name[6:9])
          q.put((i, j, model, debug))
          q.put((j, i, model, debug))
          
      elif nsets:
        for i in range(nsets):
          for j in range(nsets):
            if i != j:
              q.put((i, j, model, debug))
      
      else:
        sys.stderr.write('Uh oh! Something\'s wrong in runmultipuzzle!')
              
    except IOError: 
      sys.stderr.write('Invalid starting conditions for tree reanalysis.\n')
      usage(1)
    
    mpi.bcast((1, runpuzzle))
    if mpi.size > globals()['mpisize']:
      coordinatewait(q)
    else:
      t = threading.Thread(target=extrathread, args=(q, runpuzzle))
      t.start()
      coordinate(q)
      t.join()  
 
    print 'Finished reanalysis.'  
  
  def nofeedbackjobs(action):
    """Assign a job to an MPI thread.
    
    The threads will loop, waiting for a message, until they get a None message.
    
    Parameters:
      action -- job to be assigned

    """
    
    while True:     
      mpi.send("done", 0)
      arg, status = mpi.recv(0)
      if arg == None:
        return
      else:
        action(arg)
  
  
  def listjobs(action):
    """Assign a job to an MPI thread, passing a list back to main thread. 
    
    The threads will loop, waiting for a message, until they get a None message.

    Parameters:
      action -- job to be assigned

    """
    
    li = []
    while True:
      mpi.send(li, 0)
      arg, status = mpi.recv(0)
      #print mpi.rank, "performing", arg
      if arg == None:
        return
      else:
        li = action(arg)
  
  
  def extrathread(jobqueue, target):
    """Take care of the host's part of the job, passing a list back to main thread.
    
    Will execute the target thread while jobqueue is not empty.
    
    Parameters:
      jobqueue -- a Queue object holding jobs
      target -- function to be executed with job from jobqueue
    
    """
    
    while True:
      try:
        arg = jobqueue.get_nowait()
        target(arg)
      except Queue.Empty:
        return
  
  
  def listthread(jobqueue, target, endqueue, resultqueue):
    """Take care of the host's part of the job, passing a list to main thread.
    
    Will execute the target thread while jobqueue is not empty.
    
    Parameters:
      jobqueue -- a Queue object holding jobs
      target -- function to be executed with job from jobqueue
      endqueue -- a Queue object with intergers, pulled when job is finished
      resultqueue -- a Queue object that will store the results
      
    """

    while True:
      try:
        arg = jobqueue.get_nowait()
        li = target(arg)
        resultqueue.put(max(li[:2]) - li[1] + \
                         max(li[2:]) - li[3])
        endcount = int(endqueue.get_nowait())
        if endcount % 10 == 0:
          sys.stdout.write('[.')
        elif endcount % 10 == 9:
          sys.stdout.write('.] %d/%d\n' % (endcount + 1, arg['niter']))
        else:
          sys.stdout.write('.')    
        sys.stdout.flush()
      except Queue.Empty:
        return
   
      
  def getmpisize():
    """Determine how many processors to use.
    
    Tell processors of too high rank that it's done.

    Return: 
      Number of processors to use (int).
    
    """
    
    
    # advanced option not used
    if not os.access('_mpiusagemap_', os.R_OK):
      return mpi.size
    
    usagefile = open('_mpiusagemap_', 'r')
    usagedict = {}
    usageRE = re.compile(r'(\d+)[ \t]+(\d+)')
    for line in usagefile:
      usagematch = usageRE.match(line)
      if usagematch:
        usagedict[int(usagematch.group(1))] = int(usagematch.group(2))
    usagefile.close()
    
    # print usagedict
    
    # got my usage dict
    miss = False
    i = 0
    x = 0
    while True:
      if os.access('set%03d.seq' % x, os.F_OK):
        i = x
      else:
        if miss:
          break
        miss = True
      x += 1
    
    keys = usagedict.keys()
    keys.sort()
    keys.reverse()
    mpisize = mpi.size
    newsize = mpisize
    for key in keys:
      if key >= i:
        if usagedict[key] < mpisize and usagedict[key] > 0:
          newsize = usagedict[key]
    
    # tell processes over desired size to end
    for x in xrange(newsize, mpisize):
      msg, status = mpi.recv(x)
      mpi.send(None, x)
    
    return newsize
      
  
  def coordinatewait(jobqueue):
    """Coordinate mpi threads by handing out jobs.
    
    Parameters:
      jobqueue -- a Queue object holding jobs
      
    """
    mpisize = getmpisize()
    done = mpisize - 1
    while done > 0:
      msg, status = mpi.recv()
      try:
        arg = jobqueue.get_nowait()
        mpi.send(arg, status.source)
        #print "gave", arg, "to", key
      except Queue.Empty:
        mpi.send(None, status.source)
        done = done - 1   
            
  
  
  def coordinatelistwait(jobqueue, endqueue, resultqueue, niter):
    """Coordinate multinonparboot using mpi by handing out jobs.
    
    Result is added to resultqueue.
    
    Also coordinates multipartition, if used.
    
    Parameters:
      jobqueue -- a Queue object holding jobs
      endqueue -- a Queue object with intergers, pulled when job is finished
      resultqueue -- a Queue object that will store the results
      niter -- number of bootstrap replicates
      
    """
    
    mpisize = getmpisize()
    
    done = mpisize - 1
    while done > 0:
      msg, status = mpi.recv()         
      if len(msg) > 0:
        resultqueue.put(max(msg[:2]) - msg[1] + \
                     max(msg[2:]) - msg[3])
        endcount = int(endqueue.get_nowait())
        if endcount % 10 == 0:
          sys.stdout.write('[.')
        elif endcount % 10 == 9:
          sys.stdout.write('.] %d/%d\n' % (endcount + 1, niter))
        else:
          sys.stdout.write('.')    
        sys.stdout.flush()
        
      try:
        arg = jobqueue.get_nowait()
        mpi.send(arg, status.source)
      except Queue.Empty:
        mpi.send(None, status.source)
        done = done - 1

  
  
  
  def coordinate(jobqueue):
    """Coordinate mpi threads by handing out jobs.
    
    Parameters:
      jobqueue -- a Queue object holding jobs
    
    """
    
    mpisize = getmpisize()
    
    done = mpisize - 1
    requestdict = {}
    for x in xrange(1, mpisize):
      requestdict[x] = mpi.irecv(x)
    while done > 0:
      for key in requestdict.keys():
        if requestdict[key]:
          try:
            arg = jobqueue.get_nowait()
            mpi.send(arg, key)
            requestdict[key] = mpi.irecv(key)
            #print "gave", arg, "to", key
          except Queue.Empty:
            mpi.send(None, key)
            done = done - 1
            del requestdict[key]
          
            
      time.sleep(3)
      
    
  def coordinatelist(jobqueue, endqueue, resultqueue, niter):
    """Coordinate multinonparboot using mpi by handing out jobs.
    
    Result is added to resultqueue.
    
    Also coordinates multipartition, if used.

    Parameters:
      jobqueue -- a Queue object holding jobs
      endqueue -- a Queue object with intergers, pulled when job is finished
      resultqueue -- a Queue object that will store the results
      niter -- number of bootstrap replicates
      
    """
    mpisize = getmpisize()
    
    done = mpisize - 1
    requestdict = {}
    for x in xrange(1, mpisize):
      requestdict[x] = mpi.irecv(x)
    while done > 0:
      for key in requestdict.keys():
        if requestdict[key]:
          message = requestdict[key].message
          
          # updates the list and counters
          if len(message) > 0:
            resultqueue.put(max(message[:2]) - message[1] + \
                         max(message[2:]) - message[3])
            endcount = int(endqueue.get_nowait())
            if endcount % 10 == 0:
              sys.stdout.write('[.')
            elif endcount % 10 == 9:
              sys.stdout.write('.] %d/%d\n' % (endcount + 1, niter))
            else:
              sys.stdout.write('.')    
            sys.stdout.flush()
            
          # fetch next job
          try:
            arg = jobqueue.get_nowait()
            mpi.send(arg, key)
            requestdict[key] = mpi.irecv(key)
          except Queue.Empty:
            mpi.send(None, key)
            done = done - 1
            del requestdict[key]

      time.sleep(3)
  
  
  def processdispatcher():
    """Dispatch mpi threads.
    
    Threads will wait for a bcast, and act upon the message.  Pass the function 
    you want to execute, None to quit
    e.g. mpi.bcast(myfunc) will make all threads execute myfunc

    """
    
    while True:
      try:
        action = mpi.bcast()
      except RuntimeError:
        return
      if action:
        if action[0] == 1:
          nofeedbackjobs(action[1])
        if action[0] == 2:
          listjobs(action[1])
      else:
        mpi.finalize()
  


