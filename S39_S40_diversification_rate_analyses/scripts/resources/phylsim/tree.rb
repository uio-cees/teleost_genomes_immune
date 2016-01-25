class Tree
  private_class_method :new
  attr_reader :branch, :species, :lambda, :mu, :treeOrigin, :present, :k, :rootSplit, :Np, :sumOfSpeciesDurations, :treeReconstructed, :branchRatesAssigned, :sequencesEvolved, :traitsEvolved, :alignmentPositionsAssigned, :fossilRecordAdded, :posteriorLambda, :posteriorMu, :posteriorTreeOrigin, :substitutionModel, :samplingScheme, :NpFull, :focusGroupNpFull, :focusGroup, :focusGroupAge, :ecmFileName, :branchRatesMean, :branchRatesStandardDeviation, :branchRatesAutoCorrelation
  alias :extantSpeciesNumber :Np
  def Tree.load(fileName = "bd.dmp", verbose = true)
    raise "No file name specified!" if fileName == nil
    raise "File #{fileName} could not be found!" unless File.exists?(fileName.to_s)
    new(lambda = nil, mu = nil, treeOrigin = nil, present = nil, k = 0, rootSplit = nil, np = nil, npEach = nil, checkProbabilities = nil, algorithm = nil, verbose, threads = nil, fileName, diversityFileName = nil, fileType = "dump")
  end
  def Tree.parse(fileName = "bd.nex", fileType = "nexus", diversityFileName = nil, treeNumber = 0, verbose = true)
    raise "No file name specified!" if fileName == nil
    raise "File #{fileName} could not be found!" unless File.exists?(fileName.to_s)
    new(lambda = nil, mu = nil, treeOrigin = nil, present = nil, k = 0, rootSplit = nil, np = nil, npEach = nil, checkProbabilities = nil, algorithm = nil, verbose, threads = nil, fileName, diversityFileName, fileType, treeNumber)
  end
  def Tree.generate(lambda = nil, mu = nil, treeOrigin = nil, present = 0, k = 0, rootSplit = false, np = [0,'inf'], npEach = [0,'inf'], checkProbabilities = false, algorithm = "forward", verbose = true, threads = 1)
    new(lambda, mu, treeOrigin, present, k, rootSplit, np, npEach, checkProbabilities, algorithm, verbose, threads, fileName = nil, fileType = nil)
  end
  def Tree.quietlyGenerate(lambda = nil, mu = nil, treeOrigin = nil, present = 0, k = 0, rootSplit = false, np = [0,'inf'], npEach = [0,'inf'], checkProbabilities = false, algorithm = "forward", verbose = false, threads = 1)
    new(lambda, mu, treeOrigin, present, k, rootSplit, np, npEach, checkProbabilities, algorithm, verbose, threads, fileName = nil, fileType = nil)
  end
  def initialize(lambda = nil, mu = nil, treeOrigin = nil, present = 0, k = 0, rootSplit = false, np = [0,'inf'], npEach = [0,'inf'], checkProbabilities = false, algorithm = "forward", verbose = true, threads = 1, fileName = nil, diversityFileName = nil, fileType = nil, treeNumber = nil)
    @fileName = fileName
    @fileType = fileType
    if @fileName == nil
      startTime = Time.now
      @treeReconstructed = false
      @branchRatesAssigned = false
      @sequencesEvolved = false
      @traitsEvolved = false
      @alignmentPositionsAssigned = false
      @fossilRecordAdded = false

      if verbose
        puts
        puts "------------------------------------------------phylsim.rb | Tree generation------------------------------------------------"
      end

      # A sample size is chosen for stochastic distributions for the case that not all three parameters lambda, mu, and treeOrigin are fixed.
      sampleSize = 1000

      # The specified values of 'lambda', 'mu', and 'treeOrigin' are checked
      parameters = [lambda,mu,treeOrigin]
      parametersOk = false
      until parametersOk == true
        parameters.size.times do |p|
          priorType = []
          if parameters[p].class == Fixnum or parameters[p].class == Float
            value = parameters[p].to_f
            raise "Please specify a positive value for lambda!" if value <= 0 and p == 0
            raise "Please specify a non-negative value for mu!" if value < 0 and p == 1
            raise "Please specify a positive value for treeOrigin (we here measure past time in positive numbers)!" if value <= 0 and p == 2
            priorType << "fixed"
          elsif parameters[p].class == String
            if parameters[p] =~ /uniform\(\s*([\d.]+)\s*,\s*([\d.]+)\s*\)/
              raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'lambda'!" if p == 0 and $1.to_f < 0 or $2.to_f <= 0
              raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'mu'!" if p == 1 and $1.to_f < 0 or $2.to_f <= 0
              raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'treeOrigin'!" if p == 2 and $1.to_f < 0 or $2.to_f <= 0
              if $1 < $2
                lower = $1.to_f
                upper = $2.to_f
                raise "The lower and upper bound of the uniform distribution for 'lambda' are too close! They should be at least 0.000001" if p == 0 and upper - lower < 0.000001
                raise "The lower and upper bound of the uniform distribution for 'mu' are too close! They should be at least 0.000001" if p == 1 and upper - lower < 0.000001
                raise "The lower and upper bound of the uniform distribution for 'treeOrigin' are too close! They should be at least 0.000001" if p == 2 and upper - lower < 0.000001
                value = []
                sampleSize.times do
                  value << lower + (upper-lower)*rand
                end
                priorType << "uniform" << lower << upper
              elsif $1 == $2
                value = $1.to_f
                priorType << "fixed"
              else
                lower = $2.to_f
                upper = $1.to_f
                raise "The lower and upper bound of the uniform distribution for 'lambda' are too close! They should differ by at least 0.000001" if p == 0 and upper - lower < 0.000001
                raise "The lower and upper bound of the uniform distribution for 'mu' are too close! They should differ by at least 0.000001" if p == 1 and upper - lower < 0.000001
                raise "The lower and upper bound of the uniform distribution for 'treeOrigin' are too close! They should differ by at least 0.000001" if p == 2 and upper - lower < 0.000001
                value = []
                sampleSize.times do
                  value << lower + (upper-lower)*rand
                end
                priorType << "uniform" << lower << upper
              end
            elsif parameters[p] =~ /normal\(\s*([\d.]+)\s*,\s*([\d.]+)\s*\)/
              raise "Please specify a non-negative mean for the normal distribution of 'lambda'!" if p == 0 and $1.to_f < 0
              raise "Please specify a non-negative mean for the normal distribution of 'mu'!" if p == 1 and $1.to_f < 0
              raise "Please specify a non-negative mean for the normal distribution of 'treeOrigin'!" if p == 2 and $1.to_f < 0
              raise "Please specify a positive standard deviation of 'lambda'!" if p == 0 and $2.to_f <= 0
              raise "Please specify a positive standard deviation of 'mu'!" if p == 1 and $2.to_f <= 0
              raise "Please specify a positive standard deviation of 'treeOrigin'!" if p == 2 and $2.to_f <= 0
              mean = $1.to_f
              stdev = $2.to_f
              normDist = Rubystats::NormalDistribution.new(mean,stdev)
              value = []
              sampleSize.times do
                # This is in order to truncate the normal distribution so that it is > 0
                done = false
                until done do
                  random = normDist.rng
                  if random > 0
                    value << random
                    done = true
                  end
                end
              end
              priorType << "normal" << mean << stdev
            elsif parameters[p] =~ /exponential\(\s*([\d.]+)\s*,\s*([\d.]+)\s*\)/
              raise "Please specify a non-negative offset for the exponential distribution of 'lambda'!" if p == 0 and $1.to_f < 0
              raise "Please specify a non-negative offset for the exponential distribution of 'mu'!" if p == 1 and $1.to_f < 0
              raise "Please specify a non-negative offset for the exponential distribution of 'treeOrigin'!" if p == 2 and $1.to_f < 0
              raise "Please specify a positive mean for the exponential distribution of 'lambda'!" if p == 0 and $2.to_f <= 0
              raise "Please specify a positive mean for the exponential distribution of 'mu'!" if p == 1 and $2.to_f <= 0
              raise "Please specify a positive mean for the exponential distribution of 'treeOrigin'!" if p == 2 and $2.to_f <= 0
              offset = $1.to_f
              mean = $2.to_f
              expDist = Rubystats::ExponentialDistribution.new(1/mean)
              value = []
              sampleSize.times do
                value << offset + expDist.rng
              end
              priorType << "exponential" << offset << mean
            elsif parameters[p] =~ /lognormal\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*\)/
              raise "Please specify a non-negative offset for the lognormal distribution of 'lambda'!" if p == 0 and $1.to_f < 0
              raise "Please specify a non-negative offset for the lognormal distribution of 'mu'!" if p == 1 and $1.to_f < 0
              raise "Please specify a non-negative offset for the lognormal distribution of 'treeOrigin'!" if p == 2 and $1.to_f < 0
              raise "Please specify a positive mean for the lognormal distribution of 'lambda'!" if p == 0 and $2.to_f <= 0
              raise "Please specify a positive mean for the lognormal distribution of 'mu'!" if p == 1 and $2.to_f <= 0
              raise "Please specify a positive mean for the lognormal distribution of 'treeOrigin'!" if p == 2 and $2.to_f <= 0
              raise "Please specify a positive standard deviation for the lognormal distribution of 'lambda'" if p == 0 and $3.to_f <= 0
              raise "Please specify a positive standard deviation for the lognormal distribution of 'mu'" if p == 1 and $3.to_f <= 0
              raise "Please specify a positive standard deviation for the lognormal distribution of 'treeOrigin'" if p == 2 and $3.to_f <= 0
              offset = $1.to_f
              mean = $2.to_f
              stdev = $3.to_f
              var = stdev**2
              normalmean = Math.log(mean) - 0.5*Math.log((var/mean**2) + 1)
              normalstdev = Math.sqrt(Math.log((var/mean**2) + 1))
              normDist = Rubystats::NormalDistribution.new(normalmean,normalstdev)
              value = []
              sampleSize.times do
                value << offset + Math.exp(normDist.rng)
              end
              priorType << "lognormal" << offset << mean << stdev
            end
          elsif parameters[p].class == Array and parameters[p].size == 2
            if parameters[p][0].class == Fixnum or parameters[p][0].class == Float
              if parameters[p][1].class == Fixnum or parameters[p][1].class == Float
                parameters[p][0],parameters[p][1] = parameters[p][1],parameters[p][0] if parameters[p][0] > parameters[p][1]
                raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'lambda'!" if p == 0 and parameters[p][0] < 0 or parameters[p][1] <= 0
                raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'mu'!" if p == 1 and parameters[p][0] < 0 or parameters[p][1] <= 0
                raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'treeOrigin'!" if p == 2 and parameters[p][0] < 0 or parameters[p][1] <= 0
                if parameters[p][0] == parameters[p][1]
                  value = parameters[p][0].to_f
                  priorType << "fixed"
                else
                  lower = parameters[p][0].to_f
                  upper = parameters[p][1].to_f
                  raise "The lower and upper bound of the uniform distribution for 'lambda' are too close! They should differ by at least 0.000001" if p == 0 and upper - lower < 0.000001
                  raise "The lower and upper bound of the uniform distribution for 'mu' are too close! They should differ by at least 0.000001" if p == 1 and upper - lower < 0.000001
                  raise "The lower and upper bound of the uniform distribution for 'treeOrigin' are too close! They should differ by at least 0.000001" if p == 2 and upper - lower < 0.000001
                  value = []
                  sampleSize.times do |i|
                    value << lower + (upper-lower)*rand
                  end
                  priorType << "uniform" << lower << upper
                end
              else
                raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'lambda'!" if p == 0
                raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'mu'!" if p == 1
                raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'treeOrigin'!" if p == 2
              end
            else
              raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'lambda'!" if p == 0
              raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'mu'!" if p == 1
              raise "Please specify a non-negative value for the lower bound, and a positive value for the upper bound of 'treeOrigin'!" if p == 2
            end
          else
            raise "An invalid value for 'lambda' has been specified." if p == 0
            raise "An invalid value for 'mu' has been specified." if p == 1
            raise "An invalid value for 'treeOrigin' has been specified: #{parameters[p]}." if p == 2
          end
          if p == 0
            @lambda = value
            @lambdaType = priorType
          end
          if p == 1
            @mu = value
            @muType = priorType
          end
          if p == 2
            @treeOrigin = value
            @treeOriginType = priorType
          end
        end
        
        parametersOk = true # tentatively

        # Check that @lambda[i] != @mu[i] is true for all i. If this was the case, it would cause problems with the calculation of probabilities below.
        # If @lambda[i] == @mu[i] for a certain i, repeat the entire procedure.
        if @lambdaType[0] == "fixed" and @muType[0] != "fixed"
          @mu.size.times do |i|
            parametersOk = false if @mu[i] == @lambda
          end
        elsif @lambdaType[0] != "fixed" and @muType[0] == "fixed"
          @lambda.size.times do |i|
            parametersOk = false if @lambda[i] == @mu
          end
        elsif @lambdaType[0] != "fixed" and @muType[0] != "fixed"
          @lambda.size.times do |i|
            parametersOk = false if @lambda[i] == @mu[i]
          end
        end
      end # until parametersOk

      # The specified value of 'present' is checked
      if present.class == Fixnum or present.class == Float
        @present = present.to_f
      else
        raise "Please specify a numerical value for 'present'!"
      end

      # The specified value of 'k' is checked.
      if k.class == Fixnum or k.class == Float
        @k = k.to_f
      else
        raise "Please specify a numerical value for 'k'!"
      end

      # Make sure the following requirements are met:
      # - All values of treeOrigin must be greater than the one specified for present.
      if @treeOrigin.class == Float
        raise "Please specify a value for 'treeOrigin' that is larger than the one for 'present' (we here measure past time in positive numbers)!" if @treeOrigin <= @present
      elsif @treeOrigin.class == Array
        @treeOrigin.each do |to|
          raise "Please specify a distribution for 'treeOrigin' so that all its values are larger than the one specified for 'present' (we here measure past time in positive numbers)!" if to <= @present
        end
      else
        raise "Parameter 'treeOrigin' has an unexpected format!" # This should never be called at all.
      end
      
      # The specified value of 'rootSplit' is checked
      rootSplit = false if rootSplit == "false"
      rootSplit = true if rootSplit == "true"
      unless rootSplit == false or rootSplit == true
        raise "Please specify either 'true' or 'false' for 'rootSplit'!"
      end
      @rootSplit = rootSplit
      
      # The specified value of 'np' is checked
      if np.class == Fixnum
        @minNp = np
        @maxNp = np
      elsif np.class == String
        if parameter =~ /uniform\(\s*(\d+)\s*,\s*(\d+)\s*\)/
          if $1 < $2
            raise "The lower bound for np must be a positive number." if $1 < 0
            @minNp = $1.to_i
            @maxNp = $2.to_i
          else
            raise "The lower bound for np must be a positive number." if $2 < 0
            @minNp = $2.to_i
            @maxNp = $1.to_i
          end
        elsif parameter =~ /uniform\(\s*(\d+)\s*,\s*[Ii][Nn][Ff]\s*\)/
          raise "The lower bound for np must be a positive number." if $1 < 0
          @minNp = $1.to_i
          @maxNp = nil
        elsif parameter =~ /uniform\(\s*[Ii][Nn][Ff]\s*,\s*(\d+)\s*\)/
          raise "The lower bound for np must be a positive number." if $1 < 0
          @minNp = $1.to_i
          @maxNp = nil
        else
          raise "An invalid value for np has been specified."
        end
      elsif np.class == Array and np.size == 2
        if np[0].class == Fixnum and np[1].class == Fixnum
          np[0],np[1] = np[1],np[0] if np[0] > np[1]
          raise "Both the lower and the upper bound for np must be positive numbers." if np[0] < 0
          @minNp = np[0]
          @maxNp = np[1]
        elsif np[0].class == Fixnum and np[1] =~ /\s*[Ii][Nn][Ff]\s*/
          raise "Both the lower and the upper bound for np must be positive numbers." if np[0] < 0
          @minNp = np[0]
          @maxNp = nil
        elsif np[1].class == Fixnum and np[0] =~ /\s*[Ii][Nn][Ff]\s*/
          raise "Both the lower and the upper bound for np must be positive numbers." if np[1] < 0
          @minNp = np[1]
          @maxNp = nil
        else
          raise "An invalid value for np has been specified."
        end
      elsif np == nil
        @minNp = 0
        @maxNp = nil
      else
        raise "An invalid value for np has been specified."
      end
      
      # The specified value of 'npEach' is checked
      if @rootSplit
        if npEach.class == Fixnum
          @minNpEach = npEach
          @maxNpEach = npEach
        elsif npEach.class == String
          if parameter =~ /uniform\(\s*(\d+)\s*,\s*(\d+)\s*\)/
            if $1 < $2
              raise "The lower bound for npEach must be a positive number." if $1 < 0
              @minNpEach = $1.to_i
              @maxNpEach = $2.to_i
            else
              raise "The lower bound for npEach must be a positive number." if $2 < 0
              @minNpEach = $2.to_i
              @maxNpEach = $1.to_i
            end
          elsif parameter =~ /uniform\(\s*(\d+)\s*,\s*[Ii][Nn][Ff]\s*\)/
            raise "The lower bound for npEach must be a positive number." if $1 < 0
            @minNpEach = $1.to_i
            @maxNpEach = nil
          elsif parameter =~ /uniform\(\s*[Ii][Nn][Ff]\s*,\s*(\d+)\s*\)/
            raise "The lower bound for npEach must be a positive number." if $1 < 0
            @minNpEach = $1.to_i
            @maxNpEach = nil
          else
            raise "An invalid value for np has been specified."
          end
        elsif npEach.class == Array and npEach.size == 2
          if npEach[0].class == Fixnum and npEach[1].class == Fixnum
            npEach[0],npEach[1] = npEach[1],npEach[0] if npEach[0] > npEach[1]
            raise "Both the lower and the upper bound for npEach must be positive numbers." if npEach[0] < 0
            @minNpEach = npEach[0]
            @maxNpEach = npEach[1]
          elsif npEach[0].class == Fixnum and npEach[1] =~ /\s*[Ii][Nn][Ff]\s*/
            raise "Both the lower and the upper bound for npEach must be positive numbers." if npEach[0] < 0
            @minNpEach = npEach[0]
            @maxNpEach = nil
          elsif npEach[1].class == Fixnum and npEach[0] =~ /\s*[Ii][Nn][Ff]\s*/
            raise "Both the lower and the upper bound for npEach must be positive numbers." if npEach[1] < 0
            @minNp = npEach[1]
            @maxNp = nil
          else
            raise "An invalid value for npEach has been specified."
          end
        elsif npEach == nil
          @minNpEach = 0
          @maxNpEach = nil
        else
          raise "An invalid value for npEach has been specified."
        end
      else
        warn "WARNING: Because rootSplit has been specified as 'false', parameter npEach will be ignored." unless npEach == nil or npEach == [0,'inf']
      end
      if @rootSplit == true
        unless @maxNp == nil
          raise "The specified values for np and npEach conflict with each other." if @minNpEach > 0.5*@maxNp
        end
        unless @maxNpEach == nil
          raise "The specified values for np and npEach conflict with each other." if @maxNpEach < 0.5*@minNp
        end
      end
      
      # The specified option to check probabilities is checked
      probabilitiesSampleSize = 0
      if checkProbabilities == true
        @checkProbabilities = true
        probabilitiesSampleSize = sampleSize
      elsif checkProbabilities == false
        @checkProbabilities = false
      elsif checkProbabilities.class == Fixnum or checkProbabilities.class == Float
        if checkProbabilities > 0
          @checkProbabilities = true
          probabilitiesSampleSize = checkProbabilities.to_i
        elsif checkProbabilities == 0
          @checkProbabilities = false
        else
          raise "Please specify a valid value for parameter 'checkProbabilities'! This could be 'true' (using #{sampleSize} as the sample size), 'false', or non-negative numeric (will be used as sample size if > 0)."
        end
      elsif checkProbabilities.class == String
        if checkProbabilities =~ /[Tt][Rr][Uu][Ee]/
          @checkProbabilities = true
          probabilitiesSampleSize = sampleSize
        elsif checkProbabilities =~ /[Ff][Aa][Ll][Ss][Ee]/
          @checkProbabilities = false
        elsif checkProbabilities.to_i > 0
          @checkProbabilities = true
          probabilitiesSampleSize = checkProbabilities.to_i
        else
          raise "Please specify a valid value for parameter 'checkProbabilities'! This could be 'true' (using #{sampleSize} as the sample size), 'false', or non-negative numeric (will be used as sample size if > 0)."
        end
      else
        raise "Please specify a valid value for parameter 'checkProbabilities'! This could be 'true' (using #{sampleSize} as the sample size), 'false', or non-negative numeric (will be used as sample size if > 0)."
      end

      # The specified algorithm is checked (algorithms 'forward' and 'treesim' are accepted)
      if algorithm.class == String
        if algorithm =~ /[Ff][Oo][Rr][Ww][Aa][Rr][Dd]/
          @algorithm = "forward"
        elsif algorithm =~ /[Tt][Rr][Ee][Ee][Ss][Ii][Mm]/
          @algorithm = "treesim"
        else
          raise "Please specify either 'forward' or 'treesim' as algorithm!"
        end
      end

      # The specified number of threads is checked
      if threads.to_i > 1
        @threads = threads
      else
        @threads = 1
      end

      # Unless all three parameters lambda, mu, and treeOrigin are fixed, all will be turned into arrays of size 'sampleSize', even if one or two of them are fixed. 
      if @lambda.class == Float and @mu.class == Float and @treeOrigin.class == Float
        allFixed = true
        t = @treeOrigin-@present
      else
        allFixed = false
        if @lambda.class == Float
          temp = []
          sampleSize.times do
            temp << @lambda
          end
          @lambda = temp
        end
        if @mu.class == Float
          temp = []
          sampleSize.times do
            temp << @mu
          end
          @mu = temp
        end
        if @treeOrigin.class == Float
          temp = []
          sampleSize.times do
            temp << @treeOrigin
          end
          @treeOrigin = temp
        end
        t = []
        sampleSize.times do |a|
          t << @treeOrigin[a]-@present
        end
      end

      # Some output is printed on the chosen model, parameters, and conditions
      if verbose
        puts
        puts "Model:"
        puts "Forward continuous-time birth-death, constant speciation and extinction rates, starting with a single lineage" if @algorithm == "forward" and @rootSplit == false
        puts "Continuous-time birth-death (Stadler 2011), constant speciation and extinction rates, starting with a single lineage" if @algorithm == "treesim" and @rootSplit == false
        puts "Forward continuous-time birth-death, constant speciation and extinction rates, starting with two lineages" if @algorithm == "forward" and @rootSplit == true
        puts "Continuous-time birth-death (Stadler 2011), constant speciation and extinction rates, starting with two lineages" if @algorithm == "treesim" and @rootSplit == true
        puts
        puts "Parameters:"
        if allFixed
          puts "Speciation rate: lambda = #{@lambda}"
          puts "Extinction rate: mu = #{@mu}"
          puts "Time of tree origin: treeOrigin = #{@treeOrigin}"
        else
          puts "Speciation rate: lambda = #{@lambda[0]}" if @lambdaType[0] == "fixed"
          puts "Speciation rate: lambda ~ uniform(#{@lambdaType[1]}, #{@lambdaType[2]})" if @lambdaType[0] == "uniform"
          puts "Speciation rate: lambda ~ normal(mean = #{@lambdaType[1]}, stdev = #{@lambdaType[2]}); lambda > 0" if @lambdaType[0] == "normal"
          puts "Speciation rate: lambda ~ exponential(offset = #{@lambdaType[1]}, mean = #{@lambdaType[2]})" if @lambdaType[0] == "exponential"
          puts "Speciation rate: lambda ~ lognormal(offset = #{@lambdaType[1]}, mean = #{@lambdaType[2]}, stdev = #{@lambdaType[3]})" if @lambdaType[0] == "lognormal"
          puts "Extinction rate: mu = #{@mu[0]}" if @muType[0] == "fixed"
          puts "Extinction rate: mu ~ uniform(#{@muType[1]}, #{@muType[2]})" if @muType[0] == "uniform"
          puts "Extinction rate: mu ~ normal(mean = #{@muType[1]}, stdev = #{@muType[2]}); mu > 0" if @muType[0] == "normal"
          puts "Extinction rate: mu ~ exponential(offset = #{@muType[1]}, mean = #{@muType[2]})" if @muType[0] == "exponential"
          puts "Extinction rate: mu ~ lognormal(offset = #{@muType[1]}, mean = #{@muType[2]}, stdev = #{@muType[3]})" if @muType[0] == "lognormal"
          puts "Time of tree origin: treeOrigin = #{@treeOrigin[0]}" if @treeOriginType[0] == "fixed"
          puts "Time of tree origin: treeOrigin ~ uniform(#{@treeOriginType[1]}, #{@treeOriginType[2]})" if @treeOriginType[0] == "uniform"
          puts "Time of tree origin: treeOrigin ~ normal(mean = #{@treeOriginType[1]}, stdev = #{@treeOriginType[2]}); treeOrigin > 0" if @treeOriginType[0] == "normal"
          puts "Time of tree origin: treeOrigin ~ exponential(offset = #{@treeOriginType[1]}, mean = #{@treeOriginType[2]})" if @treeOriginType[0] == "exponential"
          puts "Time of tree origin: treeOrigin ~ lognormal(offset = #{@treeOriginType[1]}, mean = #{@treeOriginType[2]}, stdev = #{@treeOriginType[3]})" if @treeOriginType[0] == "lognormal"
        end
        puts
        puts "Condition on full tree:"
        if @maxNp == nil and @minNp == 0
          puts "none"
        elsif @maxNp == nil and @minNp > 1
          puts "extant species > #{@minNp-1}"
        elsif @maxNp != nil
          puts "#{@minNp-1} < extant species < #{@maxNp+1}"
        end
        if rootSplit
          puts
          puts "Condition on each partial tree descending from root:"
          if @maxNpEach == nil and @minNpEach == 0
            puts "none"
          elsif @maxNpEach == nil and @minNpEach > 0
            puts "extant species in partial tree > #{@minNpEach-1}"
          elsif @maxNpEach != nil
            puts "#{@minNpEach-1} < extant species in partial tree < #{@maxNpEach+1}"
          end
        end
      
        if @checkProbabilities
          # If parameters lambda, mu, and treeOrigin are fixed, then the probability that the conditions are met are directly calculated (using Kendall 1949 if mu > 0 or Yule 1925 if mu = 0).
          # If only one of the three parameters is stochastic, then the mean probability is estimated using a sample of size 'probabilitiesSampleSize'
          prob = nil
          p = []
          if @rootSplit == false and allFixed == true
            if @maxNp == nil
              if @minNp == 0
                prob = 1.0
              elsif @minNp == 1
                p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                prob = 1-p[0]
              elsif @minNp > 1
                p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))                    # p(N=0|lambda,mu,t) after Kendall 1949, same result as p(N=0|lambda,t) after Yule 1925
                p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2      # p(N=1|lambda,mu,t) after Kendall 1949, same result as p(N=0|lambda,t) after Yule 1925
                2.upto(@minNp) do |n|
                  p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                  p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                end
                prob = 1-p[0..@minNp-1].sum
              end
            elsif @maxNp == 0
              p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
              prob = p[0]
            elsif @maxNp == 1
              if @minNp == 0
                p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                prob = p[0]+p[1]
              else
                p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                prob = p[1]
              end
            elsif @maxNp > 1
              p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
              p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
              2.upto(@maxNp) do |n|
                p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
              end
              prob = p[@minNp..@maxNp].sum
            end

          elsif @rootSplit == false and allFixed == false
            if @maxNp == nil
              if @minNp == 0
              
                prob = 1.0
              
              elsif @minNp == 1
              
                # Prepare arrays for p up to the upper limit of calculation (0)
                p[0] = []
                # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                probabilitiesSampleSize.times do
                  random = rand(sampleSize)
                  p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                end
                # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                prob = []
                probabilitiesSampleSize.times do |s|
                  prob << 1-p[0][s]
                end
              
              elsif @minNp > 1
              
                # Prepare arrays for p up to the upper limit of calculation (@minNp)
                p[0] = []
                p[1] = []
                2.upto(@minNp) do |n|
                  p[n] = []
                end
                # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                probabilitiesSampleSize.times do
                  random = rand(sampleSize)
                  p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                  p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                  2.upto(@minNp) do |n|
                    p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                    p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                  end 
                end
                # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                prob = []
                probabilitiesSampleSize.times do |s|
                  sum = 0
                  0.upto(@minNp-1) do |x|
                    sum += p[x][s]
                  end
                  prob << 1-sum
                end
              
              end
            
            elsif @maxNp == 0

              # Prepare arrays for p up to the upper limit of calculation (0)
              p[0] = []
              # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
              probabilitiesSampleSize.times do |s|
                random = rand(sampleSize)
                p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
              end
              # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
              prob = []
              probabilitiesSampleSize.times do |s|
                prob << p[0][s]
              end
            
            elsif @maxNp == 1
            
              if @minNp == 0
              
                # Prepare arrays for p up to the upper limit of calculation (1)
                p[0] = []
                p[1] = []
                # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                probabilitiesSampleSize.times do |s|
                  random = rand(sampleSize)
                  p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                  p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                end
                # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                prob = []
                probabilitiesSampleSize.times do |s|
                  prob = p[0][s]+p[1][s]
                end
              
              else
              
                # Prepare arrays for p up to the upper limit of calculation (@minNp)
                p[1] = []
                # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                probabilitiesSampleSize.times do
                  random = rand(sampleSize)
                  p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                end
                # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                prob = []
                probabilitiesSampleSize.times do |s|
                  prob << p[1][s]
                end
              
              end
            elsif @maxNp > 1
            
              # Prepare arrays for p up to the upper limit of calculation (@maxNp)
              p[0] = []
              p[1] = []
              2.upto(@maxNp) do |n|
                p[n] = []
              end
              # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
              probabilitiesSampleSize.times do
                random = rand(sampleSize)
                p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                2.upto(@maxNp) do |n|
                  p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                  p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                end
              end
              # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
              prob = []
              probabilitiesSampleSize.times do |s|
                sum = 0
                @minNp.upto(@maxNp) do |x|
                  sum += p[x][s]
                end
                prob << sum
              end
            
            end

          elsif @rootSplit == true and allFixed == true

            if @maxNp == nil
              if @minNp == 0
                if @maxNpEach == nil
                  if @minNpEach == 0
                    prob = 1.0
                  elsif @minNpEach == 1
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    prob = (1-p[0])*(1-p[0])
                  elsif @minNpEach > 1
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    2.upto(@minNpEach) do |n|
                      p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                      p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                    end
                    prob = (1-p[0..@minNpEach-1].sum)*(1-p[0..@minNpEach-1].sum)
                  end
                elsif @maxNpEach == 0
                  p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                  prob = p[0]*p[0]
                elsif @maxNpEach == 1
                  if @minNpEach == 0
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    prob = p[0]*p[0] + p[0]*p[1] + p[1]*p[0] + p[1]*p[1]
                  else
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    prob = p[1]*p[1]
                  end
                elsif @maxNpEach > 1
                  p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                  p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                  2.upto(@maxNpEach) do |n|
                    p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                    p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                  end
                  prob = p[@minNpEach..@maxNpEach].sum*p[@minNpEach..@maxNpEach].sum
                end
              elsif @minNp == 1
                if @maxNpEach == nil
                  if @minNpEach == 0
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    prob = 1 - p[0]*p[0]
                  elsif @minNpEach == 1
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    prob = (1-p[0])*(1-p[0])
                  elsif @minNpEach > 1
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    2.upto(@minNpEach) do |n|
                      p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                      p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                    end
                    prob = (1-p[0..@minNpEach-1].sum)*(1-p[0..@minNpEach-1].sum)
                  end
                elsif @maxNpEach == 1
                  if @minNpEach == 0
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    prob = p[0]*p[1] + p[1]*p[0] + p[1]*p[1]
                  elsif @minNpEach == 1
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    prob = p[1]*p[1]
                  end
                elsif @maxNpEach > 1
                  if @minNpEach == 0
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    2.upto(@maxNpEach) do |n|
                      p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                      p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                    end
                    prob = p[0..@maxNpEach].sum*p[0..@maxNpEach].sum - p[0]*p[0]
                  elsif @minNpEach > 1
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    2.upto(@maxNpEach) do |n|
                      p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                      p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                    end
                    prob = p[@minNpEach..@maxNpEach].sum*p[@minNpEach..@maxNpEach].sum
                  end
                end
              elsif @minNp > 1
                if @maxNpEach == nil
                  if @minNpEach == 0
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    2.upto(@minNp) do |n|
                      p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                      p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                    end
                    negprob = 0
                    0.upto(@minNp-1) do |q1|
                      0.upto(@minNp-1-q1) do |q2|
                        negprob += p[q1]*p[q2]
                      end
                    end
                    prob = 1 - negprob
                  elsif @minNpEach == 1
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    2.upto(@minNp) do |n|
                      p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                      p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                    end
                    negprob = p[0] + p[0] - p[0]*p[0]
                    1.upto(@minNp-1) do |q1|
                      1.upto(@minNp-1-q1) do |q2|
                        negprob += p[q1]*p[q2]
                      end
                    end
                    prob = 1 - negprob
                  elsif @minNpEach > 1   # @maxNp == nil, @minNp > 1, @maxNpEach == nil, @minNpEach > 1
                    p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                    p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                    2.upto(@minNp) do |n|
                      p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                      p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                    end
                    negprob = 0
                    0.upto([@minNpEach,@minNp].max-1) do |q1|
                      0.upto([@minNpEach,@minNp].max-1) do |q2|
                        negprob += p[q1]*p[q2] if q1+q2 < @minNp or q1 < @minNp or q2 < @minNp
                      end
                    end
                    prob = 1-negprob
                  end
                elsif @maxNpEach == 1
                  p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                  prob = p[1]*p[1]
                elsif @maxNpEach > 1
                  p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                  p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                  2.upto(@maxNpEach) do |n|
                    p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                    p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
                  end
                  prob = 0
                  @minNpEach.upto(@maxNpEach) do |q1|
                    @minNpEach.upto(@maxNpEach) do |q2|
                      prob += p[q1]*p[q2] unless q1+q2 < @minNp
                    end
                  end
                end
              end
            elsif @maxNp == 0
              p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
              prob = p[0]*p[0]
            elsif @maxNp == 1
              if @minNp == 0
                if @maxNpEach == 0
                  p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                  prob = p[0]*p[0]
                else
                  p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                  p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                  prob = p[0]*p[0] + p[0]*p[1] + p[1]*p[0]
                end
              else
                p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
                p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
                prob = p[0]*p[1] + p[1]*p[0]
              end
            elsif @maxNp > 1
              p[0] = (@mu*(1-Math.exp(-t*(@lambda-@mu))))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))
              p[1] = (((@lambda-@mu)**2)*Math.exp(-t*(@lambda-@mu)))/(@lambda-@mu*Math.exp(-t*(@lambda-@mu)))**2
              2.upto(@maxNp) do |n|
                p << ((@lambda/@mu)**(n-1))*p[1]*((p[0])**(n-1)) if mu > 0                                            # p(N=n|lambda,mu,t) after Kendall 1949
                p << (Math.exp(-@lambda*t))*((1-(Math.exp(-@lambda*t)))**(n-1)) if mu == 0                            # p(N=n|lambda,t) after Yule 1925
              end
              if @minNp == 0
                if @maxNpEach == nil
                  prob = 0
                  @minNpEach.upto(@maxNp) do |q1|
                    @minNpEach.upto(@maxNp) do |q2|
                      prob += p[q1]*p[q2] unless q1+q2 > @maxNp
                    end
                  end
                elsif @maxNpEach == 0
                  prob = p[0]*p[0]
                elsif @maxNpEach == 1
                  if @minNpEach == 0
                    prob = p[0]*p[0] + p[0]*p[1] + p[1]*p[0] + p[1]*p[1]
                  else
                    prob = p[1]*p[1]
                  end
                elsif @maxNpEach > 1
                  prob = 0
                  @minNpEach.upto(@maxNpEach) do |q1|
                    @minNpEach.upto(@maxNpEach) do |q2|
                      prob += p[q1]*p[q2] unless q1+q2 > @maxNp
                    end
                  end
                end
              elsif @minNp == 1
                if @maxNpEach == nil
                  prob = 0
                  @minNpEach.upto(@maxNp) do |q1|
                    @minNpEach.upto(@maxNp) do |q2|
                      prob += p[q1]*p[q2] unless q1+q2 > @maxNp or q1+q2 < @minNp
                    end
                  end
                elsif @maxNpEach == 1
                  if @minNpEach == 0
                    prob = p[0]*p[1] + p[1]*p[0] + p[1]*p[1]
                  else
                    prob = p[1]*p[1]
                  end
                elsif @maxNpEach > 1
                  prob = 0
                  @minNpEach.upto(@maxNpEach) do |q1|
                    @minNpEach.upto(@maxNpEach) do |q2|
                      prob += p[q1]*p[q2] unless q1+q2 > @maxNp or q1+q2 < @minNp
                    end
                  end
                end
              elsif @minNp > 1
                if @maxNpEach == nil
                  prob = 0
                  @minNpEach.upto(@maxNp) do |q1|
                    @minNpEach.upto(@maxNp) do |q2|
                      prob += p[q1]*p[q2] unless q1+q2 > @maxNp or q1+q2 < @minNp
                    end
                  end
                elsif @maxNpEach == 1
                  prob = p[1]*p[1]
                elsif @maxNpEach > 1
                  prob = 0
                  @minNpEach.upto(@maxNpEach) do |q1|
                    @minNpEach.upto(@maxNpEach) do |q2|
                      prob += p[q1]*p[q2] unless q1+q2 > @maxNp or q1+q2 < @minNp
                    end
                  end
                end
              end
            end

          elsif @rootSplit == true and allFixed == false

            if @maxNp == nil
              if @minNp == 0
                if @maxNpEach == nil
                  if @minNpEach == 0

                    prob = 1.0

                  elsif @minNpEach == 1
                  
                    # Prepare arrays for p up to the upper limit of calculation (0)
                    p[0] = []
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << (1-p[0][s])*(1-p[0][s])
                    end
                  
                  elsif @minNpEach > 1

                    # Prepare arrays for p up to the upper limit of calculation (@minNpEach)
                    p[0] = []
                    p[1] = []
                    2.upto(@minNpEach) do |n|
                      p[n] = []
                    end
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                      2.upto(@minNpEach) do |n|
                        p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                        p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                      end
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      sum = 0
                      0.upto(@minNpEach-1) do |x|
                        sum += p[x][s]
                      end
                      prob << (1-sum)*(1-sum)
                    end

                  end

                elsif @maxNpEach == 0

                  # Prepare arrays for p up to the upper limit of calculation (0)
                  p[0] = []
                  # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                  probabilitiesSampleSize.times do
                    random = rand(sampleSize)
                    p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                  end
                  # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    prob << p[0][s]*p[0][s]
                  end
                
                elsif @maxNpEach == 1

                  if @minNpEach == 0

                    # Prepare arrays for p up to the upper limit of calculation (1)
                    p[0] = []
                    p[1] = []
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << p[0][s]*p[0][s] + p[0][s]*p[1][s] + p[1][s]*p[0][s] + p[1][s]*p[1][s]
                    end
                  
                  else

                    # Prepare arrays for p up to the upper limit of calculation (1)
                    p[1] = []
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << p[1][s]*p[1][s]
                    end
                  
                  end
                
                elsif @maxNpEach > 1
                
                  # Prepare arrays for p up to the upper limit of calculation (@maxNpEach)
                  p[0] = []
                  p[1] = []
                  2.upto(@maxNpEach) do |n|
                    p[n] = []
                  end
                  # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                  probabilitiesSampleSize.times do
                    random = rand(sampleSize)
                    p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                    p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                    2.upto(@maxNpEach) do |n|
                      p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                      p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                    end
                  end
                  # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    sum = 0
                    @minNpEach.upto(@maxNpEach) do |x|
                      sum += p[x][s]
                    end
                    prob << sum*sum
                  end  

                end

              elsif @minNp == 1

                if @maxNpEach == nil

                  if @minNpEach == 0

                    # Prepare arrays for p up to the upper limit of calculation (0)
                    p[0] = []
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << 1 - p[0][s]*p[0][s]
                    end

                  elsif @minNpEach == 1

                    # Prepare arrays for p up to the upper limit of calculation (0)
                    p[0] = []
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << (1-p[0][s])*(1-p[0][s])
                    end

                  elsif @minNpEach > 1

                    # Prepare arrays for p up to the upper limit of calculation (@minNpEach)
                    p[0] = []
                    p[1] = []
                    2.upto(@minNpEach) do |n|
                      p[n] = []
                    end
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                      2.upto(@minNpEach) do |n|
                        p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                        p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                      end
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      sum = 0
                      0.upto(@minNpEach-1) do |x|
                        sum += p[x][s]
                      end
                      prob << (1-sum)*(1-sum)
                    end
                  
                  end
                
                elsif @maxNpEach == 1
                
                  if @minNpEach == 0
                  
                    # Prepare arrays for p up to the upper limit of calculation (1)
                    p[0] = []
                    p[1] = []
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << p[0][s]*p[1][s] + p[1][s]*p[0][s] + p[1][s]*p[1][s]
                    end

                  end
                
                elsif @maxNpEach > 1
                
                  if @minNpEach == 0 # and maxNpEach > 1 and minNp == 1 and maxNP == nil
                  
                    # Prepare arrays for p up to the upper limit of calculation (@maxNpEach)
                    p[0] = []
                    p[1] = []
                    2.upto(@maxNpEach) do |n|
                      p[n] = []
                    end
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                      2.upto(@maxNpEach) do |n|
                        p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                        p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                      end
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      sum = 0
                      0.upto(@macNpEach) do |x|
                        sum += p[x][s]
                      end
                      prob << sum*sum - p[0][s]*p[0][s]
                    end

                  elsif @minNpEach > 0
                  
                    # Prepare arrays for p up to the upper limit of calculation (@maxNpEach)
                    p[0] = []
                    p[1] = []
                    2.upto(@maxNpEach) do |n|
                      p[n] = []
                    end
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                      2.upto(@maxNpEach) do |n|
                        p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                        p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                      end
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      sum = 0
                      @minNpEach.upto(@macNpEach) do |x|
                        sum += p[x][s]
                      end
                      prob << sum*sum
                    end

                  end
                
                end
              
              elsif @minNp > 1
              
                if @maxNpEach == nil
                
                  if @minNpEach == 0
                  
                    # Prepare arrays for p up to the upper limit of calculation (@minNp)
                    p[0] = []
                    p[1] = []
                    2.upto(@minNp) do |n|
                      p[n] = []
                    end
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                      2.upto(@minNp) do |n|
                        p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                        p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                      end
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      negsum = 0
                      0.upto(@minNp-1) do |q1|
                        0.upto(@minNp-1-q1) do |q2|
                          negsum += p[q1][s]*p[q2][s]
                        end
                      end
                      prob << 1-negsum
                    end

                  elsif @minNpEach == 1
                  
                    # Prepare arrays for p up to the upper limit of calculation (@minNp)
                    p[0] = []
                    p[1] = []
                    2.upto(@minNp) do |n|
                      p[n] = []
                    end
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                      2.upto(@minNp) do |n|
                        p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                        p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                      end
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      negsum = p[0][s] + p[0][s] - p[0][s]*p[0][s]
                      1.upto(@minNp-1) do |q1|
                        1.upto(@minNp-1-q1) do |q2|
                          negsum += p[q1][s]*p[q2][s]
                        end
                      end
                      prob << 1-negsum
                    end

                  elsif @minNpEach > 1   # @maxNp == nil, @minNp > 1, @maxNpEach == nil, @minNpEach > 1
                  
                    # Prepare arrays for p up to the upper limit of calculation (@minNp)
                    p[0] = []
                    p[1] = []
                    2.upto(@minNp) do |n|
                      p[n] = []
                    end
                    # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                    probabilitiesSampleSize.times do
                      random = rand(sampleSize)
                      p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                      p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                      2.upto(@minNp) do |n|
                        p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                        p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                      end
                    end
                    # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      negsum = 0
                      0.upto([@minNpEach,@minNp].max-1) do |q1|
                        0.upto([@minNpEach,@minNp].max-1) do |q2|
                          negsum += p[q1][s]*p[q2][s] if q1+q2 < @minNp or q1 < @minNp or q2 < @minNp
                        end
                      end
                      prob << 1-negsum
                    end
                  
                  end
                
                elsif @maxNpEach == 1
                
                  # Prepare arrays for p up to the upper limit of calculation (1)
                  p[1] = []
                  # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                  probabilitiesSampleSize.times do
                    random = rand(sampleSize)
                    p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                  end
                  # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    prob << p[1][s]*p[1][s]
                  end
                
                elsif @maxNpEach > 1

                  # Prepare arrays for p up to the upper limit of calculation (@maxNpEach)
                  p[0] = []
                  p[1] = []
                  2.upto(@maxNpEach) do |n|
                    p[n] = []
                  end
                  # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                  probabilitiesSampleSize.times do
                    random = rand(sampleSize)
                    p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                    p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                    2.upto(@maxNpEach) do |n|
                      p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949
                      p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                    end
                  end
                  # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    sum = 0
                    @minNpEach.upto(@maxNpEach) do |q1|
                      @minNpEach.upto(@maxNpEach) do |q2|
                        sum += p[q1][s]*p[q2][s] unless q1+q2 < @minNp
                      end
                    end
                    prob << sum
                  end

                end

              end

            elsif @maxNp == 0

              # Prepare arrays for p up to the upper limit of calculation (0)
              p[0] = []
              # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
              probabilitiesSampleSize.times do
                random = rand(sampleSize)
                p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
              end
              # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
              prob = []
              probabilitiesSampleSize.times do |s|
                prob << p[0][s]*p[0][s]
              end

            elsif @maxNp == 1
            
              if @minNp == 0
              
                if @maxNpEach == 0
                
                  # Prepare arrays for p up to the upper limit of calculation (0)
                  p[0] = []
                  # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                  probabilitiesSampleSize.times do
                    random = rand(sampleSize)
                    p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                  end
                  # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    prob << p[0][s]*p[0][s]
                  end
                
                else
                
                  # Prepare arrays for p up to the upper limit of calculation (1)
                  p[0] = []
                  p[1] = []
                  # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                  probabilitiesSampleSize.times do
                    random = rand(sampleSize)
                    p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                    p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                  end
                  # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    prob << p[0][s]*p[0][s] + p[0][s]*p[1][s] + p[1][s]*p[0][s]
                  end
                
                end
              
              else
              
                # Prepare arrays for p up to the upper limit of calculation (1)
                p[0] = []
                p[1] = []
                # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
                probabilitiesSampleSize.times do
                  random = rand(sampleSize)
                  p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                  p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                end
                # Now that the necessary arrays are filled, calculate the 'probabilitiesSampleSize' probabilities to match the condition and fill new array 'prob' with it.
                prob = []
                probabilitiesSampleSize.times do |s|
                  prob << p[0]*p[1] + p[1]*p[0]
                end
              
              end
            
            elsif @maxNp > 1

              # @maxNp is necessarily larger than @minNp, @minNpEach, and @maxNpEach, and thus probabilities are here calculated up to @maxNp irrespective of the other parameters.
              # Prepare arrays for p up to the upper limit of calculation (@maxNp)
              p[0] = []
              p[1] = []
              2.upto(@maxNp) do |n|
                p[n] = []
              end
              # Calculate the probabilities 'probabilitiesSampleSize' times and fill arrays with them
              probabilitiesSampleSize.times do
                random = rand(sampleSize)
                p[0] << (@mu[random]*(1-Math.exp(-t[random]*(@lambda[random]-@mu[random]))))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))
                p[1] << (((@lambda[random]-@mu[random])**2)*Math.exp(-t[random]*(@lambda[random]-@mu[random])))/(@lambda[random]-@mu[random]*Math.exp(-t[random]*(@lambda[random]-@mu[random])))**2
                2.upto(@maxNp) do |n|
                  p[n] << ((@lambda[random]/@mu[random])**(n-1))*p[1][-1]*((p[0][-1])**(n-1)) if @mu[random] > 0                          # p(N=n|lambda,mu,t) after Kendall 1949. For large values of n, this can easily lead to NaN or Infinity.
                  p[n] << (Math.exp(-@lambda[random]*t[random]))*((1-(Math.exp(-@lambda[random]*t[random])))**(n-1)) if @mu[random] == 0  # p(N=n|lambda,t) after Yule 1925
                end
              end
              if @minNp == 0
              
                if @maxNpEach == nil
                
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    sum = 0
                    @minNpEach.upto(@maxNp) do |q1|
                      @minNpEach.upto(@maxNp) do |q2|
                        sum += p[q1][s]*p[q2][s] unless q1+q2 > @maxNp
                      end
                    end
                    prob << sum
                  end
                
                elsif @maxNpEach == 0
                
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    prob << p[0][s]*p[0][s]
                  end
                
                elsif @maxNpEach == 1
                
                  if @minNpEach == 0
                  
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << p[0][s]*p[0][s] + p[0][s]*p[1][s] + p[1][s]*p[0][s] + p[1][s]*p[1][s]
                    end
                  
                  else
                  
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << p[1][s]*p[1][s]
                    end
                  
                  end
                
                elsif @maxNpEach > 1
                
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    sum = 0
                    @minNpEach.upto(@maxNpEach) do |q1|
                      @minNpEach.upto(@maxNpEach) do |q2|
                        sum += p[q1][s]*p[q2][s] unless q1+q2 > @maxNp
                      end
                    end
                    prob << sum
                  end
                
                end
              
              elsif @minNp == 1
              
                if @maxNpEach == nil
                
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    sum = 0
                    @minNpEach.upto(@maxNp) do |q1|
                      @minNpEach.upto(@maxNp) do |q2|
                        sum += p[q1][s]*p[q2][s] unless q1+q2 > @maxNp or q1+q2 < @minNp
                      end
                    end
                    prob << sum
                  end
                
                elsif @maxNpEach == 1
                
                  if @minNpEach == 0
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << p[0][s]*p[1][s] + p[1][s]*p[0][s] + p[1][s]*p[1][s]
                    end
                  
                  else
                  
                    prob = []
                    probabilitiesSampleSize.times do |s|
                      prob << p[1][s]*p[1][s]
                    end
                  
                  end
                
                elsif @maxNpEach > 1
                
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    sum = 0
                    @minNpEach.upto(@maxNpEach) do |q1|
                      @minNpEach.upto(@maxNpEach) do |q2|
                        sum += p[q1][s]*p[q2][s] unless q1+q2 > @maxNp or q1+q2 < @minNp
                      end
                    end
                    prob << sum
                  end
                
                end
              
              elsif @minNp > 1
              
                if @maxNpEach == nil
                
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    sum = 0
                    @minNpEach.upto(@maxNp) do |q1|
                      @minNpEach.upto(@maxNp) do |q2|
                        sum += p[q1][s]*p[q2][s] unless q1+q2 > @maxNp or q1+q2 < @minNp
                      end
                    end
                    prob << sum
                  end
                
                elsif @maxNpEach == 1
                
                  prob = []
                  probabilitiesSampleSize.times do |s|
                    prob << p[1][s]*p[1][s]
                  end
                
                elsif @maxNpEach > 1

                  prob = []
                  probabilitiesSampleSize.times do |s|
                    sum = 0
                    @minNpEach.upto(@maxNpEach) do |q1|
                      @minNpEach.upto(@maxNpEach) do |q2|
                        sum += p[q1][s]*p[q2][s] unless q1+q2 > @maxNp or q1+q2 < @minNp
                      end
                    end
                    prob << sum
                  end
                
                end
              
              end
            
            end

          end

          puts
          if rootSplit
            if allFixed
              puts "Probability to match conditions:"
            else
              puts "Probability to match conditions (estimated):"
            end
          else
            if allFixed
              puts "Probability to match condition:"
            else
              puts "Probability to match condition (estimated):"
            end
          end
          if prob.class == Array
            prob = prob.sum/probabilitiesSampleSize
          end
          if prob.to_s == "Infinity" or prob.to_s == "NaN" # this can happen with large values for maxNp.
            puts "Could not be calculated."
          else
            puts "#{prob.round(4)}"
            puts
            puts "Expected number of trials:"
            puts "#{(1/prob).round}"
          end

        end # if checkProbabilities

      end # if verbose
      
      # Prepare for the loop that runs until a tree is found that matches the conditions
      treeMeetsCondition = false
      dumpFileName = ".dump.txt"
      dumpFileNames = []
      @threads.times do |t|
        dumpFileNames << "#{dumpFileName.chomp(".txt")}_r#{(t.to_s.rjust(@threads)).gsub(" ","0")}.txt"
      end
      puts if verbose

      # Prepare the tree search loop.
      trial = 1
      done = false
      until done
        if verbose
          if @threads > 1
            print "\rRunning trials #{trial}-#{trial+@threads-1}..."
          else
            print "\rRunning trial #{trial}..."
          end
        end
        pids = []
        @threads.times do |t|
          Process.fork do
            trial = trial + t
# trial = trial + 1
            # Define lambda, mu, and treeOrigin for this tree finding trial.
            if allFixed == true
              trialLambda = @lambda
              trialMu = @mu
              trialTreeOrigin = @treeOrigin
            else
              random = rand(sampleSize)
              trialLambda = @lambda[random]
              trialMu = @mu[random]
              trialTreeOrigin = @treeOrigin[random]
            end

            # Initiate branch array etc
            @branch = []
            branchIdCounter = 0
            treeComplete = false
            treeTooLarge = false

            # Use the specified algorithm ("forward" or "treesim") 
            if @algorithm == "forward"

              # Initiate the branch length distribution.
              branchLengthDist = Rubystats::ExponentialDistribution.new(trialLambda+trialMu)

              # Produce first branch.
              if rand < trialLambda/(trialLambda+trialMu) # the new branch will terminate with a speciation event
                @branch << Branch.new(id = "b0", origin = trialTreeOrigin,termination = trialTreeOrigin-branchLengthDist.rng,parentId = "treeOrigin",daughterId = ["unborn","unborn"],endCause = "speciation")
              else
                @branch << Branch.new(id = "b0", origin = trialTreeOrigin,termination = trialTreeOrigin-branchLengthDist.rng,parentId = "treeOrigin",daughterId = ["none","none"],endCause = "extinction")
              end
              branchIdCounter += 1
              
              # Produce a second branch starting at @treeOrigin if the first speciation event is set to occur at tree origin.
              if @rootSplit
                if rand < trialLambda/(trialLambda+trialMu) # the new branch will terminate with a speciation event
                  @branch << Branch.new(id = "b1", origin = trialTreeOrigin,termination = trialTreeOrigin-branchLengthDist.rng,parentId = "treeOrigin",daughterId = ["unborn","unborn"],endCause = "speciation")
                else
                  @branch << Branch.new(id = "b1", origin = trialTreeOrigin,termination = trialTreeOrigin-branchLengthDist.rng,parentId = "treeOrigin",daughterId = ["none","none"],endCause = "extinction")
                end
                branchIdCounter += 1
              end
              
              # Produce all remaining branches until either all branches end with extinctions, or all branches have reached the present
              until treeComplete or treeTooLarge
                change = false
                @branch.size.times do |b| # go through all branches
                  if @branch[b].daughterId == ["unborn","unborn"] and @branch[b].termination > @present # if a branch terminated with a speciation event in the past, then add the two daughter branches
                    currentId = ["b#{branchIdCounter}","b#{branchIdCounter+1}"]
                    2.times do |d|
                      if rand < trialLambda/(trialLambda+trialMu) # the new branch will terminate with a speciation event
                        @branch << Branch.new(id = currentId[d], origin = @branch[b].termination, termination = origin-branchLengthDist.rng,parentId = @branch[b].id,daughterId = ["unborn","unborn"],endCause = "speciation")
                      else # the new branch will terminate with an extinction event
                        @branch << Branch.new(id = currentId[d], origin = @branch[b].termination, termination = origin-branchLengthDist.rng,parentId = @branch[b].id,daughterId = ["none","none"],endCause = "extinction")
                      end
                      branchIdCounter += 1
                    end
                    change = true
                    branch[b].updateDaughterId(currentId) # update the parent branch with the Ids of the two newly created daughter branches
                  end
                end
                treeComplete = true if change == false
                treeTooLarge = true if @branch.size > 10 * @maxNp unless @maxNp == nil
              end

            elsif algorithm == "treesim"

              # Initiate the current time and the counter of current branches.
              currentTime = trialTreeOrigin
              currentBranch = []

              # Produce first branch.
              @branch << Branch.new(id = "b0", origin = trialTreeOrigin, termination = nil, parentId = "treeOrigin", daughterId = ["unborn","unborn"], endCause = nil)
              branchIdCounter += 1

              # Produce a second branch starting at @treeOrigin if the first speciation event is set to occur at tree origin.
              if @rootSplit
                @branch << Branch.new(id = "b1", origin = trialTreeOrigin, termination = nil, parentId = "treeOrigin", daughterId = ["unborn","unborn"], endCause = nil)
                branchIdCounter += 1
              end

              # Add the first one or two branches to array currentBranch.
              @branch.each {|b| currentBranch << b}

              until treeComplete or treeTooLarge

                # If K>0, density-dependent speciation is assumed, with speciation rate tmpLambda = [trialLambda * (1 - m/K) , 0].max when there are m extant species.
                if @k > 0
                  tmpLambda = [trialLambda * (1 - currentBranch.size/@k) , 0].max
                else
                  tmpLambda = trialLambda
                end

                # Draw a random time step length from an exponential distribution with mean 1/(-currentBranch.size*(tmpLambda+trialMu)).
                # This is faster than using Rubystats::ExponentialDistribution.new(currentBranch.size*(tmpLambda+trialMu)).
                currentTime = Math.log(1-rand)/(-currentBranch.size*(tmpLambda+trialMu))

                if currentTime > @present

                  # Select a branch from the current branches and remove it from current branches.
                  selectedBranch = currentBranch.delete_at(rand(currentBranch.size))

                  # Update the termination time of the selected branch.
                  selectedBranch.updateTermination(currentTime)
                  if rand < tmpLambda/(tmpLambda+trialMu) # The selected branch speciates.

                    # Update the parent branch with the Ids of the two newly created daughter branches, and update its endCause.
                    currentId = ["b#{branchIdCounter}","b#{branchIdCounter+1}"]
                    selectedBranch.updateDaughterId(currentId)
                    selectedBranch.updateEndCause("speciation")

                    # Produce two daugther branches, add them to current branches and adjust the branchIdCounter.
                    2.times do |d|
                      @branch << Branch.new(id = currentId[d], origin = currentTime, termination = nil, parentId = selectedBranch.id, daughterId = ["unborn","unborn"], endCause = nil)
                      currentBranch << @branch.last
                      branchIdCounter += 1
                    end

                  else
                    # The selected branch goes extinct. Update its endCause and it's daughter IDs.
                    selectedBranch.updateEndCause("extinction")
                    selectedBranch.updateDaughterId("none","none")                    
                  end

                else
                  # When the present is reached all current branches' termination dates are updated.
                  currentBranch.each {|b| b.updateTermination(@present-0.001)}
                  treeComplete = true
                end

                # The termination criterion is reached if either the tree is complete, or too large.
                treeComplete = true if currentBranch.size == 0
                treeTooLarge = true if @branch.size > 10 * @maxNp unless @maxNp == nil
              end

            else
              raise "Invalid algorithm specified"
            end
            
            # The program comes to this point if either treeComplete or treeTooLarge. If treeTooLarge, we'ld like to skip the remainder of this loop.
            unless treeTooLarge
              # Cut all branches that terminate in the future, and at the same time count the number of extant species
              @Np = 0
              @branch.size.times do |b|
                if @branch[b].termination <= @present
                  @branch[b].cutAtPresent(@present)
                  @Np += 1
                end
              end
              
              # Check whether the tree meets the given condition
              if @maxNp == nil
                if @Np >= @minNp
                  treeMeetsCondition = true
                end
              else
                if @Np >= @minNp and @Np <= @maxNp
                  treeMeetsCondition = true
                end
              end
            end

            # Because of the above, treeMeetsCondition can only be true unless treeTooLarge. Thus if treeTooLarge, the remainder will be skipped anyway.
            if treeMeetsCondition == true # tentatively

              # Let's make sure that the array @branch is still sorted by the id numbers. This is important for the quick parent search below.
              branchesSorted = true
              (@branch.size-1).times do |b|
                if @branch[b].id[1..-1].to_i > @branch[b+1].id[1..-1].to_i
                  branchesSorted = false
                  break
                end
              end
              if branchesSorted == false
                # This is unlikely to happen anyway, and should be checked if it does.
                warn "WARNING: For some reason, branches are not sorted according to their id. Better check out why."
                until branchesSorted
                  allSortedSoFar = true
                  (@branch.size-1).times do |b|
                    if @branch[b].id[1..-1].to_i > @branch[b+1].id[1..-1].to_i
                      @branch[b],@branch[b+1] = @branch[b+1],@branch[b]
                      allSortedSoFar = false
                    end
                  end
                  branchesSorted = true if allSortedSoFar == true
                end
              end

              # Determine the progeny of each branch (needed to know whether conditions are met, and for fossil constraints).
              # First of all, set progeniesComplete to 2 for all extinct and present branches.
              @branch.size.times do |b|
                @branch[b].progeniesComplete = 2 if @branch[b].endCause != "speciation"
              end
              allProgeniesComplete = false
              until allProgeniesComplete == true do
                newlyCompleted = []
                @branch.size.times do |b|
                  # If the progeny of this branch is clear but has not been passed on to the parent yet...
                  if @branch[b].progeniesComplete == 2
                    newlyCompleted << @branch[b] if @branch[b].progenyPassedOn == false
                  end
                end
                allProgeniesComplete = true if newlyCompleted == []
                newlyCompleted.size.times do |n|
                  # Find parent, pass progeny+self on to parents progeny, add parent.progeniesComplete += 1, and change own progenyPassedOn to true
                  lower = 0
                  upper = @branch.size
                  parentFound = false
                  until parentFound == true
                    # Use the fact that @branch is sorted according to id for a quick parent search.
                    b = lower+(upper-lower)/2
                    if @branch[b].id[1..-1].to_i < newlyCompleted[n].parentId[1..-1].to_i # Comparing the index numbers (e.g. '123') of branch ids (e.g. 'b123').
                      lower = b
                    elsif @branch[b].id[1..-1].to_i > newlyCompleted[n].parentId[1..-1].to_i
                      upper = b
                    elsif @branch[b].id == newlyCompleted[n].parentId
                      parentFound = true
                      ary = newlyCompleted[n].progenyId.dup
                      ary << newlyCompleted[n].id
                      aryFlat = ary.flatten
                      aryFlat.size.times do |a|
                        @branch[b].progenyId << aryFlat[a]
                      end
                      @branch[b].progeniesComplete += 1
                      next
                    elsif newlyCompleted[n].parentId == "treeOrigin"
                      parentFound = true
                    else
                      raise "Some problem occurred with the fast parent search!" # This should never be called at all.
                    end
                  end
                  newlyCompleted[n].progenyPassedOn = true
                end
              end

              # Determine the extant progeny of every branch (needed to know whether both b0 and b1 have extant progeny).
              # At the same time, determine the extant diversity of each branch.
              @branch.size.times do |b|
                if @branch[b].extant
                  @branch[b].updateExtantDiversity(1)
                  @branch[b].updateOriginalExtantDiversity(1)
                else
                  @branch[b].progenyId.size.times do |p|
                    # Find the branch, whose id matches @branch[b].progenyId[p]. As above, we use the fact that @branch is sorted according to id for a quick search.
                    idFound = false
                    lower = 0
                    upper = @branch.size
                    until idFound == true
                      bb = lower+(upper-lower)/2
                      if @branch[bb].id[1..-1].to_i < @branch[b].progenyId[p][1..-1].to_i
                        lower = bb
                      elsif @branch[bb].id[1..-1].to_i > @branch[b].progenyId[p][1..-1].to_i
                        upper = bb
                      elsif @branch[bb].id == @branch[b].progenyId[p]
                        idFound = true
                        if @branch[bb].extant
                          @branch[b].extantProgenyId << @branch[b].progenyId[p]
                          @branch[b].updateExtantDiversity(@branch[b].extantDiversity + 1)
                          @branch[b].updateOriginalExtantDiversity(@branch[b].originalExtantDiversity + 1)
                        end
                      else 
                        raise "Some problem occurred with the fast progeny search!" # This should never be called at all.
                      end
                    end
                  end
                end
              end

              # Check whether the extant progeny of b0 and b1 is according to the tree condition.
              if @rootSplit == true
                unless @minNpEach == nil
                  if @minNpEach > 0
                    treeMeetsCondition = false if @branch[0].extantProgenyId.size < @minNpEach or @branch[1].extantProgenyId.size < @minNpEach
                  end
                end
                unless @maxNpEach == nil
                  if @maxNpEach > 0
                    treeMeetsCondition = false if @branch[0].extantProgenyId.size > @maxNpEach or @branch[1].extantProgenyId.size > @maxNpEach
                  end
                end
              end
            end # if treeMeetsCondition == true # tentatively

            # If all is good, write it up, save the dump file.
            if treeMeetsCondition
              # Fix trial parameters.
              @trial = trial
              @posteriorLambda = trialLambda
              @posteriorMu = trialMu
              @posteriorTreeOrigin = trialTreeOrigin
# Change the below commented parts to allow multithreading.
# done = true
              # Write storage array to dump file.
             dumpFile = File.new(dumpFileNames[t],"w")
             storage = [@branch,@Np,@trial,@posteriorLambda,@posteriorMu,@posteriorTreeOrigin]
             dumpFile.write(Marshal.dump(storage))
             dumpFile.close
            end # if treeMeetsCondition

            exit # Only fork exit.
          end # Process.fork do
        end # @threads.times do |t|
        Process.waitall
        dumpFileNames.each do |d|
          done = true if File.exists?(d)
        end
        trial += @threads
      end # until done

      existingDumpFileNames = []
      dumpFileNames.each do |d|
        existingDumpFileNames << d if File.exist?(d)
      end
      # Retrieve the tree information that was written to file by one of the forks, and clean up the dump and trial files.
      storage = Marshal.load(File.read(existingDumpFileNames.sample))
      @branch = storage[0]
      @Np = storage[1]
      @trial = storage[2]
      @posteriorLambda = storage[3]
      @posteriorMu = storage[4]
      @posteriorTreeOrigin = storage[5]
      existingDumpFileNames.each do |e|
        File.delete(e)
      end

      # Assign branches to species and create new species
      @species = []
      speciesCount = 0
      @branch.each do |b|
        if b.speciesId == nil
          newSpeciesId = "s#{speciesCount}"
          speciesCount += 1
          b.addSpeciesId(newSpeciesId)
          @species << Species.new(newSpeciesId,b) # This guarantees that species are not only sorted by their id, but also that the id is the same as their index in the array species.
        else
          if @species[b.speciesId[1..-1].to_i].id == b.speciesId
            @species[b.speciesId[1..-1].to_i].addBranch(b)
          else
            warn "WARNING: For some reason, the id of species #{species[b.speciesId[1..-1].to_i].id} does not match it's position in the @species array."
            # Do a fast search for the species whose id matches b.speciesId, making use of the fact that species are sorted according to their id.
            lower = -1
            upper = @species.size
            spFound = false
            until spFound
              s = lower+(upper-lower)/2
              if @species[s].id[1..-1].to_i < b.speciesId[1..-1].to_i
                lower = s
              elsif @species[s].id[1..-1].to_i > b.speciesId[1..-1].to_i
                upper = s
              elsif @species[s].id == b.speciesId
                spFound = true
                @species[s].addBranch(b)
                next
              else
                raise "Some problem occurred with the fast species search!" # This should never be called at all.
              end
            end
          end
        end
        if b.endCause == "speciation"
          selectedDaughterId = b.daughterId.sample
          # Do a fast search for branch whose id matches the selectedDaughterId.
          lower = -1
          upper = @branch.size
          dsFound = false
          until dsFound
            bb = lower+(upper-lower)/2
            if @branch[bb].id[1..-1].to_i < selectedDaughterId[1..-1].to_i
              lower = bb
            elsif @branch[bb].id[1..-1].to_i > selectedDaughterId[1..-1].to_i
              upper = bb
            elsif @branch[bb].id == selectedDaughterId
              dsFound = true
              @branch[bb].addSpeciesId(b.speciesId)
              next
            else
              raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
            end
          end
        end
      end

      # Add to each species whether it is extant or not
      @species.each do |s|
        s.addExtant(@present)
      end

      # Calculate sumOfSpeciesDurations
      @sumOfSpeciesDurations = 0
      @species.each do |s|
        @sumOfSpeciesDurations += s.duration
      end
      
      # Store the full (unreconstructed) number of extant species. '@Np.dup' is not required as Fixnums are always duplicated.
      @NpFull = @Np

      if verbose
        endTime = Time.now
        print "\r"
        print "                                                                                                                                                  "
        print "\r"
        print "Tree found after #{@trial} trials."
        unless allFixed
          puts " Successful trial settings were lambda = #{@posteriorLambda.round(5)}, mu = #{@posteriorMu.round(5)}, and treeOrigin = #{@posteriorTreeOrigin.round(2)}."
        end
        if rootSplit
          puts "Tree has #{@Np} extant species, with #{@branch[0].extantProgenyId.size} and #{@branch[1].extantProgenyId.size} extant species in the two partial trees descending from the root."
        else
          puts "Tree has #{@Np} extant species."
        end
        if endTime-startTime < 60
          puts "Time used: #{(endTime-startTime).round(2)} seconds."
        elsif endTime-startTime < 3600
          puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
        else
          puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
        end
        puts
      end
    else
      if verbose
        startTime = Time.now
        puts
        puts "----------------------------------------------------phylsim.rb | Import-----------------------------------------------------"
        puts
      end

      # Initiate nTips (needed after the subsequent if loop).
      nTips = 0

      if @fileType == "dump"
        loadTree = Marshal.load(File.read(@fileName))
        @branch = loadTree.branch
        @species = loadTree.species
        @lambda = loadTree.lambda
        @mu = loadTree.mu
        @treeOrigin = loadTree.treeOrigin
        @present = loadTree.present
        @k = loadTree.k
        @rootSplit = loadTree.rootSplit
        @Np = loadTree.Np
        @sumOfSpeciesDurations = loadTree.sumOfSpeciesDurations
        @fossilRecordAdded = loadTree.fossilRecordAdded
        @sequencesEvolved = loadTree.sequencesEvolved
        @traitsEvolved = loadTree.traitsEvolved
        @treeReconstructed = loadTree.treeReconstructed
        @branchRatesAssigned = loadTree.branchRatesAssigned
        @alignmentPositionsAssigned = loadTree.alignmentPositionsAssigned
        @posteriorLambda = loadTree.posteriorLambda
        @posteriorMu = loadTree.posteriorMu
        @posteriorTreeOrigin = loadTree.posteriorTreeOrigin
        @phylogeneticDiversity = loadTree.phylogeneticDiversity if @treeReconstructed
        @gamma = loadTree.gamma if @treeReconstructed
        @NpFull = loadTree.NpFull
        @samplingScheme = loadTree.samplingScheme
        @focusGroupNpFull = loadTree.focusGroupNpFull
        @focusGroupAge = loadTree.focusGroupAge
        @ecmFileName = loadTree.ecmFileName
        @substitutionModel = loadTree.substitutionModel
        @branchRatesMean = loadTree.branchRatesMean
        @branchRatesStandardDeviation = loadTree.branchRatesStandardDeviation
        @branchRatesAutoCorrelation = loadTree.branchRatesAutoCorrelation

      elsif @fileType == "nexus" or @fileType == "newick"
        if @fileType == "nexus"
          nexusFile = File.open(@fileName)
          nexusLines = nexusFile.readlines
          unless nexusLines[0].strip.downcase == "#nexus"
            raise "File #{@fileName} does not seem to be in nexus format!"
          end
          treeLines = []
          inTreeBlock = false
          pastTreeBlock = false
          nexusLines.size.times do |x|
            if nexusLines[x].strip.downcase == "begin trees;" and inTreeBlock == false and pastTreeBlock == false
              inTreeBlock = true
              pastTreeBlock = false
            elsif inTreeBlock and nexusLines[x].strip != "" and nexusLines[x].strip.downcase != "end;"
              treeLines << nexusLines[x]
            elsif inTreeBlock and nexusLines[x].strip.downcase == "end;"
              inTreeBlock = false
              pastTreeBlock = true
            end
          end
          treeLines.each do |tl|
            if tl.downcase.include?("translate")
              raise "\"Translate\" blocks in the Nexus file are not supported yet. Please open the tree file in Figtree and export tree again with setting \"Save as currently displayed\"."
            end
          end
          treeLine = treeLines[treeNumber]
          unless treeLine.include?("=")
            raise "File #{@fileName} does not seem to be a valid nexus file!"
          end
          treeLine.strip!.chomp!(";")
          cleanTreeString = ""
          inTreePart = false
          inComment = false
          treeLine.size.times do |x|
            if x > 0
              inTreePart = true if treeLine[x-1] == "="
            end
            if inTreePart
              if x > 0
                inComment = false if treeLine[x-1] == "]"
              end
              inComment = true if treeLine[x] == "["
              cleanTreeString << treeLine[x] unless inComment
            end
          end
        elsif @fileType == "newick"
          newickFile = File.open(@fileName,"r")
          newickLines = newickFile.readlines
          treeLine = newickLines[treeNumber]
          unless treeLine.include?("(")
            raise "File #{@fileName} does not seem to be a vaild newick file!"
          end
          if treeLine.include?("[")
            raise "Comments are not supported in newick trees!"
          end
          cleanTreeString = treeLine.strip.chomp(";")
        end
        cleanTreeString.strip!
        if cleanTreeString == ""
          raise "The tree in file #{@fileName} could not be read!"
        end
        if cleanTreeString.count("(") != cleanTreeString.count(")")
          raise "The tree in file #{@fileName} contains unequal numbers of opening and closing brackets!"
        end
        # Replace scientific notation with decimal notation.
        while cleanTreeString.match(/:(\d\.\d*[eE]\-\d+)[,\)]/)
          scientific = $1
          decimal = ("%f" % $1).sub(/\.?0*$/, "")
          puts decimal
          cleanTreeString.sub!(scientific,decimal)
        end

        # Parse the tree string and store all information temporarily in arrays.
        tmpBranchEndNodeId = []
        tmpBranchDuration = []
        tmpBranchStartNodeId = []
        numberOfInternalNodes = 0
        while cleanTreeString.match(/\(([a-zA-Z0-9_]+?):([\d\.]+?)\,([a-zA-Z0-9_]+?):([\d\.]+?)\)/)
          numberOfInternalNodes += 1
          tmpBranchEndNodeId << $1
          tmpBranchStartNodeId << "internalNode#{numberOfInternalNodes}X"
          tmpBranchDuration << $2.to_f
          tmpBranchEndNodeId << $3
          tmpBranchDuration << $4.to_f
          tmpBranchStartNodeId << "internalNode#{numberOfInternalNodes}X"
          cleanTreeString.sub!("(#{$1}:#{$2},#{$3}:#{$4})","internalNode#{numberOfInternalNodes}X")
        end
        # If we had polytomies, cleanTreeString is not fully resolved yet. Make sure this is not the case.
        if cleanTreeString.include?(",")
          raise "Apparently the parsed tree contained polytomies. These are not supported."
        end

        # Find the maximum duration between tip and root.
        rootId = tmpBranchStartNodeId[-1]
        totalLengthsToRoot = []
        tmpBranchEndNodeId.size.times do |x|
          unless tmpBranchEndNodeId[x].match(/internalNode\d+X/)
            totalLengthToRoot = tmpBranchDuration[x]
            rootReached = false
            currentIndex = 0
            currentStartNodeId = tmpBranchStartNodeId[x]
            unless currentStartNodeId == rootId
              while rootReached == false
                tmpBranchEndNodeId.size.times do |y|
                  if tmpBranchEndNodeId[y] == currentStartNodeId
                    currentIndex = y
                    currentStartNodeId = tmpBranchStartNodeId[y]
                    totalLengthToRoot += tmpBranchDuration[y]
                    rootReached = true if currentStartNodeId == rootId
                    break
                  end
                end
              end
            end
            totalLengthsToRoot << totalLengthToRoot
          end
        end
        rootAge = totalLengthsToRoot.max.round(5)

        # Prepare arrays for temporary branch format.
        tmp2BranchID = []
        tmp2BranchOrigin = []
        tmp2BranchTermination = []
        tmp2BranchParentId = []
        tmp2BranchDaughterId1 = []
        tmp2BranchDaughterId2 = []
        tmp2BranchEndCause = []
        tmp2BranchEndNodeId = []

        # Prepare the first two branches in temporary format (tmpBranchEndNodeId[-1] and tmpBranchEndNodeId[-2] are the two oldest branches).
        # If the first root branch ends in an internal node...
        if tmpBranchEndNodeId[-1].match(/internalNode\d+X/)
          tmp2BranchID << "b0"
          tmp2BranchOrigin << rootAge
          tmp2BranchTermination << rootAge-tmpBranchDuration[-1]
          tmp2BranchParentId << "treeOrigin"
          tmp2BranchDaughterId1 << "unborn"
          tmp2BranchDaughterId2 << "unborn"
          tmp2BranchEndCause << "speciation"
          tmp2BranchEndNodeId << tmpBranchEndNodeId[-1]
        # If it doesn't discriminate between two cases...
        else
          # If it's duration is the same as the root age, it extends all the way to the present.
          if rootAge - tmpBranchDuration[-1] <= 0.001
            tmp2BranchID << "b0"
            tmp2BranchOrigin << rootAge
            tmp2BranchTermination << 0.0
            tmp2BranchParentId << "treeOrigin"
            tmp2BranchDaughterId1 << "none"
            tmp2BranchDaughterId2 << "none"
            tmp2BranchEndCause << "present"
            tmp2BranchEndNodeId << tmpBranchEndNodeId[-1]
          # If it doesn't, it went extinct.
          else
            tmp2BranchID << "b0"
            tmp2BranchOrigin << rootAge
            tmp2BranchTermination << rootAge-tmpBranchDuration[-1]
            tmp2BranchParentId << "treeOrigin"
            tmp2BranchDaughterId1 << "none"
            tmp2BranchDaughterId2 << "none"
            tmp2BranchEndCause << "extinction"
            tmp2BranchEndNodeId << tmpBranchEndNodeId[-1]
          end
        end
        # Repeat the above for the second branch.
        # If the second root branch ends in an internal node...
        if tmpBranchEndNodeId[-2].match(/internalNode\d+X/)
          tmp2BranchID << "b1"
          tmp2BranchOrigin << rootAge
          tmp2BranchTermination << (rootAge-tmpBranchDuration[-2]).round(5)
          tmp2BranchParentId << "treeOrigin"
          tmp2BranchDaughterId1 << "unborn"
          tmp2BranchDaughterId2 << "unborn"
          tmp2BranchEndCause << "speciation"
          tmp2BranchEndNodeId << tmpBranchEndNodeId[-2]
        # If it doesn't discriminate between two cases...
        else
          # If it's duration is the same as the root age, it extends all the way to the present.
          if rootAge - tmpBranchDuration[-2] <= 0.001
            tmp2BranchID << "b1"
            tmp2BranchOrigin << rootAge
            tmp2BranchTermination << 0.0
            tmp2BranchParentId << "treeOrigin"
            tmp2BranchDaughterId1 << "none"
            tmp2BranchDaughterId2 << "none"
            tmp2BranchEndCause << "present"
            tmp2BranchEndNodeId << tmpBranchEndNodeId[-2]
          # If it doesn't, it went extinct.
          else
            tmp2BranchID << "b1"
            tmp2BranchOrigin << rootAge
            tmp2BranchTermination << (rootAge-tmpBranchDuration[-2]).round(5)
            tmp2BranchParentId << "treeOrigin"
            tmp2BranchDaughterId1 << "none"
            tmp2BranchDaughterId2 << "none"
            tmp2BranchEndCause << "extinction"
            tmp2BranchEndNodeId << tmpBranchEndNodeId[-2]
          end
        end

        # Find out about all remaining branches until either all branches end with extinctions, or all branches have reached the present.
        branchIdCounter = 2
        treeComplete = false
        until treeComplete
          change = false
          tmp2BranchID.size.times do |x| # go through all branches
            if tmp2BranchDaughterId1[x] == "unborn" and tmp2BranchDaughterId2[x] == "unborn" # if a branch terminated with a speciation event in the past, then add the two daughter branches
              # Find the two branches that have the same start node as this branch's end node.
              tmpBranchStartNodeId.size.times do |y|
                if tmpBranchStartNodeId[y] == tmp2BranchEndNodeId[x]
                  if tmpBranchEndNodeId[y].match(/internalNode\d+X/)
                    tmp2BranchID << "b#{branchIdCounter}"
                    tmp2BranchOrigin << tmp2BranchTermination[x]
                    tmp2BranchTermination << (tmp2BranchTermination[x]-tmpBranchDuration[y]).round(5)
                    tmp2BranchParentId << tmp2BranchID[x]
                    tmp2BranchDaughterId1 << "unborn"
                    tmp2BranchDaughterId2 << "unborn"
                    tmp2BranchEndCause << "speciation"
                    tmp2BranchEndNodeId << tmpBranchEndNodeId[y]
                  else
                    if tmp2BranchTermination[x] - tmpBranchDuration[y] <= 0.001
                      tmp2BranchID << "b#{branchIdCounter}"
                      tmp2BranchOrigin << tmp2BranchTermination[x]
                      tmp2BranchTermination << 0.0
                      tmp2BranchParentId << tmp2BranchID[x]
                      tmp2BranchDaughterId1 << "none"
                      tmp2BranchDaughterId2 << "none"
                      tmp2BranchEndCause << "present"
                      tmp2BranchEndNodeId << tmpBranchEndNodeId[y]
                    else
                      tmp2BranchID << "b#{branchIdCounter}"
                      tmp2BranchOrigin << tmp2BranchTermination[x]
                      tmp2BranchTermination << (tmp2BranchTermination[x]-tmpBranchDuration[y]).round(5)
                      tmp2BranchParentId << tmp2BranchID[x]
                      tmp2BranchDaughterId1 << "none"
                      tmp2BranchDaughterId2 << "none"
                      tmp2BranchEndCause << "extinction"
                      tmp2BranchEndNodeId << tmpBranchEndNodeId[y]
                    end
                  end
                  # Update daughter ids of temporary parent.
                  if tmp2BranchDaughterId1[x] == "unborn"
                    tmp2BranchDaughterId1[x] = "b#{branchIdCounter}"
                  else
                    tmp2BranchDaughterId2[x] = "b#{branchIdCounter}"
                  end

                  # Increase the branchIdCounter
                  branchIdCounter += 1
                  change = true
                end
              end

            end # if tmp2BranchDaughterId1[x] == "unborn"

          end # tmp2BranchID.size.times do |x|
          treeComplete = true if change == false
        end # until treeComplete

        # Fill array @branch, and at the same time, add species for terminal branches.
        @branch = []
        @species = []
        tmp2BranchID.size.times do |x|
          @branch << Branch.new(tmp2BranchID[x], tmp2BranchOrigin[x], tmp2BranchTermination[x], tmp2BranchParentId[x], [tmp2BranchDaughterId1[x],tmp2BranchDaughterId2[x]], tmp2BranchEndCause[x])
          if tmp2BranchEndNodeId[x].match(/internalNode\d+X/)
            newSpeciesId = "unknown"
            @branch.last.addSpeciesId(newSpeciesId)
          else
            newSpeciesId = tmp2BranchEndNodeId[x]
            @branch.last.addSpeciesId(newSpeciesId)
            @species << Species.new(newSpeciesId,@branch.last)
          end
        end
        @branch.each do |b|
          if b.endCause == "present"
            b.updateExtant(true)
          else
            b.updateExtant(false)
          end
        end

        # Determine the progeny of each branch (needed to know whether conditions are met, and for fossil constraints).
        # First of all, set progeniesComplete to 2 for all extinct and present branches.
        @branch.size.times do |b|
          @branch[b].progeniesComplete = 2 if @branch[b].endCause != "speciation"
        end
        allProgeniesComplete = false
        until allProgeniesComplete == true do
          newlyCompleted = []
          @branch.size.times do |b|
            # If the progeny of this branch is clear but has not been passed on to the parent yet...
            if @branch[b].progeniesComplete == 2
              newlyCompleted << @branch[b] if @branch[b].progenyPassedOn == false
            end
          end
          allProgeniesComplete = true if newlyCompleted == []
          newlyCompleted.size.times do |n|
            # Find parent, pass progeny+self on to parents progeny, add parent.progeniesComplete += 1, and change own progenyPassedOn to true
            lower = 0
            upper = @branch.size
            parentFound = false
            until parentFound == true
              # Use the fact that @branch is sorted according to id for a quick parent search.
              b = lower+(upper-lower)/2
              if @branch[b].id[1..-1].to_i < newlyCompleted[n].parentId[1..-1].to_i # Comparing the index numbers (e.g. '123') of branch ids (e.g. 'b123').
                lower = b
              elsif @branch[b].id[1..-1].to_i > newlyCompleted[n].parentId[1..-1].to_i
                upper = b
              elsif @branch[b].id == newlyCompleted[n].parentId
                parentFound = true
                ary = newlyCompleted[n].progenyId.dup
                ary << newlyCompleted[n].id
                aryFlat = ary.flatten
                aryFlat.size.times do |a|
                  @branch[b].progenyId << aryFlat[a]
                end
                @branch[b].progeniesComplete += 1
                next
              elsif newlyCompleted[n].parentId == "treeOrigin"
                parentFound = true
              else
                raise "Some problem occurred with the fast parent search!" # This should never be called at all.
              end
            end
            newlyCompleted[n].progenyPassedOn = true
          end
        end

        # Determine the extant progeny of every branch (needed to know whether both b0 and b1 have extant progeny).
        @branch.size.times do |b|
          @branch[b].progenyId.size.times do |p|
            # Find the branch, whose id matches @branch[b].progenyId[p]. As above, we use the fact that @branch is sorted according to id for a quick search.
            idFound = false
            lower = 0
            upper = @branch.size
            until idFound == true
              bb = lower+(upper-lower)/2
              if @branch[bb].id[1..-1].to_i < @branch[b].progenyId[p][1..-1].to_i
                lower = bb
              elsif @branch[bb].id[1..-1].to_i > @branch[b].progenyId[p][1..-1].to_i
                upper = bb
              elsif @branch[bb].id == @branch[b].progenyId[p]
                idFound = true
                @branch[b].extantProgenyId << @branch[b].progenyId[p] if @branch[bb].extant
              else 
                raise "Some problem occurred with the fast progeny search!" # This should never be called at all.
              end
            end
          end
        end

        # Transfer extantProgenyId to originalExtantProgenyId (using method updateExtantProgenyId() as this automatically does the transfer).
        @branch.each do |b|
          b.updateExtantProgenyId(b.extantProgenyId)
        end

        # Add to each species whether it is exant or not.
        @species.each do |s|
          s.addExtant(@present)
        end

        # Deal with the diversity file.
        unless diversityFileName == nil
          diversityFile = File.open(diversityFileName)
          diversityFileLines = diversityFile.readlines
          diversityFileSpeciesIDs = []
          diversityFileDiversities = []
          diversityFileLines.each do |l|
            diversityFileLineAry = l.split("\t")
            diversityFileSpeciesIDs << diversityFileLineAry[0]
            diversityFileDiversities << diversityFileLineAry[1].to_i
          end
          # Make sure all diversities are greater than 1.
          diversityFileDiversities.each{|d| raise "Found a diversity less than 1: #{d}!" if d < 1}
          # Make sure that species IDs in the file occur only once.
          duplicateSpeciesIDsInFile = diversityFileSpeciesIDs - diversityFileSpeciesIDs.uniq
          if duplicateSpeciesIDsInFile.size > 0
            raiseString = "The following species IDs have been found multiple times in file #{diversityFileName}:"
            raiseString.each{|d| raiseString << "#{d}, "}
            raiseString.chomp(", ")
            raiseString << "!"
            raise raiseString
          end
          # Make sure all species IDs listed in array "species" are found in the file and vice versa.
          idsMissingInFile = []
          @species.each do |s|
            unless diversityFileSpeciesIDs.include?(s.id)
              idsMissingInFile << s.id
            end
          end
          if idsMissingInFile.size > 0
            raiseString = "The following species IDs have not been found in file #{diversityFileName}:"
            idsMissingInFile.each{|d| raiseString << "#{d}, "}
            raiseString.chomp(", ")
            raiseString << "!"
            raise raiseString
          end
          idsMissingInArray = []
          diversityFileSpeciesIDs.each do |d|
            foundInSpeciesIDs = false
            @species.each{|s| foundInSpeciesIDs = true if s.id == d}
            idsMissingInArray << d if foundInSpeciesIDs == false
          end
          if idsMissingInArray.size > 0
            raiseString = "The following species IDs have been found in file #{diversityFileName}, but not in the species array:"
            idsMissingInArray.each{|d| raiseString << "#{d}, "}
            raiseString.chomp(", ")
            raiseString << "!"
            raise raiseString
          end
        end

        # Calculate the extant diversity for each branch.
        if diversityFileName == nil
          @branch.each do |b|
            if b.extant
              b.updateExtantDiversity(1)
              b.updateOriginalExtantDiversity(1)
            else
              b.updateExtantDiversity(b.extantProgenyId.size)
              b.updateOriginalExtantDiversity(b.extantProgenyId.size)
            end
          end
        else
          @branch.each do |b|
            if b.extant
              b.updateExtantDiversity(1)
              originalExtantDiversity = diversityFileDiversities[diversityFileSpeciesIDs.index(b.speciesId)]
              raise "The original extant diversity for species #{b.speciesId} was not found in file #{diversityFileName}" if originalExtantDiversity < 1
              b.updateOriginalExtantDiversity(originalExtantDiversity)
            else
              b.updateExtantDiversity(b.extantProgenyId.size)
            end
          end
          @branch.each do |b|
            # Find the entire ancestry of this branch if it's extant, and increase the extant diversity of each ancestor by that of the current branch.
            if b.extant
              ancestors = [b]
              ancestryComplete = false
              until ancestryComplete
                if ancestors.last.parentId == "treeOrigin"
                  ancestryComplete = true
                else
                  @branch.each do |bb|
                    if bb.id == ancestors.last.parentId
                      ancestors << bb
                      break
                    end
                  end
                end
              end
              # Increase the original extant diversity of each ancestor.
              ancestors[1..-1].each do |a|
                a.updateOriginalExtantDiversity(a.originalExtantDiversity + b.originalExtantDiversity)
              end
            end
          end
        end

        # Calculate sumOfSpeciesDurations.
        @sumOfSpeciesDurations = 0
        @branch.each do |b|
          @sumOfSpeciesDurations += b.duration
        end
        
        # If the tree was parsed, @Np is still nil and should be recalculated.
        if @Np == nil
          @Np = 0
          @branch.each{|b| @Np += 1 if b.extant}
        end

        # Store the full (unreconstructed) number of extant species. '@Np.dup' is not required as Fixnums are always duplicated.
        @NpFull = @Np

        # Store the tree origin.
        @treeOrigin = @branch[0].origin

      else
        raise "Unrecognized file type: #{@fileType}"
      end

      if verbose
        endTime = Time.now
        # Calculate the number of tips, regardless of whether they're extant or not.
        @branch.each{|b| nTips += 1 if b.endCause != "speciation"}
        if nTips == @Np
          puts "Tree successfully loaded from '#{@fileName}'. The tree contains #{nTips} extant tips."
        else
          puts "Tree successfully loaded from '#{@fileName}'. The tree contains #{nTips} tips, of which #{@NpFull} are extant."
        end
        if endTime-startTime < 60
          puts "Time used: #{(endTime-startTime).round(2)} seconds."
        elsif endTime-startTime < 3600
          puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
        else
          puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
        end
        puts
      end
    end
  end
  
  def kendallMoran(intervalStart, intervalEnd)
    if intervalStart < intervalEnd
      intervalStart,intervalEnd = intervalEnd,intervalStart
    end
    if intervalStart > @treeOrigin
      puts "Interval start is older than the tree origin"
    elsif intervalEnd < 0
      puts "Interval end is younger than the present"
    else
      speciesNumberAtIntervalStart = 0
      speciesNumberAtIntervalEnd = 0
      sumOfSpeciesDurationsInInterval = 0.0
      @branch.each do |b|
        if b.origin > intervalStart
          if b.termination <= intervalStart and b.termination > intervalEnd
            speciesNumberAtIntervalStart += 1
            sumOfSpeciesDurationsInInterval += intervalStart-b.termination
          elsif b.termination <= intervalEnd
            speciesNumberAtIntervalStart += 1
            speciesNumberAtIntervalEnd += 1
            sumOfSpeciesDurationsInInterval += intervalStart-intervalEnd
          end
        elsif b.origin <= intervalStart and b.origin > intervalEnd
          if b.termination <= intervalStart and b.termination > intervalEnd
            sumOfSpeciesDurationsInInterval += b.origin-b.termination
          elsif b.termination <= intervalEnd
            speciesNumberAtIntervalEnd += 1
            sumOfSpeciesDurationsInInterval += b.origin-intervalEnd
          end
        end
      end
      if speciesNumberAtIntervalEnd < speciesNumberAtIntervalStart
        raise "Species number is smaller at the end of the interval (#{speciesNumberAtIntervalEnd};#{intervalEnd}) than at the beginning (#{speciesNumberAtIntervalStart};#{intervalStart})!"
      end
      kendellmoran = (speciesNumberAtIntervalEnd-speciesNumberAtIntervalStart)/sumOfSpeciesDurationsInInterval
      kendellmoran
    end
  end

  def addFossilRecord(preservationRate, samplingGap = 0, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "----------------------------------------------------phylsim.rb | Fossil record -----------------------------------------------------"
      puts 
    end
    fossilWaitingTimeDist = Rubystats::ExponentialDistribution.new(preservationRate)
    numberOfFossils = 0
    @branch.size.times do |b|
      if samplingGap.class == Fixnum or samplingGap.class == Float
        sampling_gap_for_this_branch = samplingGap.to_f
      elsif samplingGap.class == Array and samplingGap.size == 2
        sampling_gap_for_this_branch_min = samplingGap[0]
        sampling_gap_for_this_branch_max = samplingGap[1]
        sampling_gap_for_this_branch = sampling_gap_for_this_branch_min + rand*(sampling_gap_for_this_branch_max-sampling_gap_for_this_branch_min)
      end
      fossil = []
      sumOfFossilWaitingTimes = 0
      branchTerminationReached = false
      until branchTerminationReached == true
        fossilWaitingTime = fossilWaitingTimeDist.rng
        sumOfFossilWaitingTimes += fossilWaitingTime
        if sumOfFossilWaitingTimes > @branch[b].duration
          branchTerminationReached = true
        else
          fossil << Fossil.new(branchId = @branch[b].id, age = @branch[b].origin-sumOfFossilWaitingTimes) if sumOfFossilWaitingTimes > sampling_gap_for_this_branch
          numberOfFossils += 1
        end
      end
      @branch[b].addFossils(fossil)
    end
    @fossilRecordAdded = true

    if verbose
      endTime = Time.now
      puts "\r#{numberOfFossils} fossils added with preservation rate #{preservationRate}.                 "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end

  end

  def branchReport(fileName = nil, overwrite = false, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "----------------------------------------------------phylsim.rb | Export-----------------------------------------------------"
      puts 
    end
    output = ""
    output << "Number of branches: #{@branch.size}\n"
    output << "\n"
    output << "Branches:\n"
    output << "\n"
    @branch.each do |b|
      output << "#{b.to_s}\n"
    end
    if fileName == nil
      output
    else
      unless overwrite
        # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
        if File.exists?(fileName)
          puts "Ok to replace '#{fileName}'? (y/N)"
          answer = gets
          unless answer[0].downcase == "y"
            puts "Please specify a new filename:"
            fileName = gets.strip
          end
          puts
        end
      end
      print "Reporting branches..." if verbose
      File.new(fileName,"w").write(output)
    end

    if verbose
      endTime = Time.now
      puts "\rInformation on #{@branch.size} branches reported to file #{fileName}.                 "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end

  end

  def traitReport(fileName = nil, overwrite = true, reportSpeciesTraits=true, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "-------------------------------------------------phylsim.rb | Trait report--------------------------------------------------"
      puts 
    end

    raise "Traits must first be evolved (method evolveTraits) before they can be reported!" unless @traitsEvolved
    output = []
    if reportSpeciesTraits
      speciesNumber = 0
      @branch.each{|b| speciesNumber += 1 if b.extant}
      output << "Number of species: #{speciesNumber}"
      output << ""
      output << "Species:"
    else
      output << "Number of branches: #{@branch.size}"
      output << ""
      output << "Branches:"
    end
    output << ""
    @branch.each do |b|
      if reportSpeciesTraits
        output <<  "#{b.speciesId.split("/")[-1]}\t#{b.endTrait}" if b.extant
      else
        output <<  "ID:            #{b.id}"
        output <<  "Start trait:   #{b.startTrait}"
        output <<  "End trait:     #{b.endTrait}"
        output << ""
      end
    end
    if fileName == nil
      output
    else
      unless overwrite
        # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
        if File.exists?(fileName)
          puts "Ok to replace '#{fileName}'? (y/N)"
          answer = gets
          unless answer[0].downcase == "y"
            puts "Please specify a new filename:"
            fileName = gets.strip
          end
          puts
        end
      end
      File.new(fileName,"w").puts output

      if verbose
        endTime = Time.now
        puts "\rTrait information written to '#{fileName}'.                 "
        if endTime-startTime < 60
          puts "Time used: #{(endTime-startTime).round(2)} seconds."
        elsif endTime-startTime < 3600
          puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
        else
          puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
        end
        puts
      end

    end

  end

  def speciesReport(extantOnly=false)
    output = []
    output << "Number of species: #{@species.size}"
    output << ""
    output << "Species:"
    output << ""
    if extantOnly
      @species.size.times do |s|
        if @species[s].extant
          output << "ID:            #{@species[s].id}"
          output << "Origin:        #{@species[s].origin.to_f.round(3)}"
          output << "Termination:   #{@species[s].termination.to_f.round(3)}"
          output << ""
        end
      end
    else
      @species.size.times do |s|
        output << "ID:            #{@species[s].id}"
        output << "Origin:        #{@species[s].origin.to_f.round(3)}"
        output << "Termination:   #{@species[s].termination.to_f.round(3)}"
        output << "Extant:        #{@species[s].extant}"
        output << ""
      end
    end
    output
  end

  def fossilReport(fileName = nil, overwrite = false, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "----------------------------------------------------phylsim.rb | Export-----------------------------------------------------"
      puts 
    end
    raise "In order to report fossil, these must first be assigned to all branches (method addFossilRecord)." unless @fossilRecordAdded
    output = []
    globalFossil = []
    @branch.size.times do |b|
      globalFossil << @branch[b].fossil
    end
    globalFossil.flatten!
    output << "Number of fossils: #{globalFossil.size}"
    output << ""
    globalFossil.size.times do |f|
      output << "Branch ID:         #{globalFossil[f].branchId}"
      output << "Age:               #{globalFossil[f].age.round(3)}"
      output << ""
    end
    if fileName == nil
      output
    else
      unless overwrite
        # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
        if File.exists?(fileName)
          puts "Ok to replace '#{fileName}'? (y/N)"
          answer = gets
          unless answer[0].downcase == "y"
            puts "Please specify a new filename:"
            fileName = gets.strip
          end
          puts
        end
      end
      if verbose
        puts "Reporting fossils..."
      end
      File.new(fileName,"w").puts output    
    end
  end

  def diversityReport(fileName = nil, bammFormat = false, overwrite = true, verbose = true)

    # Feedback.
    if verbose
      startTime = Time.now
      puts
      puts "-----------------------------------------------phylsim.rb | Diversity report------------------------------------------------"
      puts 
    end

    # Prepare the output array.
    output = []
    # if BAMM format is chosen, clade specific sampling probabilities will be reported according to http://bamm-project.org/advanced.html#incompsampling
    if bammFormat
      output << "1.0"
      @branch.each do |b|
        if b.extant
          output << "#{b.speciesId.split("/").last}\t#{b.speciesId.split("/").last}_clade\t#{(1/b.originalExtantDiversity.to_f).round(5)}"
        end
      end
    else
      @branch.each do |b|
        if b.extant
          output << "#{b.speciesId.split("/").last}\t#{b.originalExtantDiversity}"
        end
      end
    end

    # Return or write the array to file.
    if fileName == nil
      output
    else
      unless overwrite
        # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
        if File.exists?(fileName)
          puts "Ok to replace '#{fileName}'? (y/N)"
          answer = gets
          unless answer[0].downcase == "y"
            puts "Please specify a new filename:"
            fileName = gets.strip
          end
          puts
        end
      end
      File.new(fileName,"w").puts output    
    end

    # Feedback.
    if verbose
      endTime = Time.now
      puts "\Diversity information written to '#{fileName}'.                 "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end

  end # def diversityReport(fileName = nil, overwrite = false, verbose = true)

  def to_newick(fileName = nil, branchLengths = "duration", labels = true, plain = false, includeEmpty = true, overwrite = false, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "-------------------------------------------------phylsim.rb | Tree export---------------------------------------------------"
      puts 
    end

    # Check the branchLengths argument.
    raise "Either 'duration', 'expectedNumberOfSubstitutionsPerSite', 'actualNumberOfSubstitutions', or 'actualNumberOfSubstitutionsPerSite' must be specified for parameter 'branchLengths'!" unless ["duration","expectednumberofsubstitutionspersite","actualnumberofsubstitutions","actualnumberofsubstitutionspersite"].include?(branchLengths.downcase)
    raise "In order to export the tree with branch lengths according to the expected number of substitutions per site, branch rates must first be assigned!" if branchLengths.downcase == "expectednumberofsubstitutionspersite" and @branchRatesAssigned == false
    raise "In order to export the tree with branch lengths according to the actual number of substitutions, sequences must first be evolved!" if branchLengths.downcase == "actualnumberofsubstitutions" and @sequencesEvolved != true
    raise "In order to export the tree with branch lengths according to the actual number of substitutions per site, sequences must first be evolved!" if branchLengths.downcase == "actualnumberofsubstitutionspersite" and @sequencesEvolved != true

    # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
    if overwrite == false and fileName != nil
      if File.exists?(fileName)
        puts "Ok to replace '#{fileName}'? (y/N)"
        answer = gets
        unless answer[0].downcase == "y"
          puts "Please specify a new filename:"
          fileName = gets.strip
        end
        puts
      end
    end

    # Prepare the export.
    print "Preparing tree string in Newick format..." if verbose
    decimals = 6
    exportBranch = Marshal.load(Marshal.dump(@branch))

    # If taxa with empty sequences are to be removed, delete them from exportBranch using the same algorithm as for the reconstruction.
    removeStem = false
    unless includeEmpty
      exportBranch.each {|b| b.updateEndCause("extinction") if b.endSeq.gsub("-","") == "" and b.extant}

      # Branches with endCause 'extinction' are deleted, and if they're the only daughter of the parent, then the parent's endCause also becomes 'extinction'.
      # This is done until no branches have endCause 'extinction' anymore.
      extinction = true
      while extinction == true
        extinction = false
        exportBranch.size.times do |b|
          if exportBranch[b].endCause == "extinction"
            # We have to consider the following case. For all other reconstructions (as in the method 'reconstruct'), both sides of the root must have extant descendants.
            # This may not be the case here, and the endCause of one of the root branches itself could become 'extinction'. If this is the case, we will also have to
            # remove branches at the base of the other side of the root.
            if exportBranch[b].parentId == "treeOrigin"
              removeStem = true
            else
              # Fast parent search to find the branch that has the id that matches exportBranch[b].parentId. This is complicated by the fact that a branch can now be nil.
              lower = 0
              upper = exportBranch.size
              parentFound = false
              until parentFound == true
                # Use the fact that exportBranch is sorted according to id for a quick parent search.
                bb = lower+rand(upper-lower)
                unless exportBranch[bb] == nil
                  if exportBranch[bb].id[1..-1].to_i < exportBranch[b].parentId[1..-1].to_i
                    lower = bb
                  elsif exportBranch[bb].id[1..-1].to_i > exportBranch[b].parentId[1..-1].to_i
                    upper = bb
                  elsif exportBranch[bb].id == exportBranch[b].parentId
                    parentFound = true
                    if exportBranch[bb].daughterId.size == 2
                      daughterId = exportBranch[bb].daughterId.dup
                      daughterId.delete(exportBranch[b].id)
                      exportBranch[bb].updateDaughterId(daughterId)
                    elsif exportBranch[bb].daughterId.size == 1
                      exportBranch[bb].updateDaughterId([])
                      exportBranch[bb].updateEndCause("extinction")
                    else
                      raise "branch[#{bb}] has #{exportBranch[bb].daughterId.size} daughters!" # should never be called at all.
                    end
                    next
                  else
                    raise "Some problem occurred with the fast parent search: exportBranch[bb]: #{exportBranch[bb].id}, exportBranch[b].parentId: #{exportBranch[b].parentId}!" # This should never be called at all.
                  end
                end
              end
            end
            extinction = true
            exportBranch[b] = nil
          end
        end
        exportBranch.compact!
      end

      # Delete degree one vertices.
      unmerged = true
      while unmerged == true
        unmerged = false
        exportBranch.size.times do |b|
          unless exportBranch[b] == nil
            if exportBranch[b].daughterId.size == 1
              lower = 0
              upper = exportBranch.size
              daughterFound = false
              until daughterFound
                # Use the fact that exportBranch is sorted according to id for a quick daughter search.
                bb = lower+rand(upper-lower)
                unless exportBranch[bb] == nil
                  if exportBranch[bb].id[1..-1].to_i < exportBranch[b].daughterId[0][1..-1].to_i
                    lower = bb
                  elsif exportBranch[bb].id[1..-1].to_i > exportBranch[b].daughterId[0][1..-1].to_i
                    upper = bb
                  elsif exportBranch[bb].id == exportBranch[b].daughterId[0]
                    daughterFound = true
                    exportBranch[b].updateTermination(exportBranch[bb].termination)
                    exportBranch[b].updateDaughterId(exportBranch[bb].daughterId)
                    exportBranch[b].updateEndCause(exportBranch[bb].endCause)
                    exportBranch[b].updateExtant(exportBranch[bb].extant)
                    exportBranch[b].updateRate((exportBranch[b].duration*exportBranch[b].rate+exportBranch[bb].duration*exportBranch[bb].rate)/(exportBranch[b].duration+exportBranch[bb].duration)) if @branchRatesAssigned
                    exportBranch[b].updateDuration(exportBranch[b].origin - exportBranch[b].termination)
                    speciesIdAry = exportBranch[b].speciesId.split("/")
                    speciesIdAry1 = exportBranch[bb].speciesId.split("/")
                    speciesIdAry1.each do |s|
                      speciesIdAry << s unless speciesIdAry.include?(s)
                    end
                    speciesIdStr = ""
                    (speciesIdAry.size-1).times do |sid|
                      speciesIdStr << "#{speciesIdAry[sid]}/"
                    end
                    speciesIdStr << "#{speciesIdAry[-1]}"
                    exportBranch[b].updateSpeciesId(speciesIdStr)
                    # Find the one or two daughters of exportBranch[bb], and update their parentIds to the id of exportBranch[b], as exportBranch[bb] is about to be deleted.
                    # If exportBranch[bb] has no daughters, no parentIds need to be updated.
                    if exportBranch[bb].daughterId.size == 1
                      # Do a fast granddaughter search.
                      lower2 = 0
                      upper2 = exportBranch.size
                      granddaughterFound = false
                      until granddaughterFound
                        bbb = lower2+rand(upper2-lower2)
                        unless exportBranch[bbb] == nil
                          if exportBranch[bbb].id[1..-1].to_i < exportBranch[bb].daughterId[0][1..-1].to_i
                            lower2 = bbb
                          elsif exportBranch[bbb].id[1..-1].to_i > exportBranch[bb].daughterId[0][1..-1].to_i
                            upper2 = bbb
                          elsif exportBranch[bbb].id == exportBranch[bb].daughterId[0]
                            unless exportBranch[bbb].parentId == exportBranch[bb].id
                              raise "The daughter id of branch #{exportBranch[bb].id} is #{exportBranch[bb].daughterId[0]}, but the parent id of branch #{exportBranch[bbb].id} is #{exportBranch[bbb].parentId}!" # Should not be called at all.
                            end
                            granddaughterFound = true
                            exportBranch[bbb].updateParentId(exportBranch[b].id)
                            next
                          else
                            raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                          end
                        end
                      end
                    elsif exportBranch[bb].daughterId.size == 2 and exportBranch[bb].daughterId != ["none","none"]
                      # Do a fast search for the first granddaughter.
                      lower2 = 0
                      upper2 = exportBranch.size
                      granddaughter0Found = false
                      until granddaughter0Found
                        bbb = lower2+rand(upper2-lower2)
                        unless exportBranch[bbb] == nil
                          if exportBranch[bbb].id[1..-1].to_i < exportBranch[bb].daughterId[0][1..-1].to_i
                            lower2 = bbb
                          elsif exportBranch[bbb].id[1..-1].to_i > exportBranch[bb].daughterId[0][1..-1].to_i
                            upper2 = bbb
                          elsif exportBranch[bbb].id == exportBranch[bb].daughterId[0]
                            unless exportBranch[bbb].parentId == exportBranch[bb].id
                              raise "The daughter id of branch #{exportBranch[bb].id} is #{exportBranch[bb].daughterId[0]}, but the parent id of branch #{exportBranch[bbb].id} is #{exportBranch[bbb].parentId}!" # Should not be called at all.
                            end
                            granddaughter0Found = true
                            exportBranch[bbb].updateParentId(exportBranch[b].id)
                            next
                          else
                            raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                          end
                        end
                      end
                      # Do a fast search for the second granddaughter.
                      lower2 = 0
                      upper2 = exportBranch.size
                      granddaughter1Found = false
                      until granddaughter1Found
                        bbb = lower2+rand(upper2-lower2)
                        unless exportBranch[bbb] == nil
                          if exportBranch[bbb].id[1..-1].to_i < exportBranch[bb].daughterId[1][1..-1].to_i
                            lower2 = bbb
                          elsif exportBranch[bbb].id[1..-1].to_i > exportBranch[bb].daughterId[1][1..-1].to_i
                            upper2 = bbb
                          elsif exportBranch[bbb].id == exportBranch[bb].daughterId[1]
                            unless exportBranch[bbb].parentId == exportBranch[bb].id
                              raise "The daughter id of branch #{exportBranch[bb].id} is #{exportBranch[bb].daughterId[1]}, but the parent id of branch #{exportBranch[bbb].id} is #{exportBranch[bbb].parentId}!" # Should not be called at all.
                            end
                            granddaughter1Found = true
                            exportBranch[bbb].updateParentId(exportBranch[b].id)
                            next
                          else
                            raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                          end
                        end                          
                      end # until granddaughter1Found
                    end # if exportBranch[bb].daughterId.size == 1, elsif exportBranch[bb].daughterId.size == 2
                    exportBranch[bb] = nil
                    next
                  else
                    raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
                  end
                end
              end
              unmerged = true
            end
          end
        end
        exportBranch.compact!
      end
    end
    
    # Sort all branches according to origin.
    sortedBranch = exportBranch.sort { |a,b| b.origin <=> a.origin }
    sortedBranch.shift if removeStem

    # Determine the condition at the time of origin.
    if labels == true
      if @rootSplit == false and @treeReconstructed == false
        startLineages = 1
        if @branchRatesAssigned == false
          if branchLengths.downcase == "duration"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId}]:#{sortedBranch[0].duration.round(decimals)})"
          elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId}]:#{sortedBranch[0].expectedNumberOfSubstitutionsPerSite.round(decimals)})"
          elsif branchLengths.downcase == "actualnumberofsubstitutions"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId}]:#{sortedBranch[0].actualNumberOfSubstitutions.round(decimals)})"
          elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId}]:#{sortedBranch[0].actualNumberOfSubstitutionsPerSite.round(decimals)})"
          end
        else
          if branchLengths.downcase == "duration"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId},rate=#{sortedBranch[0].rate}]:#{sortedBranch[0].duration.round(decimals)})"
          elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId},rate=#{sortedBranch[0].rate}]:#{sortedBranch[0].expectedNumberOfSubstitutionsPerSite.round(decimals)})"
          elsif branchLengths.downcase == "actualnumberofsubstitutions"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId},rate=#{sortedBranch[0].rate}]:#{sortedBranch[0].actualNumberOfSubstitutions.round(decimals)})"
          elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId},rate=#{sortedBranch[0].rate}]:#{sortedBranch[0].actualNumberOfSubstitutionsPerSite.round(decimals)})"
          end
        end
      else
        startLineages = 2
        if @branchRatesAssigned == false
          if branchLengths.downcase == "duration"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId}]:#{sortedBranch[0].duration.round(decimals)},#{sortedBranch[1].id}[&label=#{sortedBranch[1].id},speciesId=#{sortedBranch[1].speciesId}]:#{sortedBranch[1].duration.round(decimals)})"
          elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId}]:#{sortedBranch[0].expectedNumberOfSubstitutionsPerSite.round(decimals)},#{sortedBranch[1].id}[&label=#{sortedBranch[1].id},speciesId=#{sortedBranch[1].speciesId}]:#{sortedBranch[1].expectedNumberOfSubstitutionsPerSite.round(decimals)})"
          elsif branchLengths.downcase == "actualnumberofsubstitutions"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId}]:#{sortedBranch[0].actualNumberOfSubstitutions.round(decimals)},#{sortedBranch[1].id}[&label=#{sortedBranch[1].id},speciesId=#{sortedBranch[1].speciesId}]:#{sortedBranch[1].actualNumberOfSubstitutions.round(decimals)})"
          elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId}]:#{sortedBranch[0].actualNumberOfSubstitutionsPerSite.round(decimals)},#{sortedBranch[1].id}[&label=#{sortedBranch[1].id},speciesId=#{sortedBranch[1].speciesId}]:#{sortedBranch[1].actualNumberOfSubstitutionsPerSite.round(decimals)})"
          end
        else
          if branchLengths.downcase == "duration"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId},rate=#{sortedBranch[0].rate}]:#{sortedBranch[0].duration.round(decimals)},#{sortedBranch[1].id}[&label=#{sortedBranch[1].id},speciesId=#{sortedBranch[1].speciesId},rate=#{sortedBranch[1].rate}]:#{sortedBranch[1].duration.round(decimals)})"
          elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId},rate=#{sortedBranch[0].rate}]:#{sortedBranch[0].expectedNumberOfSubstitutionsPerSite.round(decimals)},#{sortedBranch[1].id}[&label=#{sortedBranch[1].id},speciesId=#{sortedBranch[1].speciesId},rate=#{sortedBranch[1].rate}]:#{sortedBranch[1].expectedNumberOfSubstitutionsPerSite.round(decimals)})"
          elsif branchLengths.downcase == "actualnumberofsubstitutions"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId},rate=#{sortedBranch[0].rate}]:#{sortedBranch[0].actualNumberOfSubstitutions.round(decimals)},#{sortedBranch[1].id}[&label=#{sortedBranch[1].id},speciesId=#{sortedBranch[1].speciesId},rate=#{sortedBranch[1].rate}]:#{sortedBranch[1].actualNumberOfSubstitutions.round(decimals)})"
          elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
            newickString = "(#{sortedBranch[0].id}[&label=#{sortedBranch[0].id},speciesId=#{sortedBranch[0].speciesId},rate=#{sortedBranch[0].rate}]:#{sortedBranch[0].actualNumberOfSubstitutionsPerSite.round(decimals)},#{sortedBranch[1].id}[&label=#{sortedBranch[1].id},speciesId=#{sortedBranch[1].speciesId},rate=#{sortedBranch[1].rate}]:#{sortedBranch[1].actualNumberOfSubstitutionsPerSite.round(decimals)})"
          end
        end
        
      end

    else # if labels == false

      if @rootSplit == false and @treeReconstructed == false
        startLineages = 1
        if branchLengths.downcase == "duration"
          newickString = "(#{sortedBranch[0].id}:#{sortedBranch[0].duration.round(decimals)})"
        elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
          newickString = "(#{sortedBranch[0].id}:#{sortedBranch[0].expectedNumberOfSubstitutionsPerSite.round(decimals)})"
        elsif branchLengths.downcase == "actualnumberofsubstitutions"
          newickString = "(#{sortedBranch[0].id}:#{sortedBranch[0].actualNumberOfSubstitutions.round(decimals)})"
        elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
          newickString = "(#{sortedBranch[0].id}:#{sortedBranch[0].actualNumberOfSubstitutionsPerSite.round(decimals)})"
        end
      else
        startLineages = 2
        if branchLengths.downcase == "duration"
          newickString = "(#{sortedBranch[0].id}:#{sortedBranch[0].duration.round(decimals)},#{sortedBranch[1].id}:#{sortedBranch[1].duration.round(decimals)})"
        elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
          newickString = "(#{sortedBranch[0].id}:#{sortedBranch[0].expectedNumberOfSubstitutionsPerSite.round(decimals)},#{sortedBranch[1].id}:#{sortedBranch[1].expectedNumberOfSubstitutionsPerSite.round(decimals)})"
        elsif branchLengths.downcase == "actualnumberofsubstitutions"
          newickString = "(#{sortedBranch[0].id}:#{sortedBranch[0].actualNumberOfSubstitutions.round(decimals)},#{sortedBranch[1].id}:#{sortedBranch[1].actualNumberOfSubstitutions.round(decimals)})"
        elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
          newickString = "(#{sortedBranch[0].id}:#{sortedBranch[0].actualNumberOfSubstitutionsPerSite.round(decimals)},#{sortedBranch[1].id}:#{sortedBranch[1].actualNumberOfSubstitutionsPerSite.round(decimals)})"
        end
      end

    end

    newickString = newickString.dup
    # evolve the newick string. For each new pair of branches (i.e. at each speciation event) do the following
    startLineages.upto(sortedBranch.size-1) do |s|
      if (s-startLineages).even? == true

        #find sisters
        sister1 = sortedBranch[s]
        sister2 = sortedBranch[s+1]

        # Make sure the sisters are from the same parent.
        if sister1.parentId != sister2.parentId # this should not be called at all. sister1.parentId should always equal sister2.parentId.
          raise "The parent Id (#{sister1.parentId}) of sister 1 (#{sister1.id}) is different from the parent Id (#{sister2.parentId}) of sister 2 (#{sister2.id})!"
        end

        if labels == true
          findString = "#{sister1.parentId}["
          if @branchRatesAssigned == false
            if branchLengths.downcase == "duration"
              replaceString = "(#{sister1.id}[&label=#{sister1.id},speciesId=#{sister1.speciesId}]:#{sister1.duration.round(decimals)},#{sister2.id}[&label=#{sister2.id},speciesId=#{sister2.speciesId}]:#{sister2.duration.round(decimals)})["
            elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
              replaceString = "(#{sister1.id}[&label=#{sister1.id},speciesId=#{sister1.speciesId}]:#{sister1.expectedNumberOfSubstitutionsPerSite.round(decimals)},#{sister2.id}[&label=#{sister2.id},speciesId=#{sister2.speciesId}]:#{sister2.expectedNumberOfSubstitutionsPerSite.round(decimals)})["
            elsif branchLengths.downcase == "actualnumberofsubstitutions"
              replaceString = "(#{sister1.id}[&label=#{sister1.id},speciesId=#{sister1.speciesId}]:#{sister1.actualNumberOfSubstitutions.round(decimals)},#{sister2.id}[&label=#{sister2.id},speciesId=#{sister2.speciesId}]:#{sister2.actualNumberOfSubstitutions.round(decimals)})["
            elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
              replaceString = "(#{sister1.id}[&label=#{sister1.id},speciesId=#{sister1.speciesId}]:#{sister1.actualNumberOfSubstitutionsPerSite.round(decimals)},#{sister2.id}[&label=#{sister2.id},speciesId=#{sister2.speciesId}]:#{sister2.actualNumberOfSubstitutionsPerSite.round(decimals)})["
            end
          else
            if branchLengths.downcase == "duration"
              replaceString = "(#{sister1.id}[&label=#{sister1.id},speciesId=#{sister1.speciesId},rate=#{sister1.rate}]:#{sister1.duration.round(decimals)},#{sister2.id}[&label=#{sister2.id},speciesId=#{sister2.speciesId},rate=#{sister2.rate}]:#{sister2.duration.round(decimals)})["
            elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
              replaceString = "(#{sister1.id}[&label=#{sister1.id},speciesId=#{sister1.speciesId},rate=#{sister1.rate}]:#{sister1.expectedNumberOfSubstitutionsPerSite.round(decimals)},#{sister2.id}[&label=#{sister2.id},speciesId=#{sister2.speciesId},rate=#{sister2.rate}]:#{sister2.expectedNumberOfSubstitutionsPerSite.round(decimals)})["
            elsif branchLengths.downcase == "actualnumberofsubstitutions"
              replaceString = "(#{sister1.id}[&label=#{sister1.id},speciesId=#{sister1.speciesId},rate=#{sister1.rate}]:#{sister1.actualNumberOfSubstitutions.round(decimals)},#{sister2.id}[&label=#{sister2.id},speciesId=#{sister2.speciesId},rate=#{sister2.rate}]:#{sister2.actualNumberOfSubstitutions.round(decimals)})["
            elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
              replaceString = "(#{sister1.id}[&label=#{sister1.id},speciesId=#{sister1.speciesId},rate=#{sister1.rate}]:#{sister1.actualNumberOfSubstitutionsPerSite.round(decimals)},#{sister2.id}[&label=#{sister2.id},speciesId=#{sister2.speciesId},rate=#{sister2.rate}]:#{sister2.actualNumberOfSubstitutionsPerSite.round(decimals)})["
            end
          end
        else # if labels == false
          findString = "#{sister1.parentId}:"
          if branchLengths.downcase == "duration"
            replaceString = "(#{sister1.id}:#{sister1.duration.round(decimals)},#{sister2.id}:#{sister2.duration.round(decimals)}):"
          elsif branchLengths.downcase == "expectednumberofsubstitutionspersite"
            replaceString = "(#{sister1.id}:#{sister1.expectedNumberOfSubstitutionsPerSite.round(decimals)},#{sister2.id}:#{sister2.expectedNumberOfSubstitutionsPerSite.round(decimals)}):"
          elsif branchLengths.downcase == "actualnumberofsubstitutions"
            replaceString = "(#{sister1.id}:#{sister1.actualNumberOfSubstitutions.round(decimals)},#{sister2.id}:#{sister2.actualNumberOfSubstitutions.round(decimals)}):"
          elsif branchLengths.downcase == "actualnumberofsubstitutionspersite"
            replaceString = "(#{sister1.id}:#{sister1.actualNumberOfSubstitutionsPerSite.round(decimals)},#{sister2.id}:#{sister2.actualNumberOfSubstitutionsPerSite.round(decimals)}):"
          end
        end
        newickString.sub!(findString,replaceString)
      end
    end

    print "\rAdding species ids to tree string...       " if verbose
    exportBranch.size.times do |b|
      if exportBranch[b].extant or exportBranch[b].endCause == "extinction"
        if @treeReconstructed
          lastSpeciesId = exportBranch[b].speciesId.split("/")[-1]
          newickString.sub!("#{exportBranch[b].id}[","#{lastSpeciesId}[") if labels == true
          newickString.sub!("#{exportBranch[b].id}:","#{lastSpeciesId}:") if labels == false
        else
          newickString.sub!("#{exportBranch[b].id}[","#{exportBranch[b].speciesId}[") if labels == true
          newickString.sub!("#{exportBranch[b].id}:","#{exportBranch[b].speciesId}:") if labels == false
        end
      end
    end

    if fileName == nil
      if verbose
        print "\r"
      end
      newickString
    else
      print "\rWriting tree string in Newick format to '#{fileName}'..." if verbose
      file  = File.new(fileName,"w")
      output = ""
      if plain == true
        output << newickString + ";\n"
      else
        output << "#NEXUS\nbegin trees;\ntree tree1 = [&R] " + newickString + ";\n\nend;\n"
      end
      file.write(output)
      file.close
      if verbose
        endTime = Time.now
        puts "\rTree written in Newick format to '#{fileName}'.                 "
        if endTime-startTime < 60
          puts "Time used: #{(endTime-startTime).round(2)} seconds."
        elsif endTime-startTime < 3600
          puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
        else
          puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
        end
        puts
      end
    end
  end
  
  def visualize(fileName = "bd.tre", branchLengths = "duration", labels = true, plain = false, includeEmpty = true, overwrite = false, application = "/Applications/FigTree\\ v1.4.2.app", verbose = true)
    self.to_newick(fileName, branchLengths, labels, plain, includeEmpty, overwrite, verbose)
    system("open #{fileName} -a #{application}")
  end
  
  def dump(fileName = "bd.dmp", overwrite = false, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "----------------------------------------------------phylsim.rb | Export-----------------------------------------------------"
      puts
    end

    unless overwrite
      # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
      if File.exists?(fileName)
        puts "Ok to replace '#{fileName}'? (y/N)"
        answer = gets
        unless answer[0].downcase == "y"
          puts "Please specify a new filename:"
          fileName = gets.strip
        end
        puts
      end
    end

    File.new(fileName,"w").write(Marshal.dump(self))
    if verbose
      endTime = Time.now
      puts "Tree dumped to '#{fileName}'."
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts 
    end
  end

  def reconstruct(samplingScheme = "all", number = nil, focusGroup = nil, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "----------------------------------------------phylsim.rb | Tree reconstruction----------------------------------------------"
      puts
      print "Preparing tree reconstruction..."
    end

    # Handle arguments.
    @samplingScheme = samplingScheme
    @focusGroup = focusGroup

    # Handle the second argument so that only n is used.
    n = 0
    ro = 0
    unless number == nil
      if number.to_f < 1
        ro = number.to_f
      elsif number.to_f >= 1
        n = number.round
      else
        raise "The number parameter could not be read. It should either be a floating point decimal between 0 and 1 to indicate sampling fraction, or an integer >= 1 to indicate the total number of sampled extant species!"
      end
    end
    if ro > 0
      @Np.times do
        n += 1 if rand < ro
      end
    end

    # Deal with the focusGroup argument
    unless @focusGroup == nil
      focusGroupFound = false
      minExtantFGSpecies = 0
      maxExtantFGSpecies = @Np
      sampledFGSpecies = 0
      minFGOrigin = 0
      maxFGOrigin = @posteriorTreeOrigin
      # The focusGroup argument should look like this: focusGroup = "minExtantFGSpecies:2000,maxExtantFGSpecies:3000,sampledFGSpecies:579,minFGOrigin:50,maxFGOrigin:150"
      if @focusGroup.class == Array
        if @focusGroup.size == 3
          unless @focusGroup[0].class == Fixnum and @focusGroup[1].class == Fixnum and @focusGroup[2].class == Fixnum
            raise "When specifying the minimum number of extant species in the focus group, the maximum number of extant species in that group, and the number of sampled species from that group as an array, all three numbers must be integer numbers."
          end
          unless @focusGroup[0] >= 0 and @focusGroup[1] >= 0 and @focusGroup[2] >= 0
            raise "When specifying the minimum number of extant species in the focus group, the maximum number of extant species in that group, and the number of sampled species from that group as an array, all three numbers must be positive numbers."
          end
          if @focusGroup[2] > @focusGroup[0]
            raise "The number of sampled species from the focus group must be larger or equal to the minimum number of extant species of that group!"
          elsif @focusGroup[0] > @focusGroup[1]
            raise "The minimum number of extant species in the focus group must be smaller or equal to the maximum number of extant species of that group!"
          else
            minExtantFGSpecies = @focusGroup[0]
            maxExtantFGSpecies = @focusGroup[1]
            sampledFGSpecies = @focusGroup[2]
          end
        elsif @focusGroup.size == 4
          unless @focusGroup[0].class == Fixnum and @focusGroup[1].class == Fixnum and @focusGroup[2].class == Fixnum and @focusGroup[3].class == Fixnum
            raise "When specifying the minimum number of extant species in the focus group, the maximum number of extant species in that group, and the number of sampled species, and the minimum time of origin of that group as an array, all four numbers must be integer numbers."
          end
          unless @focusGroup[0] >= 0 and @focusGroup[1] >= 0 and @focusGroup[2] >= 0 and @focusGroup[3] >= 0
            raise "When specifying the minimum number of extant species in the focus group, the maximum number of extant species in that group, and the number of sampled species, and the minimum time of origin of that group as an array, all four numbers must be positive numbers."
          end
          if @focusGroup[2] > @focusGroup[0]
            raise "The number of sampled species from the focus group must be larger or equal to the minimum number of extant species of that group!"
          elsif @focusGroup[0] > @focusGroup[1]
            raise "The minimum number of extant species in the focus group must be smaller or equal to the maximum number of extant species of that group!"
          elsif @focusGroup[3] >= @posteriorTreeOrigin
            raise "The minimum time of origin of the focus group must be smaller than the origin of the full tree!"
          else
            minExtantFGSpecies = @focusGroup[0]
            maxExtantFGSpecies = @focusGroup[1]
            sampledFGSpecies = @focusGroup[2]
            minFGOrigin = @focusGroup[3]
          end
        elsif @focusGroup.size == 5
          unless @focusGroup[0].class == Fixnum and @focusGroup[1].class == Fixnum and @focusGroup[2].class == Fixnum and @focusGroup[3].class == Fixnum and @focusGroup[4].class == Fixnum
            raise "When specifying the minimum number of extant species in the focus group, the maximum number of extant species in that group, and the number of sampled species, and the minimum and maximum times of origin of that group as an array, all five numbers must be integer numbers."
          end
          unless @focusGroup[0] >= 0 and @focusGroup[1] >= 0 and @focusGroup[2] >= 0 and @focusGroup[3] >= 0 and @focusGroup[4] >= 0
            raise "When specifying the minimum number of extant species in the focus group, the maximum number of extant species in that group, and the number of sampled species, and the minimum and maximum times of origin of that group as an array, all five numbers must be positive numbers."
          end
          if @focusGroup[2] > @focusGroup[0]
            raise "The number of sampled species from the focus group must be larger or equal to the minimum number of extant species of that group!"
          elsif @focusGroup[0] > @focusGroup[1]
            raise "The minimum number of extant species in the focus group must be smaller or equal to the maximum number of extant species of that group!"
          elsif @focusGroup[3] >= @posteriorTreeOrigin
            raise "The minimum time of origin of the focus group must be smaller than the origin of the full tree!"
          elsif @focusGroup[3] >= @focusGroup[4]
            raise "The minimum time of origin of the focus group must be smaller than the maximum time of origin of that group!"
          else
            minExtantFGSpecies = @focusGroup[0]
            maxExtantFGSpecies = @focusGroup[1]
            sampledFGSpecies = @focusGroup[2]
            minFGOrigin = @focusGroup[3]
            maxFGOrigin = @focusGroup[4]
          end
        else
          raise "Please specify the focus group either as an array of three or five number that indicate minimum number of extant species, maximum number of extant species, number of sampled species, and optionally minimum and maximum times of origin of the focus group, or as a string according to the following example: \"minExtantFGSpecies:2000,maxExtantFGSpecies:3000,sampledFGSpecies:579,minFGOrigin:50,maxFGOrigin:150\"."
        end
      elsif @focusGroup.class == String
        focusGroupAry = @focusGroup.split(",")
        focusGroupAry.each do |fga|
          tmp = fga.split(":")
          raise "The string specified for the focus group parameters could not be understood!" unless tmp.size == 2
          if tmp[0].downcase == "minExtantFGSpecies"
            minExtantFGSpecies = tmp[1].strip.to_i
          elsif tmp[0].downcase == "maxExtantFGSpecies"
            maxExtantFGSpecies = tmp[1].strip.to_i
          elsif tmp[0].downcase == "sampledFGSpecies"
            sampledFGSpecies = tmp[1].strip.to_i
          elsif tmp[0].downcase = "minFGOrigin"
            minFGOrigin = tmp[1].strip.to_i
          elsif tmp[0].downcase = "maxFGOrigin"
            maxFGOrigin = tmp[1].strip.to_i
          end
        end
        if sampledFGSpecies > minExtantFGSpecies
          raise "The number of sampled species from the focus group must be smaller or equal to the minimum number of extant species of that group!"
        end
        if minExtantFGSpecies > maxExtantFGSpecies
          raise "The minimum number of extant species in the focus group must be smaller or equal to the maximum number of extant species of that group!"
        end
        if minFGOrigin >= @posteriorTreeOrigin
          raise "The minimum time of origin of the focus group must be smaller than the origin of the full tree!"
        end
        if minFGOrigin >= maxFGOrigin
          raise "The minimum time of origin of the focus group must be smaller than the maximum time of origin of that group!"
        end
      else
        raise "Please specify the focus group either as an array of three or five number that indicate minimum number of extant species, maximum number of extant species, number of sampled species, and optionally minimum and maximum times of origin of the focus group, or as a string according to the following example: \"minExtantFGSpecies:2000,maxExtantFGSpecies:3000,sampledFGSpecies:579,minFGOrigin:50,maxFGOrigin:15\"."
      end
    end
    
    # Initiate an array of all extant branches and prepare an array of sampled extant branches that is to be filled.
    extantBranch = []
    @branch.each do |b|
      extantBranch << b if b.extant
    end
    raise "extantBranch.size is not equal to @Np" unless extantBranch.size == @Np # should not be called at all.
    sampledExtantBranch = []

    # Deal with some odd arguments.
    if @samplingScheme == "all" and @focusGroup != nil
      warn "WARNING: Because the sampling scheme has been chosen to include all extant species, the parameters specified for the focus group will be ignored."
    end
    n = @Np if @samplingScheme == "all"
    if @focusGroup != nil and n < sampledFGSpecies
      raise "The specified total number of sampled species (#{n}) is smaller than the specified number of species to be sampled from the focus group (sampledFGSpecies)!"
    end
    if n >= @Np
      if verbose
        puts "\rModel:                                           "
        puts "Full sampling scheme"
        puts
        puts "Parameters:"
        puts "np_sampled/np = #{n}/#{@Np}"
        puts
      end
      sampledExtantBranch = extantBranch
    elsif n < 0
      raise "No species are sampled (n<0)!"
    elsif n == 0
      raise "No species are sampled (n=0)!"
    elsif n == 1
      if verbose
        puts "\rModel:                                           "
        puts "Birth-death-sampling_m process (Stadler 2009)"
        puts
        puts "Parameters:"
        puts "n/m = #{n}/#{@Np}"
        puts
      end
      sampledExtantBranch << extantBranch[rand(@Np)]
    else
      # a random sampling scheme (see Stadler 2009)
      if @samplingScheme == "random"
        if verbose
          puts "\rModel:                                           "
          if number.to_f < 1
            puts "Birth-death-sampling_ro process (Stadler 2009)"
          else
            puts "Birth-death-sampling_m process (Stadler 2009)"
          end
          puts
          puts "Parameters:"
          puts "ro = #{ro}" if number.to_f < 1
          puts "n/m = #{n}/#{@Np}"
          puts
          if @focusGroup != nil
            puts "Conditions on focus group:"
            puts "#{minExtantFGSpecies-1} < number of extant species in focus group < #{maxExtantFGSpecies+1}"
            puts "#{minFGOrigin} < origin of focus group < #{maxFGOrigin}"
            puts "Number of species to be sampled from focus group: #{sampledFGSpecies}"
          end
        end

        if @focusGroup != nil
          # Search for potential focus groups.
          print "Searching for potential focus groups..." if verbose
          candidateFGOriginBranches = []
          @branch.each do |b|
            if minExtantFGSpecies <= b.extantProgenyId.size and maxExtantFGSpecies >= b.extantProgenyId.size and minFGOrigin <= b.termination and maxFGOrigin >= b.termination
              candidateFGOriginBranches << b
            end
          end
          if candidateFGOriginBranches.size == 0
            if verbose
              puts "\rPotential focus groups in tree: 0 (stopping run...)"
              puts
            end
            raise RetryException.new(true), "The tree contains no clades that match the focus group criteria! Try again."
          elsif candidateFGOriginBranches.size == 1
            puts "\rPotential focus groups in tree: 1          " if verbose
          else
            puts "\rPotential focus groups in tree: #{candidateFGOriginBranches.size} (one of them is chosen at random)" if verbose
          end
          selectedFGOriginBranch = candidateFGOriginBranches.sample
          if verbose
            puts
            print "\rAdding the focus group parameter to all branches..."
          end
          @branch.each do |b|
            if selectedFGOriginBranch.progenyId.include?(b.id)
              b.addFocusGroup(focusGroup = true)
            else
              b.addFocusGroup(focusGroup = false)
            end
          end
        end

        if n == 2
          if @rootSplit
            print "Randomly sampling one species from each side of the root..." if verbose
            b0ExtantProgenyId = []
            b1ExtantProgenyId = []
            b0Found = false
            b = 0
            until b0Found
              if @branch[b].id == "b0"
                b0ExtantProgenyId = @branch[b].extantProgenyId
                b0Found = true
              end
              b += 1
            end
            b1Found = false
            b = 0
            until b1Found
              if @branch[b].id == "b1"
                b1ExtantProgenyId = @branch[b].extantProgenyId
                b1Found = true
              end
              b += 1
            end
            candidate0 = []
            candidate1 = []
            extantBranch.each do |eb|
              candidate0 << eb if b0ExtantProgenyId.include?(eb.id)
              candidate1 << eb if b1ExtantProgenyId.include?(eb.id)
            end
            sampledExtantBranch << candidate0[rand(candidate0.size)]
            sampledExtantBranch << candidate1[rand(candidate1.size)]
          else
            print "Randomly sampling two species..." if verbose
            n.times do |j|
              position = rand(extantBranch.size)
              sampledExtantBranch << extantBranch.delete_at(position)
            end
          end
          raise "Please specify a larger number of species to be sampled (n > 2) when including a focus group!" if @focusGroup != nil
        else
          if @rootSplit
            print "\rPreparing sampling...                                                   " if verbose
            b0ExtantProgenyId = []
            b1ExtantProgenyId = []
            b0Found = false
            b = 0
            until b0Found
              if @branch[b].id == "b0"
                b0ExtantProgenyId = @branch[b].extantProgenyId
                b0Found = true
              end
              b += 1
            end
            b1Found = false
            b = 0
            until b1Found
              if @branch[b].id == "b1"
                b1ExtantProgenyId = @branch[b].extantProgenyId
                b1Found = true
              end
              b += 1
            end
            if @focusGroup != nil
              # Divide array extantBranch into extantFGBranch and extantNonFGBranch, then sample sampledFGSpecies from extantFGBranch and (n-sampledFGSpecies) from extantNonFGBranch.
              extantFGBranch = []
              extantNonFGBranch = []
              extantBranch.each do |eb|
                if eb.focusGroup
                  extantFGBranch << eb
                else
                  extantNonFGBranch << eb
                end
              end
              print "\rRandomly sampling #{sampledFGSpecies} species from the focus group..." if verbose
              sampledFGSpecies.times do
                position = rand(extantFGBranch.size)
                sampledExtantBranch << extantFGBranch.delete_at(position)
              end
              # Sort all extant branches that are not in the focus group into candidate0 and candidate1.
              print "\rSampling #{n-sampledFGSpecies} more species from outside the focus group at random, but proportionally from both sides of the root..."
              candidate0 = []
              candidate1 = []
              extantNonFGBranch.each do |eb|
                candidate0 << eb if b0ExtantProgenyId.include?(eb.id)
                candidate1 << eb if b1ExtantProgenyId.include?(eb.id)
              end
              n0 = ((n-sampledFGSpecies)*(candidate0.size/(extantNonFGBranch.size).to_f)).round
              n1 = ((n-sampledFGSpecies)*(candidate1.size/(extantNonFGBranch.size).to_f)).round
              if n0 + n1 == (n-sampledFGSpecies)+1
                if rand > 0.5
                  n0 -= 1
                else
                  n1 -= 1
                end
              end
              raise "n0 + n1 is not equal to n - sampledFGSpecies" if n0 + n1 != n - sampledFGSpecies
            elsif @focusGroup == nil
              print "\rSampling #{n} species at random, but proportionally from both sides of the root..." if verbose
              candidate0 = []
              candidate1 = []
              extantBranch.each do |eb|
                candidate0 << eb if b0ExtantProgenyId.include?(eb.id)
                candidate1 << eb if b1ExtantProgenyId.include?(eb.id)
              end
              n0 = (n*(candidate0.size/@Np.to_f)).round
              n1 = (n*(candidate1.size/@Np.to_f)).round
              if n0 + n1 == n+1
                if rand > 0.5
                  n0 -= 1
                else
                  n1 -= 1
                end
              end
              raise "n0 + n1 is not equal to n" if n0 + n1 != n
            end
            n0.times do
              position = rand(candidate0.size)
              sampledExtantBranch << candidate0.delete_at(position)
            end
            n1.times do
              position = rand(candidate1.size)
              sampledExtantBranch << candidate1.delete_at(position)
            end
          else
            if @focusGroup != nil
              # Divide array extantBranch into extantFGBranch and extantNonFGBranch, then sample sampledFGSpecies from extantFGBranch and (n-sampledFGSpecies) from extantNonFGBranch.
              print "\rRandomly sampling #{sampledFGSpecies} species from the focus group..." if verbose
              extantFGBranch = []
              extantNonFGBranch = []
              extantBranch.each do |eb|
                if eb.focusGroup
                  extantFGBranch << eb
                else
                  extantNonFGBranch << eb
                end
              end
              sampledFGSpecies.times do
                position = rand(extantFGBranch.size)
                sampledExtantBranch << extantFGBranch.delete_at(position)
              end
              print "\rRandomly sampling #{n-sampledFGSpecies} more species from outside the focus group..." if verbose
              (n-sampledFGSpecies).times do
                position = rand(extantNonFGBranch.size)
                sampledExtantBranch << extantNonFGBranch.delete_at(position)
              end
            elsif @focusGroup == nil
              print "Randomly sampling #{n-1} more species..." if verbose
              n.times do
                position = rand(extantBranch.size)
                sampledExtantBranch << extantBranch.delete_at(position)
              end
            end
          end
        end

      elsif @samplingScheme == "diversified" # Implements the diversified sampling scheme of Hoehna et al. (2011).

        if verbose
          puts "\rModel:                                       "
          puts "Diversified sampling (Hoehna et al. 2011)"
          puts
          puts "Parameters:"
          puts "n/m = #{n}/#{@Np}"
          puts
          if @focusGroup != nil
            puts "Conditions on focus group:"
            puts "#{minExtantFGSpecies-1} < number of extant species in focus group < #{maxExtantFGSpecies+1}"
            puts "#{minFGOrigin} < origin of focus group < #{maxFGOrigin}"
            puts "Number of species to be sampled from focus group: #{sampledFGSpecies}"
          end
        end

        if @focusGroup != nil
          # Search for potential focus groups
          if verbose
            print "Searching for potential focus groups..."
          end
          candidateFGOriginBranches = []
          @branch.each do |b|
            if minExtantFGSpecies <= b.extantProgenyId.size and maxExtantFGSpecies >= b.extantProgenyId.size and minFGOrigin <= b.termination and maxFGOrigin >= b.termination
              candidateFGOriginBranches << b
            end
          end
          if candidateFGOriginBranches.size == 0
            if verbose
              puts "\rPotential focus groups in tree: 0 (stopping run...)"
              puts
            end
            raise RetryException.new(true), "The tree contains no clades that match the focus group criteria! Try again."
          elsif candidateFGOriginBranches.size == 1
            puts "\rPotential focus groups in tree: 1          " if verbose
          else
            puts "\rPotential focus groups in tree: #{candidateFGOriginBranches.size} (one of them is chosen at random)" if verbose
          end
          selectedFGOriginBranch = candidateFGOriginBranches.sample
          if verbose
            puts
            print "\rAdding the focus group parameter to all branches..."
          end
          @branch.each do |b|
            if selectedFGOriginBranch.progenyId.include?(b.id)
              b.addFocusGroup(focusGroup = true)
            else
              b.addFocusGroup(focusGroup = false)
            end
          end
          
        end

        print "\rPreparing sampling...                                                   " if verbose
        # Let's make sure that the array @branch is still sorted by the id numbers. This is important for the quick parent search below.
        branchesSorted = true
        (@branch.size-1).times do |b|
          if @branch[b].id[1..-1].to_i > @branch[b+1].id[1..-1].to_i
            branchesSorted = false
            break
          end
        end
        if branchesSorted == false
          # This is unlikely to happen anyway, and should be checked if it does.
          warn "WARNING: For some reason, branches are not sorted according to their id. Better check out why."
          until branchesSorted
            allSortedSoFar = true
            (@branch.size-1).times do |b|
              if @branch[b].id[1..-1].to_i > @branch[b+1].id[1..-1].to_i
                @branch[b],@branch[b+1] = @branch[b+1],@branch[b]
                allSortedSoFar = false
              end
            end
            branchesSorted = true if allSortedSoFar == true
          end
        end

        # The tree must first be reconstructed with all extant taxa.
        # Branches with endCause 'extinction' are deleted, and if they're the only daughter of the parent, then the parent's endCause also becomes 'extinction'.
        # This is done until no branches have endCause 'extinction' anymore.
        extinction = true
        while extinction == true
          extinction = false
          @branch.size.times do |b|
            if @branch[b].endCause == "extinction"
              # Fast parent search to find the branch that has the id that matches @branch[b].parentId. This is complicated by the fact that a branch can now be nil.
              lower = 0
              upper = @branch.size
              parentFound = false
              until parentFound == true
                # Use the fact that @branch is sorted according to id for a quick parent search.
                bb = lower+rand(upper-lower)
                unless @branch[bb] == nil
                if @branch[bb].id[1..-1].to_i < @branch[b].parentId[1..-1].to_i
                    lower = bb
                  elsif @branch[bb].id[1..-1].to_i > @branch[b].parentId[1..-1].to_i
                    upper = bb
                  elsif @branch[bb].id == @branch[b].parentId
                    parentFound = true
                    if @branch[bb].daughterId.size == 2
                      daughterId = @branch[bb].daughterId.dup
                      daughterId.delete(@branch[b].id)
                      @branch[bb].updateDaughterId(daughterId)
                    elsif @branch[bb].daughterId.size == 1
                      @branch[bb].updateDaughterId([])
                      @branch[bb].updateEndCause("extinction")
                    else
                      raise "branch[#{bb}] has #{@branch[bb].daughterId.size} daughters!" # should never be called at all.
                    end
                    @branch[bb].updateFossil(@branch[bb].fossil+@branch[b].fossil) if @fossilRecordAdded
                    next
                  else
                    raise "Some problem occurred with the fast parent search!" # This should never be called at all.
                  end
                end
              end
              extinction = true
              @branch[b] = nil
            end
          end
          @branch.compact!
        end
      
        # Delete degree one vertices.
        unmerged = true
        while unmerged == true
          unmerged = false
          @branch.size.times do |b|
            unless @branch[b] == nil
              if @branch[b].daughterId.size == 1
                lower = 0
                upper = @branch.size
                daughterFound = false
                until daughterFound
                  # Use the fact that @branch is sorted according to id for a quick daughter search.
                  bb = lower+rand(upper-lower)
                  unless @branch[bb] == nil
                    if @branch[bb].id[1..-1].to_i < @branch[b].daughterId[0][1..-1].to_i
                      lower = bb
                    elsif @branch[bb].id[1..-1].to_i > @branch[b].daughterId[0][1..-1].to_i
                      upper = bb
                    elsif @branch[bb].id == @branch[b].daughterId[0]
                      daughterFound = true
                      @branch[b].updateTermination(@branch[bb].termination)
                      @branch[b].updateDaughterId(@branch[bb].daughterId)
                      @branch[b].updateEndCause(@branch[bb].endCause)
                      @branch[b].updateExtant(@branch[bb].extant)
                      @branch[b].updateRate((@branch[b].duration*@branch[b].rate+@branch[bb].duration*@branch[bb].rate)/(@branch[b].duration+@branch[bb].duration)) if @branchRatesAssigned
                      @branch[b].updateDuration(@branch[b].origin - @branch[b].termination)
                      speciesIdAry = @branch[b].speciesId.split("/")
                      speciesIdAry1 = @branch[bb].speciesId.split("/")
                      speciesIdAry1.each do |s|
                        speciesIdAry << s unless speciesIdAry.include?(s)
                      end
                      speciesIdStr = ""
                      (speciesIdAry.size-1).times do |sid|
                        speciesIdStr << "#{speciesIdAry[sid]}/"
                      end
                      speciesIdStr << "#{speciesIdAry[-1]}"
                      @branch[b].updateSpeciesId(speciesIdStr)
                      @branch[b].updateFossil(@branch[b].fossil + @branch[bb].fossil) if @fossilRecordAdded
                      # Find the one or two daughters of @branch[bb], and update their parentIds to the id of @branch[b], as @branch[bb] is about to be deleted.
                      # If @branch[bb] has no daughters, no parentIds need to be updated.
                      if @branch[bb].daughterId.size == 1
                        # Do a fast granddaughter search.
                        lower2 = 0
                        upper2 = @branch.size
                        granddaughterFound = false
                        until granddaughterFound
                          bbb = lower2+rand(upper2-lower2)
                          unless @branch[bbb] == nil
                            if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[0][1..-1].to_i
                              lower2 = bbb
                            elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[0][1..-1].to_i
                              upper2 = bbb
                            elsif @branch[bbb].id == @branch[bb].daughterId[0]
                              unless @branch[bbb].parentId == @branch[bb].id
                                raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[0]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                              end
                              granddaughterFound = true
                              @branch[bbb].updateParentId(@branch[b].id)
                              next
                            else
                              raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                            end
                          end
                        end
                      elsif @branch[bb].daughterId.size == 2 and @branch[bb].daughterId != ["none","none"]
                        # Do a fast search for the first granddaughter.
                        lower2 = 0
                        upper2 = @branch.size
                        granddaughter0Found = false
                        until granddaughter0Found
                          bbb = lower2+rand(upper2-lower2)
                          unless @branch[bbb] == nil
                            if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[0][1..-1].to_i
                              lower2 = bbb
                            elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[0][1..-1].to_i
                              upper2 = bbb
                            elsif @branch[bbb].id == @branch[bb].daughterId[0]
                              unless @branch[bbb].parentId == @branch[bb].id
                                raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[0]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                              end
                              granddaughter0Found = true
                              @branch[bbb].updateParentId(@branch[b].id)
                              next
                            else
                              raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                            end
                          end
                        end
                        # Do a fast search for the second granddaughter.
                        lower2 = 0
                        upper2 = @branch.size
                        granddaughter1Found = false
                        until granddaughter1Found
                          bbb = lower2+rand(upper2-lower2)
                          unless @branch[bbb] == nil
                            if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[1][1..-1].to_i
                              lower2 = bbb
                            elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[1][1..-1].to_i
                              upper2 = bbb
                            elsif @branch[bbb].id == @branch[bb].daughterId[1]
                              unless @branch[bbb].parentId == @branch[bb].id
                                raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[1]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                              end
                              granddaughter1Found = true
                              @branch[bbb].updateParentId(@branch[b].id)
                              next
                            else
                              raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                            end
                          end                          
                        end # until granddaughter1Found
                      end # if @branch[bb].daughterId.size == 1, elsif @branch[bb].daughterId.size == 2
                      @branch[bb] = nil
                      next
                    else
                      raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
                    end
                  end
                end
                unmerged = true
              end
            end
          end
          @branch.compact!
        end

        # Now we should have a reconstructed tree of all extant lineages.
        # Make sure all branches have exactly two daughters.
        @branch.each do |b|
          raise "Branch #{b.id} has #{b.daughterId.size} daughters!" unless b.daughterId.size == 2
        end

        # Deleted branches must be removed from progenyId and extantProgenyId.
        @branch.each do |b|
          newProgenyId = []
          newExtantProgenyId = []
          progenyId = b.progenyId
          progenyId.each do |p|
            # Do another fast id search to find the branch, whose id equals p.
            lower = -1
            upper = @branch.size
            progenyFound = false
            until progenyFound
              if upper-lower == 1
                progenyFound = true
              else
                # Use the fact that @branch is sorted according to id for a quick parent search.
                bb = lower+(upper-lower)/2
                if @branch[bb].id[1..-1].to_i < p[1..-1].to_i
                  lower = bb
                elsif @branch[bb].id[1..-1].to_i > p[1..-1].to_i
                  upper = bb
                elsif @branch[bb].id == p
                  progenyFound = true
                  newProgenyId << p
                  newExtantProgenyId << p if @branch[bb].extant
                  next
                else
                  raise "Some problem occurred with the fast progeny search!" # This should never be called at all.
                end
              end
            end
          end
          b.updateProgenyId(newProgenyId)
          b.updateExtantProgenyId(newExtantProgenyId)
        end

        # The array of all extant branches needs to be filled again, as some have been renamed in the previous reconstruction.
        extantBranch = []
        @branch.each do |b|
          extantBranch << b if b.extant
        end
        raise "extantBranch.size is not equal to @Np" unless extantBranch.size == @Np # Should not be called at all.

        if verbose
          if @focusGroup != nil
            print "\rRandomly sampling one species from each side of the root, one of them from the focus group..."
          else
            print "\rRandomly sampling one species from each side of the root..."
          end
        end

        # Start with sampling two extant species, each from one side of the root to make sure the root node is included.
        extProgId0 = []
        extProgId1 = []
        if @rootSplit
          b0Found = false
          b = 0
          until b0Found
            if @branch[b].id == "b0"
              if @branch[b].extant
                extProgId0 = [@branch[b].id]
              else
                extProgId0 = @branch[b].extantProgenyId
              end
              b0Found = true
            end
            b += 1
          end
          b1Found = false
          b = 0
          until b1Found
            if @branch[b].id == "b1"
              if @branch[b].extant
                extProgId1 = [@branch[b].id]
              else
                extProgId1 = @branch[b].extantProgenyId
              end
              b1Found = true
            end
            b += 1
          end
        else
          b0Id = nil
          b1Id = nil
          b0Found = false
          b = 0
          until b0Found
            if @branch[b].id == "b0"
              b0Id = @branch[b].daughterId[0]
              b1Id = @branch[b].daughterId[1]
              b0Found = true
            end
            b += 1
          end
          newB0Found = false
          b = 0
          until newB0Found
            if @branch[b].id == b0Id
              if @branch[b].extant
                extProgId0 = [@branch[b].id]
              else
                extProgId0 = @branch[b].extantProgenyId
              end
              newB0Found = true
            end
            b += 1
          end
          newB1Found = false
          b = 0
          until newB1Found
            if @branch[b].id == b1Id
              if @branch[b].extant
                extProgId1 = [@branch[b].id]
              else
                extProgId1 = @branch[b].extantProgenyId
              end
              newB1Found = true
            end
            b += 1
          end
        end

        # If we have focus group, make sure that one of the selected ids is drawn from it. This is because the focus group must necessarily be in one of
        # the partial trees. To do that, first find out in which part of the tree the focus group is.
        if @focusGroup != nil
          # We divide the array extantBranch into extantFGBranch and extantNonFGBranch, that will facilitate things later.
          extantFGBranch = []
          extantNonFGBranch = []
          extantBranch.each do |eb|
            if eb.focusGroup
              extantFGBranch << eb
            else
              extantNonFGBranch << eb
            end
          end
          # If the first branch of array extantFGBranch is included in extProgId0, then the entire focus group must necessarily be in it.
          if extProgId0.include?(extantFGBranch[0].id)
            selectedId0 = extantFGBranch.sample.id
            selectedId1 = extProgId1.sample
          elsif extProgId1.include?(extantFGBranch[0].id)
            selectedId0 = extProgId0.sample
            selectedId1 = extantFGBranch.sample.id
          else
            raise "One extant focus group branch is not included in both progenies of the two most basal branches!"
          end
        elsif @focusGroup == nil
          selectedId0 = extProgId0.sample
          selectedId1 = extProgId1.sample
        end

        # Do another fast id search to find the extant branch, whose id equals selectedId0 (extantBranch must also be sorted according to id).
        lower = -1
        upper = extantBranch.size
        idFound = false
        until idFound
          eb = lower+(upper-lower)/2
          if extantBranch[eb].id[1..-1].to_i < selectedId0[1..-1].to_i
            lower = eb
          elsif extantBranch[eb].id[1..-1].to_i > selectedId0[1..-1].to_i
            upper = eb
          elsif extantBranch[eb].id == selectedId0
            idFound = true
            sampledExtantBranch << extantBranch.delete_at(eb)
            next
          else
            raise "Some problem occurred with the fast id search!" # This should never be called at all.
          end
        end
        # Do another fast id search to find the extant branch, whose id equals selectedId1.
        lower = -1
        upper = extantBranch.size
        idFound = false
        until idFound
          eb = lower+(upper-lower)/2
          if extantBranch[eb].id[1..-1].to_i < selectedId1[1..-1].to_i
            lower = eb
          elsif extantBranch[eb].id[1..-1].to_i > selectedId1[1..-1].to_i
            upper = eb
          elsif extantBranch[eb].id == selectedId1
            idFound = true
            sampledExtantBranch << extantBranch.delete_at(eb)
            next
          else
            raise "Some problem occurred with the fast id search!" # This should never be called at all.
          end
        end

        # Sort all branches according to their termination date.
        sortedBranch = @branch.sort { |a,b| b.termination <=> a.termination }
        sortedBranch.shift unless @rootSplit

        if @focusGroup != nil
          print "Sampling #{sampledFGSpecies-1} more species from the focus group according to the diversified sampling scheme..." if verbose
          # Extract the array sortedFGBranch from array sortedBranch.
          sortedFGBranch = []
          sortedBranch.each do |sb|
            sortedFGBranch << sb if sb.focusGroup
          end
          # Also extract the array fGBranch from @branch (fGBranch must be sorted according to id if @branch is).
          fGBranch = []
          @branch.each do |b|
            fGBranch << b if b.focusGroup
          end

          # First, apply the diversified sampling to the focus group, so that sampledFGSpecies are sampled from it.
          (sampledFGSpecies-1).times do |fgn| # We need to sampledFGSpecies-1 nodes in order to end up with sampledFGSpecies species.
            # Determine the extant progeny of both sides of the n1th split.
            extProgId0 = []
            extProgId1 = []
            # Do another fast id search to find the branch among fGBranch, whose id equals sortedFGBranch[fgn].daughterId[0].
            lower = -1
            upper = fGBranch.size
            d0Found = false
            until d0Found
              b = lower+(upper-lower)/2
              if fGBranch[b].id[1..-1].to_i < sortedFGBranch[fgn].daughterId[0][1..-1].to_i
                lower = b
              elsif fGBranch[b].id[1..-1].to_i > sortedFGBranch[fgn].daughterId[0][1..-1].to_i
                upper = b
              elsif fGBranch[b].id == sortedFGBranch[fgn].daughterId[0]
                d0Found = true
                if fGBranch[b].extant
                  extProgId0 = [fGBranch[b].id]
                else
                  extProgId0 = fGBranch[b].extantProgenyId
                end
                next
              else
                raise "Some problem occurred with the fast id search!" # This should never be called at all.
              end
            end
            # Do another fast id search to find the branch among fGBranch, whose id equals sortedFGBranch[fgn].daughterId[1].
            lower = -1
            upper = fGBranch.size
            d1Found = false
            until d1Found
              b = lower+(upper-lower)/2
              if fGBranch[b].id[1..-1].to_i < sortedFGBranch[fgn].daughterId[1][1..-1].to_i
                lower = b
              elsif fGBranch[b].id[1..-1].to_i > sortedFGBranch[fgn].daughterId[1][1..-1].to_i
                upper = b
              elsif fGBranch[b].id == sortedFGBranch[fgn].daughterId[1]
                d1Found = true
                if fGBranch[b].extant
                  extProgId1 = [fGBranch[b].id]
                else
                  extProgId1 = fGBranch[b].extantProgenyId
                end
                next
              else
                raise "Some problem occurred with the fast id search!" # This should never be called at all.
              end
            end

            # See whether any of the previously sampled extant branches are in the extant progenies of one of the sides of the split.
            extProgId0AlreadyIncluded = false
            sampledExtantBranch.size.times do |seb|
              if extProgId0.include?(sampledExtantBranch[seb].id)
                extProgId0AlreadyIncluded = true
                break
              end
            end
            extProgId1AlreadyIncluded = false
            sampledExtantBranch.size.times do |seb|
              if extProgId1.include?(sampledExtantBranch[seb].id)
                extProgId1AlreadyIncluded = true
                break
              end
            end

            # Pick at random one extant branch for both sides of the n1th split, but only if none of the previously sampled extant branches is in the extant progeny of this side of the split.
            unless extProgId0AlreadyIncluded
              selectedId0 = extProgId0.sample
              # Do another fast id search to find the extant branch, whose id equals selectedId0.
              lower = -1
              upper = extantBranch.size
              idFound = false
              until idFound
                eb = lower+(upper-lower)/2
                if extantBranch[eb].id[1..-1].to_i < selectedId0[1..-1].to_i
                  lower = eb
                elsif extantBranch[eb].id[1..-1].to_i > selectedId0[1..-1].to_i
                  upper = eb
                elsif extantBranch[eb].id == selectedId0
                  idFound = true
                  sampledExtantBranch << extantBranch[eb]
                  next
                else
                  raise "Some problem occurred with the fast id search!" # This should never be called at all.
                end
              end
            end
            unless extProgId1AlreadyIncluded
              selectedId1 = extProgId1.sample
              # Do another fast id search to find the extant branch, whose id equals selectedId1.
              lower = -1
              upper = extantBranch.size
              idFound = false
              until idFound
                eb = lower+(upper-lower)/2
                if extantBranch[eb].id[1..-1].to_i < selectedId1[1..-1].to_i
                  lower = eb
                elsif extantBranch[eb].id[1..-1].to_i > selectedId1[1..-1].to_i
                  upper = eb
                elsif extantBranch[eb].id == selectedId1
                  idFound = true
                  sampledExtantBranch << extantBranch[eb]
                  next
                else
                  raise "Some problem occurred with the fast id search!" # This should never be called at all.
                end
              end # until idFound
            end # unless extProgId1AlreadyIncluded
          end # (sampledFGSpecies-1).times do |fgn|
        end # if @focusGroup != nil

        # Whether or not there's a focus group, now fill the array sampledExtantBranches until its size is n (the specified number of species to sample).
        if verbose
          if @focusGroup != nil
            print "\rSampling #{n-sampledFGSpecies-1} more species outside the focus group according to the diversified sampling scheme..."
          else
            print "\rSampling #{n-2} more species according to the diversified sampling scheme..."
          end
        end
        n1 = 0
        until sampledExtantBranch.size == n
          if n1 > n-2 # n - 2 is the number of internal branches in a rooted tree (not counting the root stem) with n tips. If n1 becomes larger than that, sortedBranch[n1] doesn't have a daughter.
            if @focusGroup != nil
              raise "When sampling only #{sampledFGSpecies} out of #{selectedFGOriginBranch.extantProgenyId.size} species from the focus group, the tree is not large enough to sample a total of #{n} out of #{@Np} species!"
            elsif @focusGroup == nil
              raise "There seems to be a problem with the sampling scheme!"
            end
          end
          # Determine the extant progeny of both sides of the n1th split.
          extProgId0 = []
          extProgId1 = []
          # Do another fast id search to find the branch, whose id equals sortedBranch[n1].daughterId[0].
          lower = -1
          upper = @branch.size
          d0Found = false
          until d0Found
            b = lower+(upper-lower)/2
            if @branch[b].id[1..-1].to_i < sortedBranch[n1].daughterId[0][1..-1].to_i
              lower = b
            elsif @branch[b].id[1..-1].to_i > sortedBranch[n1].daughterId[0][1..-1].to_i
              upper = b
            elsif @branch[b].id == sortedBranch[n1].daughterId[0]
              d0Found = true
              if @branch[b].extant
                extProgId0 = [@branch[b].id]
              else
                extProgId0 = @branch[b].extantProgenyId
              end
              next
            else
              raise "Some problem occurred with the fast id search!" # This should never be called at all.
            end
          end
          # Do another fast id search to find the branch, whose id equals sortedBranch[n1].daughterId[1].
          lower = -1
          upper = @branch.size
          d1Found = false
          until d1Found
            b = lower+(upper-lower)/2
            if @branch[b].id[1..-1].to_i < sortedBranch[n1].daughterId[1][1..-1].to_i
              lower = b
            elsif @branch[b].id[1..-1].to_i > sortedBranch[n1].daughterId[1][1..-1].to_i
              upper = b
            elsif @branch[b].id == sortedBranch[n1].daughterId[1]
              d1Found = true
              if @branch[b].extant
                extProgId1 = [@branch[b].id]
              else
                extProgId1 = @branch[b].extantProgenyId
              end
              next
            else
              raise "Some problem occurred with the fast id search!" # This should never be called at all.
            end
          end

          # See whether any of the previously sampled extant branches are in the extant progenies of one of the sides of the split.
          extProgId0AlreadyIncluded = false
          sampledExtantBranch.size.times do |seb|
            if extProgId0.include?(sampledExtantBranch[seb].id)
              extProgId0AlreadyIncluded = true
              break
            end
          end
          extProgId1AlreadyIncluded = false
          sampledExtantBranch.size.times do |seb|
            if extProgId1.include?(sampledExtantBranch[seb].id)
              extProgId1AlreadyIncluded = true
              break
            end
          end

          candidate0 = nil
          candidate1 = nil
          # Pick at random one extant branch for both sides of the n1th split, but only if none of the previously sampled extant branches is in the extant progeny of this side of the split.
          unless extProgId0AlreadyIncluded
            selectedId0 = extProgId0.sample
            # Do another fast id search to find the extant branch, whose id equals selectedId0.
            lower = -1
            upper = extantBranch.size
            idFound = false
            until idFound
              eb = lower+(upper-lower)/2
              if extantBranch[eb].id[1..-1].to_i < selectedId0[1..-1].to_i
                lower = eb
              elsif extantBranch[eb].id[1..-1].to_i > selectedId0[1..-1].to_i
                upper = eb
              elsif extantBranch[eb].id == selectedId0
                idFound = true
                candidate0 = extantBranch.delete_at(eb)
                next
              else
                raise "Some problem occurred with the fast id search!" # This should never be called at all.
              end
            end
            sampledExtantBranch << candidate0 unless candidate0.focusGroup
          end
          unless extProgId1AlreadyIncluded
            selectedId1 = extProgId1.sample

            # Do another fast id search to find the extant branch, whose id equals selectedId1.
            lower = -1
            upper = extantBranch.size
            idFound = false
            until idFound
              eb = lower+(upper-lower)/2
              if extantBranch[eb].id[1..-1].to_i < selectedId1[1..-1].to_i
                lower = eb
              elsif extantBranch[eb].id[1..-1].to_i > selectedId1[1..-1].to_i
                upper = eb
              elsif extantBranch[eb].id == selectedId1
                idFound = true
                candidate1 = extantBranch.delete_at(eb)
                next
              else
                raise "Some problem occurred with the fast id search!" # This should never be called at all.
              end
            end
            sampledExtantBranch << candidate1 unless candidate1.focusGroup
          end
          n1 += 1
        end

      elsif @samplingScheme == "semidiversified" # A new sampling scheme (see file `semidiversified sampling.pdf' in this directory).

        if verbose
          puts "\rModel:                                       "
          puts "Semidiversified sampling (unpublished, see 'semidiversified_sampling.pdf' in this directory)"
          puts
          puts "Parameters:"
          puts "n/m = #{n}/#{@Np}"
          puts
          if @focusGroup != nil
            puts "Conditions on focus group:"
            puts "#{minExtantFGSpecies-1} < number of extant species in focus group < #{maxExtantFGSpecies+1}"
            puts "#{minFGOrigin} < origin of focus group < #{maxFGOrigin}"
            puts "Number of species to be sampled from focus group: #{sampledFGSpecies}"
          end
        end

        if @focusGroup != nil
          # Search for potential focus groups
          print "Searching for potential focus groups..." if verbose
          candidateFGOriginBranches = []
          @branch.each do |b|
            if minExtantFGSpecies <= b.extantProgenyId.size and maxExtantFGSpecies >= b.extantProgenyId.size and minFGOrigin <= b.termination and maxFGOrigin >= b.termination
              candidateFGOriginBranches << b
            end
          end
          if candidateFGOriginBranches.size == 0
            if verbose
              puts "\rPotential focus groups in tree: 0 (stopping run...)"
              puts
            end
            raise RetryException.new(true), "The tree contains no clades that match the focus group criteria! Try again."
          elsif candidateFGOriginBranches.size == 1
            puts "\rPotential focus groups in tree: 1          " if verbose
          else
            puts "\rPotential focus groups in tree: #{candidateFGOriginBranches.size} (one of them is chosen at random)" if verbose
          end
          selectedFGOriginBranch = candidateFGOriginBranches.sample
          @focusGroupAge = selectedFGOriginBranch.termination
          @focusGroupNpFull = selectedFGOriginBranch.extantProgenyId.size
          if verbose
            puts
            print "\rAdding the focus group parameter to all branches..."
          end
          @branch.each do |b|
            if selectedFGOriginBranch.progenyId.include?(b.id)
              b.addFocusGroup(focusGroup = true)
            else
              b.addFocusGroup(focusGroup = false)
            end
          end          
        end

        print "\rPreparing sampling...                                                   " if verbose
        # Let's make sure that the array @branch is still sorted by the id numbers. This is important for the quick parent search below.
        branchesSorted = true
        (@branch.size-1).times do |b|
          if @branch[b].id[1..-1].to_i > @branch[b+1].id[1..-1].to_i
            branchesSorted = false
            break
          end
        end
        if branchesSorted == false
          # This is unlikely to happen anyway, and should be checked if it does.
          warn "WARNING: For some reason, branches are not sorted according to their id. Better check out why."
          until branchesSorted
            allSortedSoFar = true
            (@branch.size-1).times do |b|
              if @branch[b].id[1..-1].to_i > @branch[b+1].id[1..-1].to_i
                @branch[b],@branch[b+1] = @branch[b+1],@branch[b]
                allSortedSoFar = false
              end
            end
            branchesSorted = true if allSortedSoFar == true
          end
        end

        # The semidiversified scheme starts off just the same as for the 'diversified' sampling scheme.
        # The tree must first be reconstructed with all extant taxa.
        # Branches with endCause 'extinction' are deleted, and if they're the only daughter of the parent, then the parent's endCause also becomes 'extinction'.
        # This is done until no branches have endCause 'extinction' anymore.
        extinction = true
        while extinction == true
          extinction = false
          @branch.size.times do |b|
            if @branch[b].endCause == "extinction"
              # Fast parent search to find the branch whose id matches @branch[b].parentId. This is complicated by the fact that a branch can now be nil.
              lower = 0
              upper = @branch.size
              parentFound = false
              until parentFound == true
                # Use the fact that @branch is sorted according to id for a quick parent search.
                bb = lower+rand(upper-lower)
                unless @branch[bb] == nil
                  if @branch[bb].id[1..-1].to_i < @branch[b].parentId[1..-1].to_i
                    lower = bb
                  elsif @branch[bb].id[1..-1].to_i > @branch[b].parentId[1..-1].to_i
                    upper = bb
                  elsif @branch[bb].id == @branch[b].parentId
                    parentFound = true
                    if @branch[bb].daughterId.size == 2
                      daughterId = @branch[bb].daughterId.dup
                      daughterId.delete(@branch[b].id)
                      @branch[bb].updateDaughterId(daughterId)
                    elsif @branch[bb].daughterId.size == 1
                      @branch[bb].updateDaughterId([])
                      @branch[bb].updateEndCause("extinction")
                    else
                      raise "branch[#{bb}] has #{@branch[bb].daughterId.size} daughters!" # should never be called at all.
                    end
                    @branch[bb].updateFossil(@branch[bb].fossil+@branch[b].fossil) if @fossilRecordAdded
                    next
                  else
                    raise "Some problem occurred with the fast parent search!" # This should never be called at all.
                  end
                end
              end
              extinction = true
              @branch[b] = nil
            end
          end
          @branch.compact!
        end

        # Delete degree one vertices.
        unmerged = true
        while unmerged == true
          unmerged = false
          @branch.size.times do |b|
            unless @branch[b] == nil
              if @branch[b].daughterId.size == 1
                lower = 0
                upper = @branch.size
                daughterFound = false
                until daughterFound
                  # Use the fact that @branch is sorted according to id for a quick daughter search.
                  bb = lower+rand(upper-lower)
                  unless @branch[bb] == nil
                    if @branch[bb].id[1..-1].to_i < @branch[b].daughterId[0][1..-1].to_i
                      lower = bb
                    elsif @branch[bb].id[1..-1].to_i > @branch[b].daughterId[0][1..-1].to_i
                      upper = bb
                    elsif @branch[bb].id == @branch[b].daughterId[0]
                      daughterFound = true
                      @branch[b].updateTermination(@branch[bb].termination)
                      @branch[b].updateDaughterId(@branch[bb].daughterId)
                      @branch[b].updateEndCause(@branch[bb].endCause)
                      @branch[b].updateExtant(@branch[bb].extant)
                      @branch[b].updateRate((@branch[b].duration*@branch[b].rate+@branch[bb].duration*@branch[bb].rate)/(@branch[b].duration+@branch[bb].duration)) if @branchRatesAssigned
                      @branch[b].updateDuration(@branch[b].origin - @branch[b].termination)
                      speciesIdAry = @branch[b].speciesId.split("/")
                      speciesIdAry1 = @branch[bb].speciesId.split("/")
                      speciesIdAry1.each do |s|
                        speciesIdAry << s unless speciesIdAry.include?(s)
                      end
                      speciesIdStr = ""
                      (speciesIdAry.size-1).times do |sid|
                        speciesIdStr << "#{speciesIdAry[sid]}/"
                      end
                      speciesIdStr << "#{speciesIdAry[-1]}"
                      @branch[b].updateSpeciesId(speciesIdStr)
                      @branch[b].updateFossil(@branch[b].fossil + @branch[bb].fossil) if @fossilRecordAdded
                      # Find the one or two daughters of @branch[bb], and update their parentIds to the id of @branch[b], as @branch[bb] is about to be deleted.
                      # If @branch[bb] has no daughters, no parentIds need to be updated.
                      if @branch[bb].daughterId.size == 1
                        # Do a fast granddaughter search.
                        lower2 = 0
                        upper2 = @branch.size
                        granddaughterFound = false
                        until granddaughterFound
                          bbb = lower2+rand(upper2-lower2)
                          unless @branch[bbb] == nil
                            if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[0][1..-1].to_i
                              lower2 = bbb
                            elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[0][1..-1].to_i
                              upper2 = bbb
                            elsif @branch[bbb].id == @branch[bb].daughterId[0]
                              unless @branch[bbb].parentId == @branch[bb].id
                                raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[0]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                              end
                              granddaughterFound = true
                              @branch[bbb].updateParentId(@branch[b].id)
                              next
                            else
                              raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                            end
                          end
                        end
                      elsif @branch[bb].daughterId.size == 2 and @branch[bb].daughterId != ["none","none"]
                        # Do a fast search for the first granddaughter.
                        lower2 = 0
                        upper2 = @branch.size
                        granddaughter0Found = false
                        until granddaughter0Found
                          bbb = lower2+rand(upper2-lower2)
                          unless @branch[bbb] == nil
                            if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[0][1..-1].to_i
                              lower2 = bbb
                            elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[0][1..-1].to_i
                              upper2 = bbb
                            elsif @branch[bbb].id == @branch[bb].daughterId[0]
                              unless @branch[bbb].parentId == @branch[bb].id
                                raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[0]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                              end
                              granddaughter0Found = true
                              @branch[bbb].updateParentId(@branch[b].id)
                              next
                            else
                              raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                            end
                          end
                        end
                        # Do a fast search for the second granddaughter.
                        lower2 = 0
                        upper2 = @branch.size
                        granddaughter1Found = false
                        until granddaughter1Found
                          bbb = lower2+rand(upper2-lower2)
                          unless @branch[bbb] == nil
                            if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[1][1..-1].to_i
                              lower2 = bbb
                            elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[1][1..-1].to_i
                              upper2 = bbb
                            elsif @branch[bbb].id == @branch[bb].daughterId[1]
                              unless @branch[bbb].parentId == @branch[bb].id
                                raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[1]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                              end
                              granddaughter1Found = true
                              @branch[bbb].updateParentId(@branch[b].id)
                              next
                            else
                              raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                            end
                          end
                        end # until granddaughter1Found
                      end # if @branch[bb].daughterId.size == 1, elsif @branch[bb].daughterId.size == 2
                      @branch[bb] = nil
                      next
                    else
                      raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
                    end
                  end
                end
                unmerged = true
              end
            end
          end
          @branch.compact!
        end

        # Now we should have a reconstructed tree of all extant lineages.
        # Make sure all branches have exactly two daughters
        @branch.each do |b|
          raise "Branch #{b.id} has #{b.daughterId.size} daughters!" unless b.daughterId.size == 2
        end

        # Deleted branches must be removed from progenyId and extantProgenyId.
        @branch.each do |b|
          newProgenyId = []
          newExtantProgenyId = []
          progenyId = b.progenyId
          progenyId.each do |p|
            # Do another fast id search to find the branch, whose id equals p.
            lower = -1
            upper = @branch.size
            progenyFound = false
            until progenyFound
              # Use the fact that @branch is sorted according to id for a quick id search.
              if upper-lower == 1
                progenyFound = true
              else
                bb = lower+(upper-lower)/2
                if @branch[bb].id[1..-1].to_i < p[1..-1].to_i
                  lower = bb
                elsif @branch[bb].id[1..-1].to_i > p[1..-1].to_i
                  upper = bb
                elsif @branch[bb].id == p
                  progenyFound = true
                  newProgenyId << p
                  newExtantProgenyId << p if @branch[bb].extant
                  next
                else
                  raise "Some problem occurred with the fast progeny search!" # This should never be called at all.
                end
              end
            end
          end
          b.updateProgenyId(newProgenyId)
          b.updateExtantProgenyId(newExtantProgenyId)
        end

        # The array of all extant branches needs to be filled again, as some have been renamed in the previous reconstruction.
        extantBranch = []
        @branch.each do |b|
          extantBranch << b if b.extant
        end
        raise "extantBranch.size is not equal to @Np" unless extantBranch.size == @Np # should not be called at all.

        if verbose
          if @focusGroup != nil
            print "\rRandomly sampling one species from each side of the root, one of them from the focus group..."
          else
            print "\rRandomly sampling one species from each side of the root..."
          end
        end
        # Start with sampling two extant species, each from one side of the root to make sure the root node is included.
        extProgId0 = []
        extProgId1 = []
        if @rootSplit
          b0Found = false
          b = 0
          until b0Found
            if @branch[b].id == "b0"
              if @branch[b].extant
                extProgId0 = [@branch[b].id]
              else
                extProgId0 = @branch[b].extantProgenyId
              end
              b0Found = true
            end
            b += 1
          end
          b1Found = false
          b = 0
          until b1Found
            if @branch[b].id == "b1"
              if @branch[b].extant
                extProgId1 = [@branch[b].id]
              else
                extProgId1 = @branch[b].extantProgenyId
              end
              b1Found = true
            end
            b += 1
          end
        else
          b0Id = nil
          b1Id = nil
          b0Found = false
          b = 0
          until b0Found
            if @branch[b].id == "b0"
              b0Id = @branch[b].daughterId[0]
              b1Id = @branch[b].daughterId[1]
              b0Found = true
            end
            b += 1
          end
          newB0Found = false
          b = 0
          until newB0Found
            if @branch[b].id == b0Id
              if @branch[b].extant
                extProgId0 = [@branch[b].id]
              else
                extProgId0 = @branch[b].extantProgenyId
              end
              newB0Found = true
            end
            b += 1
          end
          newB1Found = false
          b = 0
          until newB1Found
            if @branch[b].id == b1Id
              if @branch[b].extant
                extProgId1 = [@branch[b].id]
              else
                extProgId1 = @branch[b].extantProgenyId
              end
              newB1Found = true
            end
            b += 1
          end
        end

        # If we have focus group, make sure that one of the selected ids is drawn from it. This is because the focus group must necessarily be in one of
        # the partial trees. To do that, first find out in which part of the tree the focus group is.
        if @focusGroup != nil
          # We divide the array extantBranch into extantFGBranch and extantNonFGBranch, that will facilitate things later.
          extantFGBranch = []
          extantNonFGBranch = []
          extantBranch.each do |eb|
            if eb.focusGroup
              extantFGBranch << eb
            else
              extantNonFGBranch << eb
            end
          end
          # If the first branch of array extantFGBranch is included in extProgId0, then the entire focus group must be in it.
          if extProgId0.include?(extantFGBranch[0].id)
            selectedId0 = extantFGBranch.sample.id
            selectedId1 = extProgId1.sample
          elsif extProgId1.include?(extantFGBranch[0].id)
            selectedId0 = extProgId0.sample
            selectedId1 = extantFGBranch.sample.id
          else
            raise "One extant focus group branch is not included in both progenies of the two most basal branches!"
          end
        elsif @focusGroup == nil
          selectedId0 = extProgId0.sample
          selectedId1 = extProgId1.sample
        end

        # Do another fast id search to find the extant branch, whose id equals selectedId0 (extantBranch must also be sorted according to id).
        lower = -1
        upper = extantBranch.size
        idFound = false
        until idFound
          eb = lower+(upper-lower)/2
          if extantBranch[eb].id[1..-1].to_i < selectedId0[1..-1].to_i
            lower = eb
          elsif extantBranch[eb].id[1..-1].to_i > selectedId0[1..-1].to_i
            upper = eb
          elsif extantBranch[eb].id == selectedId0
            idFound = true
            sampledExtantBranch << extantBranch.delete_at(eb)
            next
          else
            raise "Some problem occurred with the fast id search!" # This should never be called at all.
          end
        end
        # Do another fast id search to find the extant branch, whose id equals selectedId1.
        lower = -1
        upper = extantBranch.size
        idFound = false
        until idFound
          eb = lower+(upper-lower)/2
          if extantBranch[eb].id[1..-1].to_i < selectedId1[1..-1].to_i
            lower = eb
          elsif extantBranch[eb].id[1..-1].to_i > selectedId1[1..-1].to_i
            upper = eb
          elsif extantBranch[eb].id == selectedId1
            idFound = true
            sampledExtantBranch << extantBranch.delete_at(eb)
            next
          else
            raise "Some problem occurred with the fast id search!" # This should never be called at all.
          end
        end
        
        # Create an array for all branches that are not extant.
        extinctBranch = []
        @branch.each do |b|
          extinctBranch << b
        end
        extinctBranch.delete_if { |b| b.extant }

        if @focusGroup != nil
          print "\rSampling #{sampledFGSpecies-1} more species from the focus group according to the semidiversified sampling scheme..." if verbose
          # Create an array for all focus group branches that are not extant.
          extinctFGBranch = []
          extinctNonFGBranch = []
          extinctBranch.each do |ecb|
            extinctFGBranch << ecb
            extinctNonFGBranch << ecb
          end
          extinctFGBranch.delete_if { |b| b.focusGroup == false }
          extinctNonFGBranch.delete_if { |b| b.focusGroup == true }
        end

        # First apply the semidiversified sampling scheme to the focus group if we have one. Start of by making sure that the root of the focus group is sampled.
        if @focusGroup != nil
          # Find the two oldest members of array selectedFGOriginBranch.
          oldestFGBranches = [extinctFGBranch[0],extinctFGBranch[1]]
          oldestFGBranches.sort! { |a,b| b.origin <=> a.origin }
          fGRootBranchesSampled = 0
          extinctFGBranch.each do |ecb|
            if ecb.origin > oldestFGBranches[1].origin
              oldestFGBranches[1] = ecb
              oldestFGBranches.sort! { |a,b| b.origin <=> a.origin }
            end
          end

          until sampledExtantBranch.size > sampledFGSpecies # Because one species from the other side of the root is in there already.
            # Pick a node at random, with uniform probability (as nodes are not explicitely used here, we consider their parental branches instead, and use their termination date as the node age).
            if fGRootBranchesSampled < 2
              selectedBranch = oldestFGBranches[fGRootBranchesSampled] # Make sure the root of the focus group is sampled.
              fGRootBranchesSampled += 1
            else
              selectedBranch = extinctFGBranch.sample
            end
            # Determine the probability p(t) that the selected node is sampled.
            tFGRoot = selectedFGOriginBranch.termination
            if fGRootBranchesSampled < 2
              pT = 1.0
            else
              pT = selectedBranch.termination/tFGRoot
            end
            # sample this node with probability pT
            if rand < pT
              # Determine the extant progeny of both sides of the n1th split.
              extProgId0 = []
              # Do another fast id search to find the branch, whose id equals selectedBranch.daughterId[0].
              lower = 0
              upper = @branch.size
              d0Found = false
              until d0Found
                b = lower+(upper-lower)/2
                if @branch[b].id[1..-1].to_i < selectedBranch.daughterId[0][1..-1].to_i
                  lower = b
                elsif @branch[b].id[1..-1].to_i > selectedBranch.daughterId[0][1..-1].to_i
                  upper = b
                elsif @branch[b].id == selectedBranch.daughterId[0]
                  d0Found = true
                  if @branch[b].extant
                    extProgId0 = [@branch[b].id]
                  else
                    extProgId0 = @branch[b].extantProgenyId
                  end
                  next
                else
                  raise "Some problem occurred with the fast id search!" # This should never be called at all.
                end
              end
              extProgId1 = []
              # Do another fast id search to find the branch, whose id equals selectedBranch.daughterId[1].
              lower = 0
              upper = @branch.size
              d1Found = false
              until d1Found
                b = lower+(upper-lower)/2
                if @branch[b].id[1..-1].to_i < selectedBranch.daughterId[1][1..-1].to_i
                  lower = b
                elsif @branch[b].id[1..-1].to_i > selectedBranch.daughterId[1][1..-1].to_i
                  upper = b
                elsif @branch[b].id == selectedBranch.daughterId[1]
                  d1Found = true
                  if @branch[b].extant
                    extProgId1 = [@branch[b].id]
                  else
                    extProgId1 = @branch[b].extantProgenyId
                  end
                  next
                else
                  raise "Some problem occurred with the fast id search!" # This should never be called at all.
                end
              end

              # See whether any of the previously sampled extant branches are in the extant progenies of one of the sides of the split.
              extProgId0AlreadyIncluded = false
              sampledExtantBranch.size.times do |seb|
                if extProgId0.include?(sampledExtantBranch[seb].id)
                  extProgId0AlreadyIncluded = true
                  break
                end
              end
              extProgId1AlreadyIncluded = false
              sampledExtantBranch.size.times do |seb|
                if extProgId1.include?(sampledExtantBranch[seb].id)
                  extProgId1AlreadyIncluded = true
                  break
                end
              end

              # Pick at random one extant branch for both sides of the n1th split, but only if none of the previously sampled extant branches is in the extant progeny of this side of the split.
              unless extProgId0AlreadyIncluded
                selectedId0 = extProgId0[rand(extProgId0.size)]
                # Do another fast id search to find the branch, whose id equals selectedId0.
                lower = -1
                upper = extantBranch.size
                idFound = false
                until idFound
                  eb = lower+(upper-lower)/2
                  if extantBranch[eb].id[1..-1].to_i < selectedId0[1..-1].to_i
                    lower = eb
                  elsif extantBranch[eb].id[1..-1].to_i > selectedId0[1..-1].to_i
                    upper = eb
                  elsif extantBranch[eb].id == selectedId0
                    idFound = true
                    sampledExtantBranch << extantBranch[eb]
                    next
                  else
                    raise "Some problem occurred with the fast id search!" # This should never be called at all.
                  end
                end
              end
              unless extProgId1AlreadyIncluded
                selectedId1 = extProgId1[rand(extProgId1.size)]
                # Do another fast id search to find the branch, whose id equals selectedId0.
                lower = -1
                upper = extantBranch.size
                idFound = false
                until idFound
                  eb = lower+(upper-lower)/2
                  if extantBranch[eb].id[1..-1].to_i < selectedId1[1..-1].to_i
                    lower = eb
                  elsif extantBranch[eb].id[1..-1].to_i > selectedId1[1..-1].to_i
                    upper = eb
                  elsif extantBranch[eb].id == selectedId1
                    idFound = true
                    sampledExtantBranch << extantBranch[eb]
                    next
                  else
                    raise "Some problem occurred with the fast id search!" # This should never be called at all.
                  end
                end
              end
            end # if rand < pT
          end # until sampledExtantBranch.size > sampledFGSpecies
          while sampledExtantBranch.size > sampledFGSpecies+1
            sampledExtantBranch = sampledExtantBranch[0..-2]
          end
        end

        if verbose
          if @focusGroup != nil
            print "\rSampling #{n-sampledFGSpecies-1} more species outside the focus group according to the semidiversified sampling scheme..."
          else
            print "\rSampling #{n-2} more species according to the semidiversified sampling scheme..."
          end
        end
        # Now apply the semidiversified sampling scheme to the entire tree (without the focus group if we have one).
        until sampledExtantBranch.size > n-1

          # Pick a node at random, with uniform probability (as nodes are not explicitely used here, we consider their parental branches instead, and use their termination date as the node age).
          if @focusGroup != nil
            selectedBranch = extinctNonFGBranch.sample
          else
            selectedBranch = extinctBranch.sample
          end
          # Determine the probability p(t) that the selected node is sampled.
          if @rootSplit
            tRoot = @posteriorTreeOrigin
          else
            b = 0
            b0Found = false
            until b0Found
              if @branch[b].id == "b0"
                b0Found = true
                tRoot = @branch[b].termination
              end
              b += 1
            end  
          end
          pT = selectedBranch.termination/tRoot
          # sample this node with probability pT
          if rand < pT
            # Determine the extant progeny of both sides of the n1th split.
            extProgId0 = []
            # Do another fast id search to find the branch, whose id equals selectedBranch.daughterId[0].
            lower = 0
            upper = @branch.size
            d0Found = false
            until d0Found
              b = lower+(upper-lower)/2
              if @branch[b].id[1..-1].to_i < selectedBranch.daughterId[0][1..-1].to_i
                lower = b
              elsif @branch[b].id[1..-1].to_i > selectedBranch.daughterId[0][1..-1].to_i
                upper = b
              elsif @branch[b].id == selectedBranch.daughterId[0]
                d0Found = true
                if @branch[b].extant
                  extProgId0 = [@branch[b].id]
                else
                  extProgId0 = @branch[b].extantProgenyId
                end
                next
              else
                raise "Some problem occurred with the fast id search!" # This should never be called at all.
              end
            end
            extProgId1 = []
            # Do another fast id search to find the branch, whose id equals selectedBranch.daughterId[1].
            lower = 0
            upper = @branch.size
            d1Found = false
            until d1Found
              b = lower+(upper-lower)/2
              if @branch[b].id[1..-1].to_i < selectedBranch.daughterId[1][1..-1].to_i
                lower = b
              elsif @branch[b].id[1..-1].to_i > selectedBranch.daughterId[1][1..-1].to_i
                upper = b
              elsif @branch[b].id == selectedBranch.daughterId[1]
                d1Found = true
                if @branch[b].extant
                  extProgId1 = [@branch[b].id]
                else
                  extProgId1 = @branch[b].extantProgenyId
                end
                next
              else
                raise "Some problem occurred with the fast id search!" # This should never be called at all.
              end
            end

            # See whether any of the previously sampled extant branches are in the extant progenies of one of the sides of the split.
            extProgId0AlreadyIncluded = false
            sampledExtantBranch.size.times do |seb|
              if extProgId0.include?(sampledExtantBranch[seb].id)
                extProgId0AlreadyIncluded = true
                break
              end
            end
            extProgId1AlreadyIncluded = false
            sampledExtantBranch.size.times do |seb|
              if extProgId1.include?(sampledExtantBranch[seb].id)
                extProgId1AlreadyIncluded = true
                break
              end
            end

            # Pick at random one extant branch for both sides of the n1th split, but only if none of the previously sampled extant branches is in the extant progeny of this side of the split.
            unless extProgId0AlreadyIncluded
              selectedId0 = extProgId0[rand(extProgId0.size)]
              # Do another fast id search to find the branch, whose id equals selectedId0.
              lower = -1
              upper = extantBranch.size
              idFound = false
              until idFound
                eb = lower+(upper-lower)/2
                if extantBranch[eb].id[1..-1].to_i < selectedId0[1..-1].to_i
                  lower = eb
                elsif extantBranch[eb].id[1..-1].to_i > selectedId0[1..-1].to_i
                  upper = eb
                elsif extantBranch[eb].id == selectedId0
                  idFound = true
                  sampledExtantBranch << extantBranch[eb]
                  next
                else
                  raise "Some problem occurred with the fast id search!" # This should never be called at all.
                end
              end
            end
            unless extProgId1AlreadyIncluded
              selectedId1 = extProgId1[rand(extProgId1.size)]
              # Do another fast id search to find the branch, whose id equals selectedId0.
              lower = -1
              upper = extantBranch.size
              idFound = false
              until idFound
                eb = lower+(upper-lower)/2
                if extantBranch[eb].id[1..-1].to_i < selectedId1[1..-1].to_i
                  lower = eb
                elsif extantBranch[eb].id[1..-1].to_i > selectedId1[1..-1].to_i
                  upper = eb
                elsif extantBranch[eb].id == selectedId1
                  idFound = true
                  sampledExtantBranch << extantBranch[eb]
                  next
                else
                  raise "Some problem occurred with the fast id search!" # This should never be called at all.
                end
              end
            end

          end # if rand < pT
        end # while sampledExtantBranch < n
        while sampledExtantBranch.size > n
          sampledExtantBranch = sampledExtantBranch[0..-2]
        end
      else
        raise "Sampling scheme could not be recognized. It should be either 'all', 'random', 'diversified', or 'semidiversified'!"
      end # if sampling == "all", elsif sampling == "random", elsif sampling == "diversified", elsif sampling == "semidiversified", else...
      
    end # if n > @nP, elsif n == 0, elsif n == 1, else...

    if verbose
      print "\r                                                                                                                    "
      print "\rFinalizing reconstruction..." 
    end
    @branch.each do |b|
      sampled = false
      # See whether one of the branches in sampledExtantBranch has the id of this branch.
      sampledExtantBranch.each do |seb|
        if seb.id == b.id
          sampled = true
          break
        end
      end
      # The endCause of all nonsampled extant branches becomes 'extinction' instead of 'present', so that they can be treated as if they had gone extinct infinitesimally short before the present.
      b.updateEndCause("extinction") if b.extant and sampled == false
    end

    # Branches with endCause 'extinction' are deleted, and if they're the only daughter of the parent, then the parent's endCause also becomes 'extinction'.
    # This is done until no branches have endCause 'extinction' anymore.
    extinction = true
    while extinction == true
      extinction = false
      @branch.size.times do |b|
        if @branch[b].endCause == "extinction"
          # Fast parent search to find the branch that has the id that matches @branch[b].parentId. This is complicated by the fact that a branch can now be nil.
          lower = 0
          upper = @branch.size
          parentFound = false
          until parentFound == true
            # Use the fact that @branch is sorted according to id for a quick parent search.
            bb = lower+rand(upper-lower)
            unless @branch[bb] == nil
              if @branch[bb].id[1..-1].to_i < @branch[b].parentId[1..-1].to_i
                lower = bb
              elsif @branch[bb].id[1..-1].to_i > @branch[b].parentId[1..-1].to_i
                upper = bb
              elsif @branch[bb].id == @branch[b].parentId
                parentFound = true
                if @branch[bb].daughterId.size == 2
                  daughterId = @branch[bb].daughterId.dup
                  daughterId.delete(@branch[b].id)
                  @branch[bb].updateDaughterId(daughterId)
                elsif @branch[bb].daughterId.size == 1
                  @branch[bb].updateDaughterId([])
                  @branch[bb].updateEndCause("extinction")
                else
                  raise "branch[#{bb}] has #{@branch[bb].daughterId.size} daughters!" # should never be called at all.
                end
                @branch[bb].updateFossil(@branch[bb].fossil+@branch[b].fossil) if @fossilRecordAdded
                next
              else
                raise "Some problem occurred with the fast parent search!" # This should never be called at all.
              end
            end
          end
          extinction = true
          @branch[b] = nil
        end
      end
      @branch.compact!
    end

    # Delete degree one vertices.
    unmerged = true
    while unmerged == true
      unmerged = false
      @branch.size.times do |b|
        unless @branch[b] == nil
          if @branch[b].daughterId.size == 1
            lower = 0
            upper = @branch.size
            daughterFound = false
            until daughterFound
              # Use the fact that @branch is sorted according to id for a quick daughter search.
              bb = lower+rand(upper-lower)
              unless @branch[bb] == nil
                if @branch[bb].id[1..-1].to_i < @branch[b].daughterId[0][1..-1].to_i
                  lower = bb
                elsif @branch[bb].id[1..-1].to_i > @branch[b].daughterId[0][1..-1].to_i
                  upper = bb
                elsif @branch[bb].id == @branch[b].daughterId[0]
                  daughterFound = true
                  @branch[b].updateTermination(@branch[bb].termination)
                  @branch[b].updateDaughterId(@branch[bb].daughterId)
                  @branch[b].updateEndCause(@branch[bb].endCause)
                  @branch[b].updateExtant(@branch[bb].extant)
                  @branch[b].updateRate((@branch[b].duration*@branch[b].rate+@branch[bb].duration*@branch[bb].rate)/(@branch[b].duration+@branch[bb].duration)) if @branchRatesAssigned
                  @branch[b].updateDuration(@branch[b].origin - @branch[b].termination)
                  speciesIdAry = @branch[b].speciesId.split("/")
                  speciesIdAry1 = @branch[bb].speciesId.split("/")
                  speciesIdAry1.each do |s|
                    speciesIdAry << s unless speciesIdAry.include?(s)
                  end
                  speciesIdStr = ""
                  (speciesIdAry.size-1).times do |sid|
                    speciesIdStr << "#{speciesIdAry[sid]}/"
                  end
                  speciesIdStr << "#{speciesIdAry[-1]}"
                  @branch[b].updateSpeciesId(speciesIdStr)
                  @branch[b].updateFossil(@branch[b].fossil + @branch[bb].fossil) if @fossilRecordAdded
                  # Find the one or two daughters of @branch[bb], and update their parentIds to the id of @branch[b], as @branch[bb] is about to be deleted.
                  # If @branch[bb] has no daughters, no parentIds need to be updated.
                  if @branch[bb].daughterId.size == 1
                    # Do a fast granddaughter search.
                    lower2 = 0
                    upper2 = @branch.size
                    granddaughterFound = false
                    until granddaughterFound
                      bbb = lower2+rand(upper2-lower2)
                      unless @branch[bbb] == nil
                        if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[0][1..-1].to_i
                          lower2 = bbb
                        elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[0][1..-1].to_i
                          upper2 = bbb
                        elsif @branch[bbb].id == @branch[bb].daughterId[0]
                          unless @branch[bbb].parentId == @branch[bb].id
                            raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[0]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                          end
                          granddaughterFound = true
                          @branch[bbb].updateParentId(@branch[b].id)
                          next
                        else
                          raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                        end
                      end
                    end
                  elsif @branch[bb].daughterId.size == 2 and @branch[bb].daughterId != ["none","none"]
                    # Do a fast search for the first granddaughter.
                    lower2 = 0
                    upper2 = @branch.size
                    granddaughter0Found = false
                    until granddaughter0Found
                      bbb = lower2+rand(upper2-lower2)
                      unless @branch[bbb] == nil
                        if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[0][1..-1].to_i
                          lower2 = bbb
                        elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[0][1..-1].to_i
                          upper2 = bbb
                        elsif @branch[bbb].id == @branch[bb].daughterId[0]
                          unless @branch[bbb].parentId == @branch[bb].id
                            raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[0]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                          end
                          granddaughter0Found = true
                          @branch[bbb].updateParentId(@branch[b].id)
                          next
                        else
                          raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                        end
                      end
                    end
                    # Do a fast search for the second granddaughter.
                    lower2 = 0
                    upper2 = @branch.size
                    granddaughter1Found = false
                    until granddaughter1Found
                      bbb = lower2+rand(upper2-lower2)
                      unless @branch[bbb] == nil
                        if @branch[bbb].id[1..-1].to_i < @branch[bb].daughterId[1][1..-1].to_i
                          lower2 = bbb
                        elsif @branch[bbb].id[1..-1].to_i > @branch[bb].daughterId[1][1..-1].to_i
                          upper2 = bbb
                        elsif @branch[bbb].id == @branch[bb].daughterId[1]
                          unless @branch[bbb].parentId == @branch[bb].id
                            raise "The daughter id of branch #{@branch[bb].id} is #{@branch[bb].daughterId[1]}, but the parent id of branch #{@branch[bbb].id} is #{@branch[bbb].parentId}!" # Should not be called at all.
                          end
                          granddaughter1Found = true
                          @branch[bbb].updateParentId(@branch[b].id)
                          next
                        else
                          raise "Some problem occurred with the fast granddaughter search!" # This should never be called at all.
                        end
                      end
                    end # until granddaughter1Found
                  end # if @branch[bb].daughterId.size == 1, elsif @branch[bb].daughterId.size == 2
                  @branch[bb] = nil
                  next
                else
                  raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
                end
              end              
            end
            unmerged = true
          end
        end
      end
      @branch.compact!
    end

    # If rootSplit == false, delete the root branch and update @treeOrigin
    unless @rootSplit
      @branch.size.times do |b|
        @branch[b] = nil if @branch[b].origin == @posteriorTreeOrigin
      end
      @branch.compact!
      origin = 0
      @branch.each do |b|
        origin = b.origin if b.origin > origin
      end
      @treeOrigin = origin
    end
    
    # Update @Np and add @phylogeneticDiversity
    n = 0
    pd = 0
    @branch.each do |b|
      n += 1 if b.extant
      pd += b.duration
    end
    @Np = n
    @phylogeneticDiversity = pd

    # Add the sampled extant progeny to each branch. At the same time, recalculate the extant diversity of each branch.
    newExtantBranch = [] # The array extantBranch that has been used above should not be used anymore after degree-one vertices have been deleted.
    @branch.each {|b| newExtantBranch << b if b.extant}
    @branch.each do |b|
      if b.extant
        b.updateExtantProgenyId([])
        b.updateExtantDiversity(1)
      else
        newExtantProgenyId = []
        newExtantBranch.each do |neb|
          raise "Extant branch #{neb.id} is not included in the array of all branches!" unless @branch.include?(neb)
          newExtantProgenyId << neb.id if b.progenyId.include?(neb.id)
        end
        b.updateExtantProgenyId(newExtantProgenyId)
        b.updateExtantDiversity(newExtantProgenyId.size)
      end
    end

    # Calculate the gamma statistic of Pybus & Harvey (2000)
    n = @Np
    if n > 2 # gamma is only defined if n > 2 (see Pybus & Harvey 2000)
      # the array g holds all internode distances. It is initiated as [0,0] cause Pybus & Harvey only use g_2 to g_n.
      g = [0,0] 
      nodeAge = []
      @branch.each do |b|
        nodeAge << b.origin
      end
      nodeAge.sort!.uniq!.reverse!
      raise "The number of node ages (#{nodeAge.size}) does not match the number of extant species-1 (#{@Np-1})!" unless nodeAge.size == @Np-1
      (nodeAge.size-1).times do |na|
        g << nodeAge[na]-nodeAge[na+1]
      end
      g << nodeAge[-1]
      tau = 0
      2.upto(n) do |j|
        tau += j*g[j]
      end
      sumblock = 0
      2.upto(n-1) do |i|
        2.upto(i) do |k|
          sumblock += k*g[k]
        end
      end
      @gamma = ((1/(n-2).to_f)*sumblock - tau/2.0)/(tau*Math.sqrt( 1/(12.0*(n-2)) ) ) # this may need to be tested with pure birth runs, which should result in a mean gamma of 0.\
    else
      @gamma = nil
    end

    # Assign alignment positions with or without focus group.
    print "\rAssigning alignment positions...                                                                               " if verbose
    alignmentPositionsAssigned = false # Otherwise the method assignAlignmentPositions will raise an error.
    verboseHere = verbose
    self.assignAlignmentPositions(verbose = false)

    if verboseHere
      puts "\rSuccessfully reconstructed the tree after sampling #{@Np} species according to the #{@samplingScheme} sampling scheme."
      endTime = Time.now
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end

    @treeReconstructed = true
    self
  end
  
  def introduceAgeError(meanAgeChangePerStep = 0.1, targetedMeanAgeError = 0.3, fileName = "changes.txt", verbose = true)

    # Feedback.
    if verbose
      startTime = Time.now
      puts
      puts "---------------------------------------------phylsim.rb | Age error introduction---------------------------------------------"
      puts
    end

    # Set the mean age error to 0.
    meanAgeError = 0

    # Memorize the mean age before error introduction.
    meanAgeBeforeErrorIntroduction = self.meanAge

    # Make sure the selected mean age change per step is a reasonable value.
    unless meanAgeChangePerStep > 0 and meanAgeChangePerStep < 0.5
      raise "The specified mean age change should be between 0 and 0.5 (the change is measured relative to branch lengths)!"
    end

    # Feedback.
    if verbose
      puts "Mean age change per step: #{meanAgeChangePerStep}"
      puts "Targeted mean age error: #{targetedMeanAgeError}"
      puts
    end

    # Initiate a loop count.
    loopcount = 0

    # Start a loop that only ends when the mean age error is equal to or greater than the targeted mean age error.
    while meanAgeError < targetedMeanAgeError

      # Increase the loop count.
      loopcount += 1

      # Pick a branch at random, but use only internal branches (those that end in speciation events).
      selectedBranch = nil
      until selectedBranch
        selectedBranch = @branch.sample
        selectedBranch = nil unless selectedBranch.endCause == "speciation"
      end

      # Determine the age change at this step.
      ageChange = 1 + rand*2*meanAgeChangePerStep

      # Determine the direction of the age change.
      if self.meanAge < meanAgeBeforeErrorIntroduction
        ageChangeDirection = "older"
      else
        ageChangeDirection = "younger"
      end

      # Change the duration of the branch by changing its termination.
      if ageChangeDirection == "older"
        newBranchDuration = (selectedBranch.duration/ageChange)*((selectedBranch.duration+selectedBranch.termination)/((selectedBranch.duration/ageChange)+selectedBranch.termination))
      else
        newBranchDuration = selectedBranch.duration*ageChange*((selectedBranch.duration+selectedBranch.termination)/((selectedBranch.duration*ageChange)+selectedBranch.termination))
      end
      newTermination = selectedBranch.origin-newBranchDuration
      ageChangeForProgeny = newTermination/selectedBranch.termination
      selectedBranch.updateTerminationAfterErrorIntroduction(newTermination)

      # Adjust the progeny of the selected branch accordingly.
      @branch.each do |b|
        if selectedBranch.progenyId.include?(b.id)
          b.updateOriginAfterErrorIntroduction(b.origin*ageChangeForProgeny)
          b.updateTerminationAfterErrorIntroduction(b.termination*ageChangeForProgeny)
        end
      end

      # Determine the current mean age error.
      meanAgeError = self.meanAgeError

      # Feedback.
      print "\rCurrent mean age error: #{meanAgeError.round(5)}..."
    end

    # Make sure that the termination for each branch is younger than the origin.
    @branch.each do |b|
      raise "The branch termination is older than the origin for branch #{b.id}: termination: #{b.termination}; origin: #{b.origin}!" if b.termination > b.origin
    end

    if verbose
      puts "\rMean age before error introduction: #{meanAgeBeforeErrorIntroduction.round(5)}"
      puts "Mean age after error introduction: #{self.meanAge.round(5)}"
      puts "Mean age error: #{meanAgeError.round(5)} was obtained after a total of #{loopcount} changes."
      puts "Number of underestimated branches: #{self.numberOfUnderestimatedBranches}."
      puts "Number of overestimated branches: #{self.numberOfOverestimatedBranches}."
      puts "Number of unchanged branches: #{self.numberOfUnchangedBranches}"
    end

    # Write table with changes to file.
    unless fileName == nil
      changes_file = File.new(fileName,"w")
      changes_string = "Root\t#{@branch[0].origin}\t#{@branch[0].originBeforeErrorIntroduction}\n"
      @branch.each do |b|
        if b.endCause == "speciation"
          changes_string << "#{b.id}\t#{b.termination}\t#{b.terminationBeforeErrorIntroduction}\n"
        end
      end
      changes_file.write(changes_string)
      if verbose
        puts "Wrote changes to file #{fileName}."
      end
    end

    if verbose
      endTime = Time.now
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end

  end

  def meanAgeError
    ageErrors = []
    @branch.each{|b| ageErrors << b.ageError}
    ageErrors.mean
  end

  def meanAge
    ages = []
    @branch.each{|b| ages << (b.origin+b.termination)/2.0}
    ages.mean
  end

  def numberOfOverestimatedBranches
    numberOfOverestimatedBranches = 0
    @branch.each do |b|
      numberOfOverestimatedBranches += 1 if b.origin+b.termination > b.originBeforeErrorIntroduction+b.terminationBeforeErrorIntroduction
    end
    numberOfOverestimatedBranches
  end

  def numberOfUnderestimatedBranches
    numberOfUnderestimatedBranches = 0
    @branch.each do |b|
      numberOfUnderestimatedBranches += 1 if b.origin+b.termination < b.originBeforeErrorIntroduction+b.terminationBeforeErrorIntroduction
    end
    numberOfUnderestimatedBranches
  end

  def numberOfUnchangedBranches
    numberOfUnchangedBranches = 0
    @branch.each do |b|
      numberOfUnchangedBranches += 1 if b.origin+b.termination == b.originBeforeErrorIntroduction+b.terminationBeforeErrorIntroduction
    end
    numberOfUnchangedBranches
  end

  def phylogeneticDiversity
    if @treeReconstructed
      @phylogeneticDiversity
    else
      raise "The phylogenetic diversity is only available after the tree has been reconstructed. This tree is not yet reconstructed. For complete trees, use the sum of species durations instead!"
    end
  end
  
  def gamma
    if @treeReconstructed
      @gamma
    else
      raise "The gamma statistic should only be calculated from a reconstructed tree. This tree is not yet reconstructed!"
    end
  end
  
  def assignBranchRates(mean = 0.01, standardDeviation = 0.005, autoCorrelation = false, nStepsPerTimeUnit = 1000, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "-------------------------------------------------phylsim.rb | Branch rates--------------------------------------------------"
      puts
    end
    @branchRatesMean = mean
    @branchRatesStandardDeviation = standardDeviation
    @branchRatesAutoCorrelation = autoCorrelation
    # Assign substitution rates to all branches:
    # - if autoCorrelation == false, the uncorrelated lognormal clock model of Drummond et al. (2006) is used.
    # - if autoCorrelation > 0, the CIR process model is used to infer node rates, as described in Lepage et al. (2006,2007).
    #   In this case, the value of sRAutoCorrelation will be used as the decorrelation time (1/b in Lepage et al. 2006 and 1/Theta in Lepage et al. 2007)
    #   The rate of the branch between nodes i and j is calculated according to the approximation given in Lepage et al. (2007): (R(i)+R(j))/2
    variance = standardDeviation**2

    raise "A positive value must be specified for parameter 'mean'!" if mean <= 0
    raise "A non-negative value must be specified for parameter 'standardDeviation'!" if variance < 0
    raise "The number of steps per time unit must be an positive integer!" unless nStepsPerTimeUnit.class == Fixnum and nStepsPerTimeUnit > 0

    if autoCorrelation == false or autoCorrelation == nil or autoCorrelation == 0 or variance == 0
      if verbose
        puts "Model:"
        puts "Uncorrelated lognormal (UCLN) distribution of branch rates (Drummond et al. 2006)"
        puts
        puts "Parameters:"
        puts "Mean: #{mean}"
        puts "Standard deviation: #{standardDeviation}"
        puts "Variance: #{(standardDeviation**2).round(4)}"
        puts
      end
      if variance == 0
        @branch.size.times do |b|
          @branch[b].addStartRate(mean)
          @branch[b].addEndRate(mean)
          @branch[b].addMeanRate(mean)
        end
      else
        # From the parameters 'mean' and 'variance' that refer to the lognormal distribution, calculate corresponding
        # parameters mu and sigma of the the non-logarithmized normal distribution.
        mu = Math.log(mean) - 0.5*Math.log(variance.to_f/mean**2 + 1)
        sigma = Math.sqrt(Math.log(variance.to_f/mean**2 + 1))
        normDist = Rubystats::NormalDistribution.new(mu,sigma)
        @branch.size.times do |b|
          rate = Math.exp(normDist.rng)
          @branch[b].addStartRate(rate)
          @branch[b].addEndRate(rate)
          @branch[b].addMeanRate(rate)
        end
      end
    elsif autoCorrelation.class == Float or autoCorrelation.class == Fixnum
      if autoCorrelation > 0
        if verbose
          puts "Model:"
          puts "Branch rate autocorrelation according to the Cox-Ingersoll-Ross (CIR) process (Lepage et al. 2006, 2007)"
          puts
          puts "Parameters:"
          puts "Stationary mean: #{mean}"
          puts "Standard deviation: #{standardDeviation}"
          puts "Variance: #{(standardDeviation**2).round(8)}"
          puts "Autocorrelation: #{autoCorrelation}"
          puts 
          print "Assigning substitution rates to all branches..."
        end
        # The following variables are named according to Lepage et al. (2006; p.224).
        a = mean.to_f                              # a is the stationary mean, and chosen to match mean.
        b = 1/autoCorrelation.to_f                 # b determines how fast the autocovariance goes to 0 with t.
        sigmaSquare = (b*variance)/(2*a)           # sigmaSquare is chosen so that the stationary variance 2*a*sigmaSquare/b matches the specified sRVariance (see Lepage et al. 2006; p. 224).
        warn "WARNING: The CIR process does not have a stationary distribution, because the standardDeviation is equal or greater than twice the mean" if 2*a*b <= sigmaSquare
        df = (4*a*b)/sigmaSquare                   # degrees of freedom for all pdfs of the non-central chi-square distributions for R(t) (see Lepage et al. 2006; p. 224).

        @r = RSRuby.instance
        parent = nil

        r0 = a                                     # The startRate of the root is chosen to match the specified mean of the stationary distribution sRMean (= a).
        @branch[0].addStartRate(r0)
        nSteps = (nStepsPerTimeUnit * @branch[0].duration).to_i + 1
        t = @branch[0].duration/nSteps.to_f
        constant = (2*b)/(sigmaSquare*(1-Math.exp(-b*t)))  # is ct on http://en.wikipedia.org/wiki/CIR_process
        stepRates = [r0]
        meanRate = 0
        nSteps.times do
          ncp = (4*b*stepRates.last*Math.exp(-b*t))/(sigmaSquare*(1-Math.exp(-b*t)))
          stepRates << (@r.rchisq(1, df, ncp))/(2*constant)
          meanRate += ((stepRates[-2]+stepRates[-1])/2.0)/nSteps.to_f
        end
        @branch[0].addEndRate(stepRates.last)
        @branch[0].addMeanRate(meanRate)
        startBranch = 1

        if @rootSplit == true
          @branch[1].addStartRate(r0)
          nSteps = (nStepsPerTimeUnit * @branch[1].duration).to_i + 1
          t = @branch[1].duration/nSteps.to_f
          constant = (2*b)/(sigmaSquare*(1-Math.exp(-b*t)))  # is ct on http://en.wikipedia.org/wiki/CIR_process
          stepRates = [r0]
          meanRate = 0
          nSteps.times do
            ncp = (4*b*stepRates.last*Math.exp(-b*t))/(sigmaSquare*(1-Math.exp(-b*t)))
            stepRates << (@r.rchisq(1, df, ncp))/(2*constant)
            meanRate += ((stepRates[-2]+stepRates[-1])/2.0)/nSteps.to_f
          end
          @branch[1].addEndRate(stepRates.last)
          @branch[1].addMeanRate(meanRate)
          startBranch = 2
        end

        startBranch.upto(@branch.size-1) do |b_| # here called 'b_' instead the usual 'b', as 'b' is preoccupied
          print "\rAssigning substitution rates to all branches... (#{b_}/#{@branch.size})" if verbose
          # Find parent.
          lower = 0
          upper = @branch.size
          parentFound = false
          until parentFound == true
            # Use the fact that @branch is sorted according to id for a quick parent search.
            p = lower+(upper-lower)/2
            if @branch[p].id[1..-1].to_i < @branch[b_].parentId[1..-1].to_i # Comparing the index numbers (e.g. '123') of branch ids (e.g. 'b123').
              lower = p
            elsif @branch[p].id[1..-1].to_i > @branch[b_].parentId[1..-1].to_i
              upper = p
            elsif @branch[p].id == @branch[b_].parentId
              parentFound = true
              parent = @branch[p]
              next
            else
              raise "Some problem occurred with the fast parent search!" # This should never be called at all.
            end
          end
          r0 = parent.endRate
          @branch[b_].addStartRate(r0)
          nSteps = (nStepsPerTimeUnit * @branch[b_].duration).to_i + 1
          t = @branch[b_].duration/nSteps.to_f
          constant = (2*b)/(sigmaSquare*(1-Math.exp(-b*t)))  # is ct on http://en.wikipedia.org/wiki/CIR_process
          stepRates = [r0]
          meanRate = 0
          nSteps.times do
            ncp = (4*b*stepRates.last*Math.exp(-b*t))/(sigmaSquare*(1-Math.exp(-b*t)))
            stepRates << (@r.rchisq(1, df, ncp))/(2*constant)
            meanRate += ((stepRates[-2]+stepRates[-1])/2.0)/nSteps.to_f
          end
          @branch[b_].addEndRate(stepRates.last)
          @branch[b_].addMeanRate(meanRate)
        end
      else
        raise "Parameter autoCorrelation must be either a positive value or the default setting (autoCorrelation = false)"
      end
    elsif autoCorrelation == true
      raise "Parameter autoCorrelation must be either a positive value or the default setting (autoCorrelation = false)"
    end

    # Set variable @branchRatesAssigned.
    @branchRatesAssigned = true

    # Report time consumed.
    if verbose
      endTime = Time.now
      puts "\rAssigned substitution rates to all branches.                                               "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
  end

  def modifyBranchRates(rateMultiplierMean = 1, rateMultiplierStandardVariation = 0, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "-------------------------------------------------phylsim.rb | Branch rates--------------------------------------------------"
      puts
    end
    raise "In order to modify branch rates, these must first be assigned!" unless @branchRatesAssigned
    raise "A positive number must be specified for parameter 'rateMultiplierMean'!" unless [Fixnum,Float].include?(rateMultiplierMean.class) and rateMultiplierMean > 0
    raise "A non-negative number must be specified for parameter 'rateMultiplierStandardVariation'!" unless [Fixnum,Float].include?(rateMultiplierStandardVariation.class) and rateMultiplierStandardVariation >= 0
    if rateMultiplierStandardVariation == 0
      @branch.each {|b| b.updateRate(b.rate * rateMultiplierMean)}
    else
      # From the parameters 'rateMultiplierMean' and 'rateMultiplierStandardVariation' that refer to the lognormal distribution, calculate corresponding
      # parameters mu and sigma of the the non-logarithmized normal distribution.
      rateMultiplierVariance = rateMultiplierStandardVariation**2.0
      mu = Math.log(rateMultiplierMean) - 0.5*Math.log(rateMultiplierVariance.to_f/rateMultiplierMean**2 + 1)
      sigma = Math.sqrt(Math.log(rateMultiplierVariance.to_f/rateMultiplierMean**2 + 1))
      normDist = Rubystats::NormalDistribution.new(mu,sigma)
      @branch.each {|b| b.updateRate(b.rate * Math.exp(normDist.rng))}
    end
    if verbose
      endTime = Time.now
      puts "\rMultiplied substitution rates of all branches with #{rateMultiplierMean.round(3)}+/-#{rateMultiplierStandardVariation.round(3)}. Mean branch rate is now: #{self.meanBranchRate(verbose = false).round(8)}"
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
  end

  def meanBranchRate(verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "-----------------------------------------------------phylsim.rb | Info------------------------------------------------------"
      puts 
    end

    branchRateSum = 0
    @branch.each {|b| branchRateSum += b.rate}
    meanBranchRate = branchRateSum.to_f/@branch.size.to_f

    if verbose
      endTime = Time.now
      puts "\rMean branch rate calculated: #{meanBranchRate}."
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
    return meanBranchRate

  end

  def evolveTraitsBM(initial = 0, rate = 1, elevatedRate = nil, elevatedPeriod = nil, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "-----------------------------------------------phylsim.rb | Trait evolution-------------------------------------------------"
      puts
    end

    if @traitsEvolved
      if verbose
        puts "WARNING: Traits have previously been added to this tree. These will be replaced."
        puts
      end
      # Delete all previously generated traits
      branch.each do |b|
        b.addStartTrait(nil)
        b.addEndTrait(nil)
      end
    end

    # Check the arguments:
    # Check the rate.
    if rate <= 0
      raise "Positive values must be specified for parameter 'rate'!"
    end

    # Check the elevated rate.
    unless elevatedRate == nil
      if elevatedRate <= 0
        raise "Positive values must be specified for parameter 'elevatedRate'!"
      end
    end

    # Check the elevated period.
    unless elevatedRate == nil
      if elevatedPeriod == nil
        raise "If an elevated rate is given, a period of elevated diversification must also be specified!"
      else
        if elevatedPeriod.class == Array and elevatedPeriod.size == 2
          elevatedStart = elevatedPeriod.max
          elevatedEnd = elevatedPeriod.min
          raise "Both the start and end of the elevated diversification period must be non-negative numbers!" if elevatedEnd < 0
        else
          raise "An array with two non-negative numbers must be specified for the period of elevated diversification!"
        end
      end
    end

    if verbose
      puts "Model:"
      puts "Brownian motion trait evolution"
      puts
      puts "Brownian motion initial = #{initial}"
      puts "Brownian motion rate = #{rate}"
      unless elevatedRate == nil
        puts "Elevated period of diversification = #{elevatedStart}-#{elevatedEnd}"
        puts "Elevated Brownian motion rate = #{elevatedRate}"
      end
      puts
    end

    print "Generating an ancestral trait state..." if verbose
    # Use the value given for parameter 'initial' as the ancestral state and assign it to the first branch if @rootSplit == false, or to the first two branches if @rootSplit == true.
    if @rootSplit == false
      @branch[0].addStartTrait(trait = initial)
      branchesNow = [@branch[0]]
    else
      @branch[0].addStartTrait(trait = initial)
      @branch[1].addStartTrait(trait = initial)
      branchesNow = [@branch[0],@branch[1]]
    end

    print "\rEvolving the ancestral trait along the tree..." if verbose
    # Evolve the start trait throughout the tree, using the specified rate (and elevated rate).
    change = true
    count = 0
    until change == false
      change = false
      branchesNext = []
      branchesNow.size.times do |b|
        print "\rEvolving the ancestral trait along the tree... (#{count}/#{@branch.size} branches)" if verbose
        branchStart = branchesNow[b].origin
        branchStartTrait = branchesNow[b].startTrait
        branchEnd = branchesNow[b].termination
        if elevatedRate == nil or branchEnd >= elevatedStart or branchStart <= elevatedEnd
          sigma = Math.sqrt((branchStart-branchEnd)*rate)
          normDist = Rubystats::NormalDistribution.new(branchStartTrait,sigma)
          branchEndTrait = normDist.rng
        elsif branchStart > elevatedStart and branchEnd < elevatedStart and branchEnd >= elevatedEnd
          sigma = Math.sqrt((branchStart-elevatedStart)*rate)
          normDist = Rubystats::NormalDistribution.new(branchStartTrait,sigma)
          branchInBetweenTrait = normDist.rng
          sigma = Math.sqrt((elevatedStart-branchEnd)*elevatedRate)
          normDist = Rubystats::NormalDistribution.new(branchInBetweenTrait,sigma)
          branchEndTrait = normDist.rng
        elsif branchStart <= elevatedStart and branchStart > elevatedEnd and branchEnd < elevatedEnd
          sigma = Math.sqrt((branchStart-elevatedEnd)*elevatedRate)
          normDist = Rubystats::NormalDistribution.new(branchStartTrait,sigma)
          branchInBetweenTrait = normDist.rng
          sigma = Math.sqrt((elevatedEnd-branchEnd)*rate)
          normDist = Rubystats::NormalDistribution.new(branchInBetweenTrait,sigma)
          branchEndTrait = normDist.rng
        elsif branchStart > elevatedStart and branchEnd < elevatedEnd
          sigma = Math.sqrt((branchStart-elevatedStart)*rate)
          normDist = Rubystats::NormalDistribution.new(branchStartTrait,sigma)
          branchInBetweenTrait1 = normDist.rng
          sigma = Math.sqrt((elevatedStart-elevatedEnd)*elevatedRate)
          normDist = Rubystats::NormalDistribution.new(branchInBetweenTrait1,sigma)
          branchInBetweenTrait2 = normDist.rng
          sigma = Math.sqrt((elevatedEnd-branchEnd)*rate)
          normDist = Rubystats::NormalDistribution.new(branchInBetweenTrait2,sigma)
          branchEndTrait = normDist.rng
        elsif branchStart <= elevatedStart and branchEnd >= elevatedEnd
          sigma = Math.sqrt((branchStart-branchEnd)*elevatedRate)
          normDist = Rubystats::NormalDistribution.new(branchStartTrait,sigma)
          branchEndTrait = normDist.rng
        else
          raise "This should never happen!"
        end
        raise "branchEndTrait is nil!" if branchEndTrait == nil
        branchesNow[b].addEndTrait(branchEndTrait)
        count += 1

        if branchesNow[b].daughterId[0][1..-1].to_i > 0 # This means only if daughter ids are not 'none' or 'unborn'.
          # Do another fast id search to find the two daughters.
          lower = -1
          upper = @branch.size
          d0Found = false
          until d0Found
            bb = lower+(upper-lower)/2
            if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[0][1..-1].to_i
              lower = bb
            elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[0][1..-1].to_i
              upper = bb
            elsif @branch[bb].id == branchesNow[b].daughterId[0]
              d0Found = true
              branchesNext << @branch[bb]
              @branch[bb].addStartTrait(branchesNow[b].endTrait)
              next
            else
              raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
            end
          end
        end
        if branchesNow[b].daughterId[1][1..-1].to_i > 0 # This means only if they daughter ids are not 'none' or 'unborn'.
          # Do another fast id search to find the two daughters.
          lower = -1
          upper = @branch.size
          d1Found = false
          until d1Found
            bb = lower+(upper-lower)/2
            if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[1][1..-1].to_i
              lower = bb
            elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[1][1..-1].to_i
              upper = bb
            elsif @branch[bb].id == branchesNow[b].daughterId[1]
              d1Found = true
              branchesNext << @branch[bb]
              @branch[bb].addStartTrait(branchesNow[b].endTrait)
              next
            else
              raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
            end
          end
        end
      end
      branchesNow = branchesNext
      change = true if branchesNext.size > 0
      
    end # until change == false

    @traitsEvolved = true
    if verbose
      endTime = Time.now
      puts "\rTraits added to all branches.                                                               "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
  end

  def evolveTraitsOU(initial = 0, optimum = 0, limits = nil, rate = 1, selectionStrength = 1, elevatedRate = nil, elevatedPeriod = nil, nStepsPerTimeUnit = 1000, verbose = true)

    # The optimum is parameter Theta, and selection strength is parameter alpha in Butler & King (2004).
    # All other parameters are identical to those in evolveTraitsBM.

    if verbose
      startTime = Time.now
      puts
      puts "-----------------------------------------------phylsim.rb | Trait evolution-------------------------------------------------"
      puts
    end

    if @traitsEvolved
      if verbose
        puts "WARNING: Traits have previously been added to this tree. These will be replaced."
        puts
      end
      # Delete all previously generated traits
      branch.each do |b|
        b.addStartTrait(nil)
        b.addEndTrait(nil)
      end
    end

    # Check the arguments:
    # Check the optimum.
    if optimum.class == Array
      optimum.each do |o|
        raise "All optima must be specified as integers or floats!" unless o.class == Fixnum or o.class == Float
      end
    else
      raise "The optimum must be specified as integer or float!" unless optimum.class == Fixnum or optimum.class == Float
    end

    # Check the rate.
    if rate <= 0
      raise "Positive values must be specified for parameter 'rate'!"
    end

    # Check the limit.
    if limits == nil
      lower_limit = -Float::INFINITY
      upper_limit = Float::INFINITY      
    else
      if limits.class == Array and limits.size == 2
        if limits[0] == nil or limits[0] == 'inf'
          if limits[1].class == Fixnum or limits[1].class == Float
            sole_limit = limits[1]
          else
            raise "At least one of the two values given for parameter 'limits' must be a number!"
          end
        elsif limits[1] == nil or limits[1] == 'inf'
          if limits[0].class == Fixnum or limits[0].class == Float
            sole_limit = limits[0]
          else
            raise "At least one of the two values given for parameter 'limits' must be a number!"
          end
        end
        if sole_limit == nil
          lower_limit = limits.min.to_f
          upper_limit = limits.max.to_f
        else
          if sole_limit > optimum
            lower_limit = -Float::INFINITY
            upper_limit = sole_limit.to_f
          else
            lower_limit = sole_limit.to_f
            upper_limit = Float::INFINITY
          end
        end
        raise "The minimum and maximum values given for parameter 'limits' must be different!" if lower_limit == upper_limit
        if optimum.class == Fixnum or optimum.class == Float
          raise "The optimum lies outside the limits!" if lower_limit > optimum or upper_limit < optimum
        elsif optimum.class == Array
          raise "At least one of the optima lies outside the limits!" if lower_limit > optimum.min or upper_limit < optimum.max
        end
      else
        raise "Either nil, or an array of size two should be given for parameter 'limits'!"
      end
    end

    # Check the selection strength.
    if selectionStrength < 0
      raise "Positive values must be specified for parameter 'selectionStrength'!"
    end

    # Check the elevated rate.
    unless elevatedRate == nil
      if elevatedRate <= 0
        raise "Positive values must be specified for parameter 'elevatedRate'!"
      end
    end

    # Check the elevated period.
    unless elevatedRate == nil
      if elevatedPeriod == nil
        raise "If an elevated rate is given, a period of elevated diversification must also be specified!"
      else
        if elevatedPeriod.class == Array and elevatedPeriod.size == 2
          elevatedStart = elevatedPeriod.max
          elevatedEnd = elevatedPeriod.min
          raise "Both the start and end of the elevated diversification period must be non-negative numbers!" if elevatedEnd < 0
        else
          raise "An array with two non-negative numbers must be specified for the period of elevated diversification!"
        end
      end
    end

    # Check the number of steps per time unit.
    if nStepsPerTimeUnit <= 0
      raise "Positive values must be specified for parameter 'nStepsPerTimeUnit'!"
    end

    if verbose
      puts "Model:"
      puts "Ornstein-Uhlenbeck trait evolution"
      puts
      puts "Ornstein-Uhlenbeck initial = #{initial}"
      if optimum.class == Array
        print "Ornstein-Uhlenbeck optima = "
        optimum[0..-2].each {|o| print "#{o}, "}
        puts "#{optimum[-1]}"
      else
        puts "Ornstein-Uhlenbeck optimum = #{optimum}"
      end
      unless limits == nil
        puts "Limits = #{lower_limit}-#{upper_limit}"
      end
      puts "Ornstein-Uhlenbeck selection strength = #{selectionStrength}"
      puts "Ornstein-Uhlenbeck rate = #{rate}"
      unless elevatedRate == nil
        puts "Elevated period of diversification = #{elevatedStart}-#{elevatedEnd}"
        puts "Elevated Ornstein-Uhlenbeck rate = #{elevatedRate}"
      end
      puts
    end

    print "Generating an ancestral trait state..." if verbose
    # Use the value given for parameter 'initial' as the ancestral state and assign it to the first branch if @rootSplit == false, or to the first two branches if @rootSplit == true.
    if @rootSplit == false
      @branch[0].addStartTrait(trait = initial)
      branchesNow = [@branch[0]]
    else
      @branch[0].addStartTrait(trait = initial)
      @branch[1].addStartTrait(trait = initial)
      branchesNow = [@branch[0],@branch[1]]
    end

    print "\rEvolving the ancestral trait along the tree..." if verbose
    # Evolve the start trait throughout the tree, using the specified rate (and elevated rate).
    change = true
    count = 0
    until change == false
      change = false
      branchesNext = []
      branchesNow.size.times do |b|
        print "\rEvolving the ancestral trait along the tree... (#{count}/#{@branch.size} branches)" if verbose
        branchStart = branchesNow[b].origin
        branchEnd = branchesNow[b].termination
        branchSectionStart = branchStart
        branchSectionEnd = branchSectionStart - 1/nStepsPerTimeUnit.to_f
        branchSectionStartTrait = branchesNow[b].startTrait
        while branchSectionEnd > branchEnd
          # Find nearest optimum.
          if optimum.class == Array
            optimum_distances = []
            optimum.each {|o| optimum_distances << (o - branchSectionStartTrait).abs}
            nearest_optimum = optimum[optimum_distances.index(optimum_distances.min)]
          else
            nearest_optimum = optimum
          end
          if elevatedRate == nil or branchSectionEnd >= elevatedStart or branchSectionStart <= elevatedEnd
            sigma = Math.sqrt((branchSectionStart-branchSectionEnd)*rate)
            normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionEndTrait = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-branchSectionEnd) + normDist.rng
              valid_trait_value_found = true if branchSectionEndTrait > lower_limit and branchSectionEndTrait < upper_limit
            end
          elsif branchSectionStart > elevatedStart and branchSectionEnd < elevatedStart and branchSectionEnd >= elevatedEnd
            sigma = Math.sqrt((branchSectionStart-elevatedStart)*rate)
            normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionInBetweenTrait = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-elevatedStart) + normDist.rng
              valid_trait_value_found = true if branchSectionInBetweenTrait > lower_limit and branchSectionInBetweenTrait < upper_limit
            end
            sigma = Math.sqrt((elevatedStart-branchSectionEnd)*elevatedRate)
            normDist = Rubystats::NormalDistribution.new(branchSectionInBetweenTrait,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionEndTrait = selectionStrength*(nearest_optimum-branchSectionInBetweenTrait)*(elevatedStart-branchSectionEnd) + normDist.rng
              valid_trait_value_found = true if branchSectionEndTrait > lower_limit and branchSectionEndTrait < upper_limit
            end
          elsif branchSectionStart <= elevatedStart and branchSectionStart > elevatedEnd and branchSectionEnd < elevatedEnd
            sigma = Math.sqrt((branchSectionStart-elevatedEnd)*elevatedRate)
            normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionInBetweenTrait = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-elevatedEnd) + normDist.rng
              valid_trait_value_found = true if branchSectionInBetweenTrait > lower_limit and branchSectionInBetweenTrait < upper_limit
            end
            sigma = Math.sqrt((elevatedEnd-branchSectionEnd)*rate)
            normDist = Rubystats::NormalDistribution.new(branchSectionInBetweenTrait,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionEndTrait = selectionStrength*(nearest_optimum-branchSectionInBetweenTrait)*(elevatedEnd-branchSectionEnd) + normDist.rng
              valid_trait_value_found = true if branchSectionEndTrait > lower_limit and branchSectionEndTrait < upper_limit
            end
          elsif branchSectionStart > elevatedStart and branchSectionEnd < elevatedEnd
            sigma = Math.sqrt((branchSectionStart-elevatedStart)*rate)
            normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionInBetweenTrait1 = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-elevatedStart) + normDist.rng
              valid_trait_value_found = true if branchSectionInBetweenTrait1 > lower_limit and branchSectionInBetweenTrait1 < upper_limit
            end
            sigma = Math.sqrt((elevatedStart-elevatedEnd)*elevatedRate)
            normDist = Rubystats::NormalDistribution.new(branchSectionInBetweenTrait1,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionInBetweenTrait2 = selectionStrength*(nearest_optimum-branchSectionInBetweenTrait1)*(elevatedStart-elevatedEnd) + normDist.rng
              valid_trait_value_found = true if branchSectionInBetweenTrait2 > lower_limit and branchSectionInBetweenTrait2 < upper_limit
            end
            sigma = Math.sqrt((elevatedEnd-branchSectionEnd)*rate)
            normDist = Rubystats::NormalDistribution.new(branchSectionInBetweenTrait2,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionEndTrait = selectionStrength*(nearest_optimum-branchSectionInBetweenTrait2)*(elevatedEnd-branchSectionEnd) + normDist.rng
              valid_trait_value_found = true if branchSectionEndTrait > lower_limit and branchSectionEndTrait < upper_limit
            end
          elsif branchSectionStart <= elevatedStart and branchSectionEnd >= elevatedEnd
            sigma = Math.sqrt((branchSectionStart-branchSectionEnd)*elevatedRate)
            normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
            valid_trait_value_found = false
            until valid_trait_value_found == true
              branchSectionEndTrait = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-branchSectionEnd) + normDist.rng
              valid_trait_value_found = true if branchSectionEndTrait > lower_limit and branchSectionEndTrait < upper_limit
            end
          else
            raise "This should never happen!"
          end
          raise "branchSectionEndTrait is nil!" if branchSectionEndTrait == nil
          branchSectionStart = branchSectionEnd
          branchSectionEnd = branchSectionStart - 1/nStepsPerTimeUnit.to_f
          branchSectionStartTrait = branchSectionEndTrait
        end
        # Find nearest optimum.
        if optimum.class == Array
          optimum_distances = []
          optimum.each {|o| optimum_distances << (o - branchSectionStartTrait).abs}
          nearest_optimum = optimum[optimum_distances.index(optimum_distances.min)]
        else
          nearest_optimum = optimum
        end
        if elevatedRate == nil or branchEnd >= elevatedStart or branchSectionStart <= elevatedEnd
          sigma = Math.sqrt((branchSectionStart-branchEnd)*rate)
          normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchEndTrait = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-branchEnd) + normDist.rng
            valid_trait_value_found = true if branchEndTrait > lower_limit and branchEndTrait < upper_limit
          end
        elsif branchSectionStart > elevatedStart and branchEnd < elevatedStart and branchEnd >= elevatedEnd
          sigma = Math.sqrt((branchSectionStart-elevatedStart)*rate)
          normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchSectionInBetweenTrait = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-elevatedStart) + normDist.rng
            valid_trait_value_found = true if branchSectionInBetweenTrait > lower_limit and branchSectionInBetweenTrait < upper_limit
          end
          sigma = Math.sqrt((elevatedStart-branchEnd)*elevatedRate)
          normDist = Rubystats::NormalDistribution.new(branchSectionInBetweenTrait,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchEndTrait = selectionStrength*(nearest_optimum-branchSectionInBetweenTrait)*(elevatedStart-branchEnd) + normDist.rng
            valid_trait_value_found = true if branchEndTrait > lower_limit and branchEndTrait < upper_limit
          end
        elsif branchSectionStart <= elevatedStart and branchSectionStart > elevatedEnd and branchEnd < elevatedEnd
          sigma = Math.sqrt((branchSectionStart-elevatedEnd)*elevatedRate)
          normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchSectionInBetweenTrait = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-elevatedEnd) + normDist.rng
            valid_trait_value_found = true if branchSectionInBetweenTrait > lower_limit and branchSectionInBetweenTrait < upper_limit
          end
          sigma = Math.sqrt((elevatedEnd-branchEnd)*rate)
          normDist = Rubystats::NormalDistribution.new(branchSectionInBetweenTrait,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchEndTrait = selectionStrength*(nearest_optimum-branchSectionInBetweenTrait)*(elevatedEnd-branchEnd) + normDist.rng
            valid_trait_value_found = true if branchEndTrait > lower_limit and branchEndTrait < upper_limit
          end
        elsif branchSectionStart > elevatedStart and branchEnd < elevatedEnd
          sigma = Math.sqrt((branchSectionStart-elevatedStart)*rate)
          normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchSectionInBetweenTrait1 = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-elevatedStart) + normDist.rng
            valid_trait_value_found = true if branchSectionInBetweenTrait1 > lower_limit and branchSectionInBetweenTrait1 < upper_limit
          end
          sigma = Math.sqrt((elevatedStart-elevatedEnd)*elevatedRate)
          normDist = Rubystats::NormalDistribution.new(branchSectionInBetweenTrait1,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchSectionInBetweenTrait2 = selectionStrength*(nearest_optimum-branchSectionInBetweenTrait1)*(elevatedStart-elevatedEnd) + normDist.rng
            valid_trait_value_found = true if branchSectionInBetweenTrait2 > lower_limit and branchSectionInBetweenTrait2 < upper_limit
          end
          sigma = Math.sqrt((elevatedEnd-branchEnd)*rate)
          normDist = Rubystats::NormalDistribution.new(branchSectionInBetweenTrait2,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchEndTrait = selectionStrength*(nearest_optimum-branchSectionInBetweenTrait2)*(elevatedEnd-branchEnd) + normDist.rng
            valid_trait_value_found = true if branchEndTrait > lower_limit and branchEndTrait < upper_limit
          end
        elsif branchSectionStart <= elevatedStart and branchEnd >= elevatedEnd
          sigma = Math.sqrt((branchSectionStart-branchEnd)*elevatedRate)
          normDist = Rubystats::NormalDistribution.new(branchSectionStartTrait,sigma)
          valid_trait_value_found = false
          until valid_trait_value_found == true
            branchEndTrait = selectionStrength*(nearest_optimum-branchSectionStartTrait)*(branchSectionStart-branchEnd) + normDist.rng
            valid_trait_value_found = true if branchEndTrait > lower_limit and branchEndTrait < upper_limit
          end
        else
          raise "This should never happen!"
        end
        raise "branchEndTrait is nil!" if branchEndTrait == nil
        branchesNow[b].addEndTrait(branchEndTrait)
        count += 1

        if branchesNow[b].daughterId[0][1..-1].to_i > 0 # This means only if daughter ids are not 'none' or 'unborn'.
          # Do another fast id search to find the two daughters.
          lower = -1
          upper = @branch.size
          d0Found = false
          until d0Found
            bb = lower+(upper-lower)/2
            if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[0][1..-1].to_i
              lower = bb
            elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[0][1..-1].to_i
              upper = bb
            elsif @branch[bb].id == branchesNow[b].daughterId[0]
              d0Found = true
              branchesNext << @branch[bb]
              @branch[bb].addStartTrait(branchesNow[b].endTrait)
              next
            else
              raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
            end
          end
        end
        if branchesNow[b].daughterId[1][1..-1].to_i > 0 # This means only if they daughter ids are not 'none' or 'unborn'.
          # Do another fast id search to find the two daughters.
          lower = -1
          upper = @branch.size
          d1Found = false
          until d1Found
            bb = lower+(upper-lower)/2
            if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[1][1..-1].to_i
              lower = bb
            elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[1][1..-1].to_i
              upper = bb
            elsif @branch[bb].id == branchesNow[b].daughterId[1]
              d1Found = true
              branchesNext << @branch[bb]
              @branch[bb].addStartTrait(branchesNow[b].endTrait)
              next
            else
              raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
            end
          end
        end
      end
      branchesNow = branchesNext
      change = true if branchesNext.size > 0
      
    end # until change == false

    @traitsEvolved = true
    if verbose
      endTime = Time.now
      puts "\rTraits added to all branches.                                                               "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
  end

  def evolveSequences(sequenceLength, substitutionModel = "ecm", threads = 1, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "----------------------------------------------phylsim.rb | Sequence evolution-----------------------------------------------"
      puts
    end

    if @sequencesEvolved
      if verbose
        puts "WARNING: Nucleotide sequences have previously been added to this tree. These will be replaced."
        puts
      end
      # Delete all previously generated sequences
      branch.each do |b|
        b.addStartSeq(nil)
        b.addEndSeq(nil)
      end
    end

    # Make sure substitution rates have been assigned previously.
    raise "In order to add nucleotide sequences, branch rates must first be assigned!" unless @branchRatesAssigned

    # Check the arguments:
    # Check the specified sequence length
    if sequenceLength <= 0
      raise "Positive values must be specified for parameter 'sequenceLength'!"
    end

    # Check the substitution model (currently, only the unrestricted empirircal codon model of Kosiol et al. 2007 and the GTR model are supported)
    @ecmFileName = "#{$libPath}resources/ECMunrest.dat"
    if substitutionModel =~ /^[Ee][Cc][Mm]$/
      @substitutionModel = ECM.new(@ecmFileName)
      if verbose
        puts "Model:"
        puts "Empirical codon model (ECM) for protein sequence evolution (Kosiol et al. 2007)"
        puts
      end
    elsif substitutionModel =~ /^[Gg][Tt][Rr]([Gg][Aa][Mm][Mm][Aa])?\((.+)\)$/
      string = $2
      @substitutionModel = GTR.new(string)
      if verbose
        puts "Model:"
        puts "General time reversible (GTR) model of sequence evolution (Tavare 1986)"
        puts
      end
    else
      raise "The specified substitution model is unknown. It should be either 'ecm(string)' or 'gtr(string)', whereby 'string' should either specify the name of a file containing all parameter values, or, if gtr is chosen, a series of parameter values in the order freq(a), freq(c), freq(g), rate(a<->c), rate(a<->g), rate(a<->t), rate(c<->g), rate(c<->t), e.g. 'gtr(0.241373,0.184797,0.256406,0.890985,2.499690,0.755603,0.899543,2.512762)'!"
    end

    if verbose
      puts "Parameters:"
      if @substitutionModel.class == ECM
        puts "61 codon frequencies and 1830 substitution rates as specified in '#{@ecmFileName}'"
      elsif @substitutionModel.class == GTR
        puts "(drawn from normal distributions that approximate the model-averaged estimates of Arbiza et al. 2011)" if ["mammal","mammals","mammalia","vertebrate","vertebrates","vertebrata"].include?(string.downcase)
        puts "Frequency of adenine: pi(A) = #{@substitutionModel.pi[0].round(5)}"
        puts "Frequency of cytosine: pi(C) = #{@substitutionModel.pi[1].round(5)}"
        puts "Frequency of guanine: pi(G) = #{@substitutionModel.pi[2].round(5)}"
        puts "Frequency of thymine: pi(T) = #{@substitutionModel.pi[3].round(5)}"
        puts "Relative substitution rate adenine to cytosine: rate(A<->C) = #{@substitutionModel.s[0].round(5)}"
        puts "Relative substitution rate adenine to guanine: rate(A<->G) = #{@substitutionModel.s[1].round(5)}"
        puts "Relative substitution rate adenine to thymine: rate(A<->T) = #{@substitutionModel.s[2].round(5)}"
        puts "Relative substitution rate cytosine to guanine: rate(C<->G) = #{@substitutionModel.s[3].round(5)}"
        puts "Relative substitution rate cytosine to thymine: rate(C<->T) = #{@substitutionModel.s[4].round(5)}"
        puts "Relative substitution rate guanine to thymine: rate(G<->T) = #{@substitutionModel.s[5].round(5)}"
        puts "Shape parameter of gamma distribution for among-site rate variation: alpha = #{@substitutionModel.alpha.round(5)}" unless @substitutionModel.alpha == nil
      else
        raise "The class of the substitution model is neither 'ECM' nor 'GTR'!"
      end
      puts
    end

    # Check the number of threads
    if threads.class == Fixnum and threads > 1
      @threads = threads
    else
      @threads = 1
    end

    print "Generating an ancestral nucleotide sequence of #{sequenceLength} bp..." if verbose
    # Randomly generate a root sequence of a given length, and assign it to the first branch if @rootSplit == false, or to the first two branches if @rootSplit == true.
    rootSeq = @substitutionModel.randomSeq(seqLength = sequenceLength)
    if @rootSplit == false
      @branch[0].addStartSeq(seq = rootSeq)
      branchesNow = [@branch[0]]
    else
      @branch[0].addStartSeq(seq = rootSeq)
      @branch[1].addStartSeq(seq = rootSeq)
      branchesNow = [@branch[0],@branch[1]]
    end
    
    print "\rEvolving the ancestral nucleotide sequence along the tree..." if verbose
    # Evolve the start sequence throughout the tree, according to the assigned branch rates.
    if @threads == 1
      change = true
      count = 0
      until change == false
        change = false
        branchesNext = []
        branchesNow.size.times do |b|
          print "\rEvolving the ancestral nucleotide sequence along the tree... (#{count}/#{@branch.size} branches)" if verbose
          returnArray = @substitutionModel.evolveSeq(startSeq = branchesNow[b].startSeq,expectedNumberOfSubstitutionsPerSite =  branchesNow[b].expectedNumberOfSubstitutionsPerSite)
          branchesNow[b].addEndSeq(returnArray[0])
          branchesNow[b].addActualNumberOfSubstitutions(returnArray[1])
          branchesNow[b].addActualNumberOfSubstitutionsPerSite(returnArray[1]/sequenceLength.to_f)
          count += 1

          if branchesNow[b].daughterId[0][1..-1].to_i > 0 # This means only if they daughter ids are not 'none' or 'unborn'.
            # Do another fast id search to find the two daughters.
            lower = -1
            upper = @branch.size
            d0Found = false
            until d0Found
              bb = lower+(upper-lower)/2
              if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[0][1..-1].to_i
                lower = bb
              elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[0][1..-1].to_i
                upper = bb
              elsif @branch[bb].id == branchesNow[b].daughterId[0]
                d0Found = true
                branchesNext << @branch[bb]
                @branch[bb].addStartSeq(branchesNow[b].endSeq)
                next
              else
                raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
              end
            end
          end
          if branchesNow[b].daughterId[1][1..-1].to_i > 0 # This means only if they daughter ids are not 'none' or 'unborn'.
            # Do another fast id search to find the two daughters.
            lower = -1
            upper = @branch.size
            d1Found = false
            until d1Found
              bb = lower+(upper-lower)/2
              if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[1][1..-1].to_i
                lower = bb
              elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[1][1..-1].to_i
                upper = bb
              elsif @branch[bb].id == branchesNow[b].daughterId[1]
                d1Found = true
                branchesNext << @branch[bb]
                @branch[bb].addStartSeq(branchesNow[b].endSeq)
                next
              else
                raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
              end
            end
          end
        end
        branchesNow = branchesNext
        change = true if branchesNext.size > 0
        
      end # until change == false

    elsif @threads.class == Fixnum and @threads > 1

      # Memorize properties of the array @branch, so we can later make sure that the array is still the same.
      branchSizeBefore = @branch.size
      idsBefore = []
      @branch.each {|b| idsBefore << b.id}

      # Prepare the single-threaded sequence evolution until branchesNext.size == @threads
      change = true
      branchesPrevious = []
      count = 0
      branchesNext = []
      until change == false or branchesNext.size == @threads do
        change = false
        branchesNow.size.times do |b|
          print "\rEvolving the ancestral nucleotide sequence along the tree... (#{count}/#{@branch.size})" if verbose
          returnArray = @substitutionModel.evolveSeq(startSeq = branchesNow[b].startSeq,expectedNumberOfSubstitutionsPerSite =  branchesNow[b].expectedNumberOfSubstitutionsPerSite)
          branchesNow[b].addEndSeq(returnArray[0])
          branchesNow[b].addActualNumberOfSubstitutions(returnArray[1])
          branchesNow[b].addActualNumberOfSubstitutionsPerSite(returnArray[1]/sequenceLength.to_f)
          count += 1

          if branchesNow[b].daughterId[0][1..-1].to_i > 0 # This means only if they daughter ids are not 'none' or 'unborn'.
            # Do another fast id search to find the two daughters.
            lower = -1
            upper = @branch.size
            d0Found = false
            until d0Found
              bb = lower+(upper-lower)/2
              if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[0][1..-1].to_i
                lower = bb
              elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[0][1..-1].to_i
                upper = bb
              elsif @branch[bb].id == branchesNow[b].daughterId[0]
                d0Found = true
                branchesNext << @branch[bb]
                @branch[bb].addStartSeq(branchesNow[b].endSeq)
                next
              else
                raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
              end
            end
          end
          if branchesNow[b].daughterId[1][1..-1].to_i > 0 # This means only if they daughter ids are not 'none' or 'unborn'.
            # Do another fast id search to find the two daughters.
            lower = -1
            upper = @branch.size
            d1Found = false
            until d1Found
              bb = lower+(upper-lower)/2
              if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[1][1..-1].to_i
                lower = bb
              elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[1][1..-1].to_i
                upper = bb
              elsif @branch[bb].id == branchesNow[b].daughterId[1]
                d1Found = true
                branchesNext << @branch[bb]
                @branch[bb].addStartSeq(branchesNow[b].endSeq)
                next
              else
                raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
              end
            end
          end
        end
        if count <= 2
          branchesNow.each {|bn| branchesPrevious << bn}
        else
          branchesPrevious << branchesNow[0] # When count > 2, there can only be one item in branchesNow anyway.
        end

        if branchesNext.size > 0
          change = true
          if branchesNext.size < @threads
            # Find the branch among branchesNext that has the largest progeny.
            maxProgenySize = 0
            branchesNext.each {|bn| maxProgenySize = bn.extantProgenyId.size if bn.extantProgenyId.size > maxProgenySize}
            # Now take the first branch among branchesNext that has the maxProgenySize, remove it from branchesNext and put it into branchesNow.
            branchesNext.size.times do |bn|
              if branchesNext[bn].extantProgenyId.size == maxProgenySize
                branchesNow = [branchesNext[bn]]
                branchesNext[bn] = nil
                branchesNext.compact!
                break
              end
            end
          end
        end

      end # until change == false or branchesNext.size == @threads

      # This point is reached if either change == false (unlikely, only if the number of threads is large, but the number of branches is very small) or branchesNext.size == @threads.
      # If it's for the latter, and change is still true, then we must now split all branches into individual bunches for each thread.
      if change == true

        # Make sure @threads is the same as branchesNext.size
        raise "For some reason, the size of array 'branchesNext' is not the same as the number of threads!" unless branchesNext.size == @threads

        # First of all, calculate process priorities. To do so, we need to know the size of the sampled extant progenies of each branch in array 'branchesNext'.
        threadPriorities = []
        maxProgenySize = branchesNext[0].extantProgenyId.size
        minProgenySize = branchesNext[0].extantProgenyId.size
        branchesNext.each do |bn|
          maxProgenySize = bn.extantProgenyId.size if bn.extantProgenyId.size > maxProgenySize
          minProgenySize = bn.extantProgenyId.size if bn.extantProgenyId.size < minProgenySize
        end

        branchesNext.each do |bn|
          if minProgenySize == maxProgenySize
            threadPriorities << 20
          else
            threadPriorities << ((((bn.extantProgenyId.size-minProgenySize)/(maxProgenySize-minProgenySize).to_f)**2)*(-20) + 20).round
          end
        end

        # Define file names to which parts of the @branch array will be written.
        branchPackageDumpFileNames = []
        @threads.times do |t|
          branchPackageDumpFileNames << ".branchPackage#{(t+1).to_s.rjust((@threads+1).to_s.length).gsub(" ","0")}.dmp"
        end

        # Define file names to which the threadCount will be written.
        threadCountFileNames = []
        @threads.times do |t|
          threadCountFileNames << ".threadCount#{(t+1).to_s.rjust((@threads+1).to_s.length).gsub(" ","0")}.txt"
        end

        @threads.times do |t|
          Process.fork do
            Process.setpriority(Process::PRIO_PROCESS, 0, threadPriorities[t])

            # Prepare.
            branchesNow = [branchesNext[t]]
            change = true
            threadCount = 0

            # Create the threadCountFile.
            threadCountFile = File.new(threadCountFileNames[t],"w")
            threadCountFile.write("0,r")
            threadCountFile.close

            until change == false
              change = false
              branchesNext = []
              branchesNow.size.times do |b|
                returnArray = @substitutionModel.evolveSeq(startSeq = branchesNow[b].startSeq,expectedNumberOfSubstitutionsPerSite = branchesNow[b].expectedNumberOfSubstitutionsPerSite)
                branchesNow[b].addEndSeq(returnArray[0])
                branchesNow[b].addActualNumberOfSubstitutions(returnArray[1])
                branchesNow[b].addActualNumberOfSubstitutionsPerSite(returnArray[1]/sequenceLength.to_f)
                threadCount += 1
                if branchesNow[b].daughterId[0][1..-1].to_i > 0 # This means only if they daughter ids are not 'none' or 'unborn'.
                  # Do another fast id search to find the two daughters.
                  lower = -1
                  upper = @branch.size
                  d0Found = false
                  until d0Found
                    bb = lower+(upper-lower)/2
                    if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[0][1..-1].to_i
                      lower = bb
                    elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[0][1..-1].to_i
                      upper = bb
                    elsif @branch[bb].id == branchesNow[b].daughterId[0]
                      d0Found = true
                      branchesNext << @branch[bb]
                      @branch[bb].addStartSeq(branchesNow[b].endSeq)
                      next
                    else
                      raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
                    end
                  end
                end
                if branchesNow[b].daughterId[1][1..-1].to_i > 0 # This means only if they daughter ids are not 'none' or 'unborn'.
                  # Do another fast id search to find the two daughters.
                  lower = -1
                  upper = @branch.size
                  d1Found = false
                  until d1Found
                    bb = lower+(upper-lower)/2
                    if @branch[bb].id[1..-1].to_i < branchesNow[b].daughterId[1][1..-1].to_i
                      lower = bb
                    elsif @branch[bb].id[1..-1].to_i > branchesNow[b].daughterId[1][1..-1].to_i
                      upper = bb
                    elsif @branch[bb].id == branchesNow[b].daughterId[1]
                      d1Found = true
                      branchesNext << @branch[bb]
                      @branch[bb].addStartSeq(branchesNow[b].endSeq)
                      next
                    else
                      raise "Some problem occurred with the fast daughter search!" # This should never be called at all.
                    end
                  end
                end
              end
              branchesNow = branchesNext
              change = true if branchesNext.size > 0

              # Write the threadCount to the respective threadCountFile.
              threadCountFile = File.open(threadCountFileNames[t],"w")
              threadCountFile.write("#{threadCount},r")
              threadCountFile.close

            end  # until change == false

            # Delete all branches from array @branch that have no endSeq. Just to speed up writing, reading, and merging.
            @branch.size.times do |b|
              @branch[b] = nil if @branch[b].endSeq == nil
            end
            @branch.compact!

            # Use Marshal.dump to write the remaining branches to a dump file.
            branchPackageDumpFile = File.new(branchPackageDumpFileNames[t],"w")
            branchPackageDumpFile.write(Marshal.dump(@branch))

            # Write the threadCount to the respective threadCountFile.
            threadCountFile = File.open(threadCountFileNames[t],"w")
            threadCountFile.write("#{threadCount},f")
            threadCountFile.close

            # This process can now quit.
            exit
          end # Process.fork
        end # @threads.times

        # Allow some time until the first threadCount files have been written.
        sleep 0.2

        # Continuously monitor the threadCount files, get the threadCount and see whether processes are still running.
        allProcessesFinished = false
        until allProcessesFinished
          sleep 0.1
          threadCounts = []
          runCodes = []
          threadCountFileNames.each do |t|
            tcf = File.open(t,"r")
            tcfString = tcf.read
            tcf.close
            tcfAry = tcfString.split(",")
            threadCounts << tcfAry[0].to_i
            runCodes << tcfAry[1]
          end
          processesRunning = 0
          runCodes.each do |rc|
            processesRunning += 1 if rc == "r"
          end
          allProcessesFinished = true unless runCodes.include?("r")
          if verbose
            if @threads - processesRunning == 0
              print "\rEvolving the ancestral nucleotide sequence along the tree... (#{count+threadCounts.sum}/#{branchSizeBefore} branches)"
            elsif @threads - processesRunning == 1
              print "\rEvolving the ancestral nucleotide sequence along the tree... (#{count+threadCounts.sum}/#{branchSizeBefore} branches - 1/#{@threads} process finished)"
            else
              print "\rEvolving the ancestral nucleotide sequence along the tree... (#{count+threadCounts.sum}/#{branchSizeBefore} branches - #{@threads - processesRunning}/#{@threads} processes finished)"
            end
          end
        end

        # Read all dump files.
        print "\r                                                                                                                     " if verbose
        print "\rMerging processes..." if verbose
        branchPackages = []
        branchPackageDumpFileNames.each do |bp|
          branchPackages << Marshal.load(File.read(bp))
        end

        # Merge all branch packages into array @branch again.
        @branch = []
        ids = []
        branchesPrevious.each do |bpv|
          @branch << bpv
          ids << bpv.id
        end
        branchPackages.each do |bp|
          bp.each do |b|
            unless ids.include?(b.id)
              unless b.endSeq == nil
                @branch << b
                ids << b.id
              end
            end
          end
        end

        # Sort the @branch array again according to id.
        branchesSorted = false
        until branchesSorted
          allSortedSoFar = true
          (@branch.size-1).times do |b|
            if @branch[b].id[1..-1].to_i > @branch[b+1].id[1..-1].to_i
              @branch[b],@branch[b+1] = @branch[b+1],@branch[b]
              allSortedSoFar = false
            end
          end
          branchesSorted = true if allSortedSoFar == true
        end

        # Make sure the array @branch has the same properties that it had before the multi-threaded sequence evolution operation.
        raise "After sequence evolution the array '@branch' is not as large as it was before (now: #{@branch.size}, before: #{branchSizeBefore})!" unless @branch.size == branchSizeBefore
        idsAfter = []
        @branch.each do |b|
          idsAfter << b.id
        end
        raise "After sequence evolution, the array '@branch' has different ids than it had before!" unless idsAfter == idsBefore

        # Clean up: Delete files branchPackageDumpFileNames and threadCountFileNames.
        branchPackageDumpFileNames.each {|bf| File.delete(bf) if File.exists?(bf)}
        threadCountFileNames.each {|tf| File.delete(tf) if File.exists?(tf)}

      end

    else
      raise "The specified number of threads can't be understood!"
    end # if @threads == 1, elsif @threads.class == Fixnum and @threads > 1, else, ...

    @sequencesEvolved = true
    @sequencesMasked = false
    if verbose
      endTime = Time.now
      puts "\rNucleotide sequences added to all branches.                                                               "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end

  end

  def assignAlignmentPositions(verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "--------------------------------------------------phylsim.rb | Sequences----------------------------------------------------"
      puts 
    end

    # Assigning alignment positions.
    print "Assigning alignment positions to all branches..." if verbose
    if @focusGroup # If there is a focus group, all focus group branches should come first in the alignment.

      # Find the total number of focus group branches (including ancestral branches).
      totalFocusGroupSize = 0
      @branch.each {|b| totalFocusGroupSize += 1 if b.focusGroup}
      raise "No focus group branches found even though a focus group seems to exist!" if totalFocusGroupSize == 0
      
      # Fill two arrays, one with position number for the focus group, one with position number for all remaining branches.
      focusGroupPositionNumbers = []
      totalFocusGroupSize.times do |i|
        focusGroupPositionNumbers << i
      end
      nonFocusGroupPositionNumbers = []
      totalFocusGroupSize.upto(@branch.size - 1) do |i|
        nonFocusGroupPositionNumbers << i
      end

      # Make sure, the sum of both array sizes matches the size of array @branch.
      raise "The sizes of focus group and non-focus group position arrays do not sum to the size of the array including all branches" unless @branch.size == focusGroupPositionNumbers.size + nonFocusGroupPositionNumbers.size

      # Shuffle both arrays.
      focusGroupPositionNumbers.shuffle!
      nonFocusGroupPositionNumbers.shuffle!

      # Assign numbers from both arrays to all branches.
      @branch.each do |b|
        if b.focusGroup
          b.addAlignmentPosition(alignmentPosition = focusGroupPositionNumbers.shift)
        else
          b.addAlignmentPosition(alignmentPosition = nonFocusGroupPositionNumbers.shift)
        end
      end

    else # If there's no focus group, all positions can be assigned at random.
      
      # Fill just one array with position numbers.
      positionNumbers = []
      @branch.size.times do |i|
        positionNumbers << i
      end

      # Shuffle this array.
      positionNumbers.shuffle!

      # Assign numbers from this array to all branches.
      @branch.each {|b| b.addAlignmentPosition(alignmentPosition = positionNumbers.shift)}

    end

    # Setting variable @alignmentPositionsAssigned
    @alignmentPositionsAssigned = true

    # Report time consumed.
    if verbose
      endTime = Time.now
      puts "\rAlignment positions assigned to all branches.                                                               "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
  end

  def maskSequences(fileName, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "--------------------------------------------------phylsim.rb | Sequences----------------------------------------------------"
      puts 
    end

    # Make sure this tree has sequences at all.
    raise "Sequences must be evolved before they can be masked!" unless @sequencesEvolved

    # Make sure sequences of this tree have not been removed since they've last been evolved.
    raise "Sequences of this tree have previously been removed. Before removing further sequences, they must be evolved again!" if @sequencesMasked

    # Make sure alignment positions have been assigned before masking sequences.
    self.assignAlignmentPositions unless @alignmentPositionsAssigned

    # Make sure the mask file exists.
    until File.exists?(fileName)
      puts "The file '#{fileName}' could not be found. Please specify a new file name:"
      fileName = gets
    end

    # Read the mask file.
    print "Reading '#{fileName}'..." if verbose
    file = File.open(fileName)
    lines = file.readlines
    file.close
    maskIds = []
    maskSequences = []

    # Check the file format.
    if lines[0] =~ /^>.*/
      # Assume aligned fasta format.
      lines.size.times do |l|
        if lines[l][0] == ">"
          maskIds << lines[l].strip[1..-1]
          maskSequences << ""
        elsif lines[l].strip != ""
          maskSequences.last << lines[l]
        end
      end
      claimedNumberOfTaxa = maskIds.size
      claimedNumberOfCharacters = maskSequences[0].length

    elsif lines[0].downcase =~ /^#nexus/
      # Assume nexus format.
      claimedNumberOfTaxa = 0
      claimedNumberOfCharacters = 0
      lines.size.times do |l|
        lines[l].strip!
        lines[l].gsub!("\t","  ")
      end
      matrix = false
      lines.size.times do |l|
        if lines[l].downcase =~ /dimensions ntax=(\d+) nchar=(\d+)/
          claimedNumberOfTaxa = $1.to_i
          claimedNumberOfCharacters = $2.to_i
        end
        matrix = true if lines[l-1].strip[0..5].downcase == "matrix"
        matrix = false if lines[l].strip == ";"
        if matrix == true
          linesAry = lines[l].strip.split(" ")
          raise "A line of the Nexus data matrix in file '#{fileName}' could not be read correctly: #{lines[l]}: linesAry.size = #{linesAry.size}" unless linesAry.size == 2
          maskIds << linesAry[0]
          maskSequences << linesAry[1]
        end      
      end
      raise "The number of taxa could not be read from file '#{fileName}'!" if claimedNumberOfTaxa == 0
      raise "The number of characters could not be read from file '#{fileName}'!" if claimedNumberOfCharacters == 0

    elsif lines[0] =~ /^(\d+) (\d+)/
      # Assume phylip format.
      claimedNumberOfTaxa = $1.to_i
      claimedNumberOfCharacters = $2.to_i
      1.upto(lines.size-1) do |l|
        unless lines[l].strip == ""
          linesAry = lines[l].strip.split(" ")
          maskIds << linesAry[0]
          maskSequences << linesAry[1]
        end
      end

    else
      raise "File format of '#{fileName}' could not be recognized!"
    end

    # Make sure the number of mask ids and mask sequences is ok.
    raise "The number of ids does not equal the number of sequences in the mask alignment!" unless maskIds.size == maskSequences.size
    raise "The number of sequences found in the mask alignment (#{maskSequences.size}) does not match the number of taxa specified in '#{fileName}'(#{claimedNumberOfTaxa})!" unless maskSequences.size == claimedNumberOfTaxa

    # Make sure the sequences are all of the length that has been specified in the nexus or phylip file.
    maskSequences.each do |s|
      raise "The length of sequence #{maskIds[maskSequences.index(s)]} differs from the number of characters specified in '#{fileName}'" unless s.length == claimedNumberOfCharacters
    end

    # Make sure the number of sequences in the mask file matches the number of extant branches in the (reconstructed) tree.
    raise "The number of sequences in the mask alignment does not match the number of extant species!" unless maskIds.size == @Np

    # Make sure the length of sequences in the mask alignment matches that of branches
    raise "The length of sequences in the mask alignment does not match that of branches!" unless @branch[0].endSeq.length == maskSequences[0].length

    # Create an array of only extant branches which is sorted according to their alignment position.
    alignmentBranch = []
    @branch.each {|b| alignmentBranch[b.alignmentPosition] = b if b.extant}
    alignmentBranch.compact!
    
    # Make sure the size of the array alignmentBranch matches the size of sequences in the mask alignment.
    raise "The size of the array with sorted extant branches (#{alignmentBranch.size}) does not match that of sequences in the mask alignment (#{maskSequences.size})!" unless alignmentBranch.size == maskSequences.size

    # Create a temporary sequence from the end sequence of each branch in array alignmentBranch, then copy all gaps and missing data from the mask sequence onto it, and update endSeq of that branch.
    alignmentBranch.size.times do |b|
      tempSequence = alignmentBranch[b].endSeq
      tempSequence.length.times do |p|
        tempSequence[p] = "-" if maskSequences[b][p] == "-"
        tempSequence[p] = "?" if maskSequences[b][p] == "?"
        tempSequence[p] = "n" if maskSequences[b][p] == "n"
      end
      alignmentBranch[b].updateEndSeq(tempSequence)
    end

    # Set variable @sequencesMasked.
    @sequencesMasked = true

    # Report time consumed.
    if verbose
      endTime = Time.now
      puts "\rEnd sequences of extant branches have been masked according to the alignment in '#{fileName}'."
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
  end

  def length(unit = "duration", application = "/usr/bin/raxmlLight-PTHREADS-SSE3", threads = 7, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "-----------------------------------------------------phylsim.rb | Info------------------------------------------------------"
      puts 
    end

    # Check the unit argument.
    raise "Either 'duration', 'expectedNumberOfSubstitutionsPerSite', 'actualNumberOfSubstitutions', or 'actualNumberOfSubstitutionsPerSite' must be specified for parameter 'unit'!" unless ["duration","expectednumberofsubstitutionspersite","actualnumberofsubstitutions","actualnumberofsubstitutionspersite","inferrednumberofsubstitutionspersite"].include?(unit.downcase)
    raise "In order to measure the tree length in expected substitutions per site, branch rates must first be assigned!" if unit.downcase == "expectednumberofsubstitutionspersite" and @branchRatesAssigned == false
    raise "In order to measure the tree length in actual substitutions, sequences must first be evolved!" if unit.downcase == "actualnumberofsubstitutions" and @sequencesEvolved != true
    raise "In order to measure the tree length actual substitutions per site, sequences must first be evolved!" if unit.downcase == "actualnumberofsubstitutionspersite" and @sequencesEvolved != true
    raise "In order to infer the tree length actual substitutions per site, sequences must first be evolved!" if unit.downcase == "inferrednumberofsubstitutionspersite" and @sequencesEvolved != true

    if unit.downcase == "duration"
      treeLength = 0
      @branch.each do |b|
        treeLength += b.duration
      end
    elsif unit.downcase == "expectednumberofsubstitutionspersite"
      treeLength = 0
      @branch.each do |b|
        treeLength += b.expectedNumberOfSubstitutionsPerSite
      end
    elsif unit.downcase == "actualnumberofsubstitutions"
      treeLength = 0
      @branch.each do |b|
        treeLength += b.actualNumberOfSubstitutions
      end
    elsif unit.downcase == "actualnumberofsubstitutionspersite"
      treeLength = 0
      @branch.each do |b|
        treeLength += b.actualNumberOfSubstitutionsPerSite
      end
    elsif unit.downcase == "inferrednumberofsubstitutionspersite"

      # Check the application argument.
      until File.exists?(application)
        puts "Application '#{application}' could not be found. Please specify a new application name:"
        application = gets
      end

      # Define the name of a temporary folder for the analysis.
      random = rand(100000000000000)
      tempDirName = ".tmp" + random.to_s
      until File.exists?(tempDirName) == false
        random = rand(100000000000000)
        tempDirName = ".tmp" + random
      end

      # Initiate treeLength.
      treeLength = 0

      # Create the temporary folder.
      Dir.mkdir(tempDirName)

      # Move into the temporary folder.
      Dir.chdir(tempDirName)
      
      # Define the name of temporary input and output files.
      tempInFileName = "in.phy"
      tempTreeFileName = "in.tre"
      tempOutFileName = "out.txt"

      # Write an alignment of all sequences in phylip format.
      print "Writing sequence alignment in phylip format to a temporary folder..." if verbose
      self.writeSequences(fileName = tempInFileName, format = "phylip", false, false, false, false)

      # Run the tree inference software and read its result.
      resultString = ""
      if application.downcase.include?("raxml")
        print "\rRunning #{application} to infer the number of substitutions per site..." if verbose
        self.to_newick(fileName = tempTreeFileName, branchLengths = "actualNumberOfSubstitutionsPerSite", false, true, false, false, false)
        system("#{application} -s #{tempInFileName} -n #{tempOutFileName} -T #{threads} -m GTRGAMMA -t #{tempTreeFileName} > /dev/null")
        resultFileName = "RAxML_result.out.txt"
        raise "File #{resultFileName} does not exist!" unless File.exists?(resultFileName)
        resultFile = File.open(resultFileName)
        resultString = resultFile.read
        resultString.strip!
        resultString.chomp!(";")
        resultAry = resultString.split(":")
        resultAry.each {|i| treeLength += $1.to_f if i =~ /^(\d+\.\d+)[,\)].*$/}
      else
        raise "Application #{application} is not supported!"
      end

      # Delete all files in the temporary directory.
      Dir.entries(".").each {|f| File.delete(f) unless File.directory?(f)}

      # Move out of the temporary folder into the previous working directory.
      Dir.chdir("../")

      # Remove the temporary directory.
      Dir.rmdir(tempDirName)

    end

    if verbose
      endTime = Time.now
      puts "\rTree length inferred: #{treeLength.round(3)}.                                                                               "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end

    # Return the resulting tree length.
    treeLength
  end

  def sPScore(includeEmpty = false, match = 1, mismatch = -2, gapOpening = -3, gapExtension = -2, threads = 7, verbose = true)
    # SP stands for 'sum of all pairs' (Wang et al. 1994)
    if verbose
      startTime = Time.now
      puts
      puts "-----------------------------------------------------phylsim.rb | Info------------------------------------------------------"
      puts 
    end

    # Feedback
    print "\rCalculating SP score..." if verbose

    # Make sure that nucleotide sequences have been added to all branches.
    raise "In order to calculate the SP score, nucleotide sequences must first be added to all branches (method 'evolveSequences')!" unless @sequencesEvolved

    # Make sure, match, mismatch, gapOpening, and gapExtension are either Fixnum or Float.
    raise "Please specify a number for parameter match!" unless [Fixnum,Float].include?(match.class)
    raise "Please specify a number for parameter mismatch!" unless [Fixnum,Float].include?(mismatch.class)
    raise "Please specify a number for parameter gapOpening!" unless [Fixnum,Float].include?(gapOpening.class)
    raise "Please specify a number for parameter gapExtension!" unless [Fixnum,Float].include?(gapExtension.class)

    # Check the number of threads
    if threads.class == Fixnum and threads > 1
      @threads = threads
    else
      @threads = 1
    end

    # Collect sequences of extant branches.
    sequences = []
    if includeEmpty
      @branch.each {|b| sequences << b.endSeq.dup if b.extant}
    else
      @branch.each {|b| sequences << b.endSeq.dup if b.extant and b.endSeq.gsub("-","") != ""}
    end

    # Get the alignment length.
    alignmentLength = sequences[0].length

    # Define file names to which the threadCount will be written.
    threadSPScoreFileNames = []
    @threads.times do |t|
      threadSPScoreFileNames << ".threadSPScore#{(t+1).to_s.rjust((@threads+1).to_s.length).gsub(" ","0")}.txt"
    end

    # Define which part of the alignment should be analysed by each thread.
    threadStartPos = [0]
    threadEndPos = []
    if @threads > 1
      if (alignmentLength/@threads)*@threads == alignmentLength
        partLength = alignmentLength/@threads
      else
        partLength = alignmentLength/@threads + 1
      end
      1.upto(@threads-1) do |t|
        threadEndPos << partLength*t - 1
        threadStartPos << partLength*t
      end
    end
    threadEndPos << alignmentLength - 1

    # Define bases and simplify all sequences
    bases = ["a","c","g","t"]
    sequences.size.times do |s|
      sequences[s] = sequences[s].gsub("?","n")
      sequences[s] = sequences[s].gsub("n","n")
      sequences[s] = sequences[s].gsub("k","n")
      sequences[s] = sequences[s].gsub("m","n")
      sequences[s] = sequences[s].gsub("r","n")
      sequences[s] = sequences[s].gsub("y","n")
      sequences[s] = sequences[s].gsub("s","n")
      sequences[s] = sequences[s].gsub("w","n")
      sequences[s] = sequences[s].gsub("b","n")
      sequences[s] = sequences[s].gsub("v","n")
      sequences[s] = sequences[s].gsub("h","n")
      sequences[s] = sequences[s].gsub("d","n")
      sequences[s] = sequences[s].gsub("x","n")
    end

    # Start multi-threaded score counting.
    sPScore = 0
    @threads.times do |t|
      Process.fork do

        # Calculate partial SP score.
        threadSPScore = 0
        0.upto(sequences.size-2) do |i|
          i.upto(sequences.size-1) do |j|
            threadStartPos[t].upto(threadEndPos[t]) do |p|
              if sequences[i][p] == "-" or sequences[j][p] == "-"
                if p == 0
                  threadSPScore += gapOpening
                else
                  if sequences[i][p] == "-" and sequences[i][p-1] == "-"
                    threadSPScore += gapExtension
                  elsif sequences[j][p] == "-" and sequences[j][p-1] == "-"
                    threadSPScore += gapExtension
                  else
                    threadSPScore += gapOpening
                  end
                end
              elsif sequences[i][p] == "n" or sequences[j][p] == "n"
                threadSPScore += mismatch
              elsif bases.include?(sequences[i][p]) and bases.include?(sequences[j][p])
                if sequences[i][p] == sequences[j][p]
                  threadSPScore += match
                elsif
                  threadSPScore += mismatch
                end
              else
                raise "Unexpected base combination found: '#{sequences[i][p]}' and '#{sequences[j][p]}'!"
              end
            end
          end
        end

        # Write resulting score to file.
        threadSPScoreFile = File.new(threadSPScoreFileNames[t],"w")
        threadSPScoreFile.write(threadSPScore)
        threadSPScoreFile.close

        # The thread can now exit
        exit
      end
    end

    # Wait for all processes to finish.
    Process.waitall

    # Read the partial sp scores from all result files and add it to the overall sp score.
    @threads.times do |t|
      threadSPScoreFile = File.open(threadSPScoreFileNames[t])
      threadSPScore = threadSPScoreFile.read.strip.to_i
      threadSPScoreFile.close
      sPScore += threadSPScore
    end

    # Clean up.
    threadSPScoreFileNames.each {|t| File.delete(t) if File.exists?(t)}

    if verbose
      endTime = Time.now
      puts "\rSP score inferred: #{sPScore}.                                                                               "
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end

    # Return the resulting SP-score
    sPScore
  end

  def writeTraits(fileName = "traits.txt", includeAncestral = false, overwrite = false, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "----------------------------------------------------phylsim.rb | Export-----------------------------------------------------"
      puts 
    end

    # Make sure that traits have been added to all branches.
    raise "In order to export traits, these must first be added to all branches (method 'evolveTraits')!" unless @traitsEvolved

    # Make sure that alignment positions have been added to all branches. If not, do so now.
    self.assignAlignmentPositions(verbose = false) unless @alignmentPositionsAssigned

    # Check the arguments.
    unless overwrite
      # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
      if File.exists?(fileName)
        puts "Ok to replace '#{fileName}'? (y/N)"
        answer = gets
        unless answer[0].downcase == "y"
          puts "Please specify a new filename:"
          fileName = gets.strip
        end
        puts
      end
    end

    # Feedback.
    print "Writing traits..." if verbose

    # Create ids and trait arrays and fill them.
    ids = []
    traits = []
    @branch.each do |b|
      if b.extant == false and includeAncestral == false
        ids[b.alignmentPosition] = nil
        traits[b.alignmentPosition] = nil
      else
        if b.extant
          ids[b.alignmentPosition] = b.speciesId.split("/")[-1]
        else
          ids[b.alignmentPosition] = b.id
        end
        traits[b.alignmentPosition] = b.endTrait
      end
    end
    ids.compact!
    traits.compact!

    # Prepare output.
    output = "Species\tTrait\n"
    ids.size.times do |i|
      output << "#{ids[i]}\t#{traits[i]}\n"
    end

    # Write output.
    file = File.open(fileName,"w")
    file.write(output)
    file.close

    # Write time consumed.
    if verbose
      endTime = Time.now
      puts "\rTraits written to #{fileName}."
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
  end

  def writeSequences(fileName = "sequences.fasta", format = "fasta", includeAncestral = false, includeEmpty = true, overwrite = false, verbose = true)
    if verbose
      startTime = Time.now
      puts
      puts "----------------------------------------------------phylsim.rb | Export-----------------------------------------------------"
      puts 
    end

    # Make sure that nucleotide sequences have been added to all branches.
    raise "In order to export to fasta format, nucleotide sequences must first be added to all branches (method 'evolveSequences')!" unless @sequencesEvolved

    # Make sure that alignment positions have been added to all branches. If not, do so now.
    self.assignAlignmentPositions(verbose = false) unless @alignmentPositionsAssigned

    # Make sure sequences are not masked and ancestral sequences are included.
    raise "After extant sequences have been masked, including ancestral sequences in the alignment may lead to conflicts!" if includeAncestral and @sequencesMasked

    # Check the arguments.
    unless overwrite
      # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
      if File.exists?(fileName)
        puts "Ok to replace '#{fileName}'? (y/N)"
        answer = gets
        unless answer[0].downcase == "y"
          puts "Please specify a new filename:"
          fileName = gets.strip
        end
        puts
      end
    end

    # Feedback.
    print "Writing nucleotide sequences..." if verbose

    # Create ids and sequences arrays and fill them.
    ids = []
    seqs = []
    @branch.each do |b|
      if b.extant == false and includeAncestral == false
        ids[b.alignmentPosition] = nil
        seqs[b.alignmentPosition] = nil
      elsif b.endSeq.gsub("-","").gsub("?","").gsub("n","").gsub("N","") == "" and includeEmpty == false
        ids[b.alignmentPosition] = nil
        seqs[b.alignmentPosition] = nil
      else
        if b.extant
          ids[b.alignmentPosition] = b.speciesId.split("/")[-1]
        else
          ids[b.alignmentPosition] = b.id
        end
        seqs[b.alignmentPosition] = b.endSeq
      end
    end
    ids.compact!
    seqs.compact!

    # Prepare output.
    if format.downcase == "fasta"
      output = ""
      ids.size.times do |i|
        output << ">#{ids[i]}\n"
        output << "#{seqs[i].slice!(0..59)}\n" while seqs[i].length > 0
      end
    elsif format.downcase == "nexus"
      maxIdLength = 0
      ids.each do |i|
        maxIdLength = i.length if i.length > maxIdLength
      end
      output = "#NEXUS\n"
      output << "begin data;\n"
      output << "  dimensions ntax=#{ids.size} nchar=#{seqs[0].length};\n"
      output << "  format datatype = dna gap = - missing = ?;\n"
      output << "  matrix\n"
      ids.size.times do |i|
        output << "  #{ids[i].ljust(maxIdLength+2)}#{seqs[i]}\n"
      end
      output << "  ;\n"
      output << "end;\n"      
    elsif format.downcase == "phylip"
      maxIdLength = 0
      ids.each do |i|
        maxIdLength = i.length if i.length > maxIdLength
      end
      output = "#{ids.size} #{seqs[0].length}\n"
      ids.size.times do |i|
        output << "#{ids[i].ljust(maxIdLength+2)}#{seqs[i]}\n"
      end
    else
      raise "Format '#{format}' is not known!"
    end

    # Write output.
    file = File.open(fileName,"w")
    file.write(output)
    file.close

    # Write time consumed.
    if verbose
      endTime = Time.now
      puts "\rNucleotide sequences written in #{format.downcase} format to #{fileName}."
      if endTime-startTime < 60
        puts "Time used: #{(endTime-startTime).round(2)} seconds."
      elsif endTime-startTime < 3600
        puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
      else
        puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
      end
      puts
    end
  end

  def info(fileName = nil, overwrite = false, verbose = true)
    # Just a report of all the settings that have been used to produce the current tree.

    if verbose
      startTime = Time.now
      puts
      puts "-----------------------------------------------------phylsim.rb | Info------------------------------------------------------"
      puts
    end

    # Check the arguments.
    unless overwrite or fileName == nil
      # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
      if File.exists?(fileName)
        puts "Ok to replace '#{fileName}'? (y/N)"
        answer = gets
        unless answer[0].downcase == "y"
          puts "Please specify a new filename:"
          fileName = gets.strip
        end
        puts
      end
    end

    # Prepare the output.
    output = "Tree generation:\n"
    output << " lambda:                                #{@posteriorLambda.round(4)}\n"
    output << " mu:                                    #{@posteriorMu.round(4)}\n"
    output << " treeOrigin:                            #{@posteriorTreeOrigin.round(4)}\n"
    output << " Number of extant species:              #{@NpFull}\n"
    output << " Sum of species durations:              #{@sumOfSpeciesDurations.round(4)}\n"
    output << "\n"
    if @treeReconstructed
      output << "Tree reconstruction:\n"
      output << " Sampling scheme:                       #{@samplingScheme}\n"
      output << " Number of sampled species:             #{@Np}\n"
      if @focusGroup != nil
        output << " Number of sampled focus group species: #{@focusGroup[2]}\n"
        output << " Number of extant focus group species:  #{@focusGroupNpFull}\n"
        output << " Focus group age:                       #{@focusGroupAge.round(4)}\n"
      end
      output << " Sampled phylogenetic diversity:        #{self.phylogeneticDiversity.round(4)}\n"
      output << " Gamma statistic:                       #{@gamma.round(4)}\n"
      output << "\n"
    end
    if @branchRatesAssigned
      output << "Branch rates:\n"
      output << " Mean:                                  #{@branchRatesMean}\n"
      output << " Standard deviation:                    #{@branchRatesStandardDeviation}\n"
      unless @branchRatesAutoCorrelation == nil or @branchRatesAutoCorrelation == false or @branchRatesAutoCorrelation == 0
        output << " Autocorrelation:                       #{@branchRatesAutoCorrelation}\n"
      end
      output << "\n"
    end
    if @sequencesEvolved
      output << "Sequence evolution:\n"
      output << " Substitution model class:              #{@substitutionModel.class}\n"
      if @substitutionModel.class == ECM
        output << " 61 codon frequencies and 1830 substitution rates as specified in '#{@ecmFileName}'\n"
      elsif @substitutionModel.class == GTR
        output << " pi(A):                                 #{@substitutionModel.pi[0].round(4)}\n"
        output << " pi(C):                                 #{@substitutionModel.pi[1].round(4)}\n"
        output << " pi(G):                                 #{@substitutionModel.pi[2].round(4)}\n"
        output << " pi(T):                                 #{@substitutionModel.pi[3].round(4)}\n"
        output << " rate(A<->C):                           #{@substitutionModel.s[0].round(4)}\n"
        output << " rate(A<->G):                           #{@substitutionModel.s[1].round(4)}\n"
        output << " rate(A<->T):                           #{@substitutionModel.s[2].round(4)}\n"
        output << " rate(C<->G):                           #{@substitutionModel.s[3].round(4)}\n"
        output << " rate(C<->T):                           #{@substitutionModel.s[4].round(4)}\n"
        output << " rate(G<->T):                           #{@substitutionModel.s[5].round(4)}\n"
        output << " alpha:                                 #{@substitutionModel.alpha.round(4)}\n"
      else
        raise "Parameters of this substitution model class can't yet be written to info file!"
      end
      output << "\n"
    end

    # Write the output to screen or file.
    if fileName == nil
      output
    else
      File.new(fileName,"w").puts output
      if verbose
        endTime = Time.now
        puts "\rModel information written to #{fileName}."
        if endTime-startTime < 60
          puts "Time used: #{(endTime-startTime).round(2)} seconds."
        elsif endTime-startTime < 3600
          puts "Time used: #{((endTime-startTime)/60.0).round(2)} minutes."
        else
          puts "Time used: #{((endTime-startTime)/3600.0).round(2)} hours."
        end
        puts
      end
    end

  end

  def outgroup(fileName = "outgroup.txt", overwrite = false)
    unless overwrite
      # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
      if File.exists?(fileName)
        puts "Ok to replace '#{fileName}'? (y/N)"
        answer = gets
        unless answer[0].downcase == "y"
          puts "Please specify a new filename:"
          fileName = gets.strip
        end
        puts
      end
    end

    # Prepare the outgroup string.
    outgroupString = ""
    if @branch[0].extant
      outgroupString << "#{@branch[0].speciesId.split("/")[-1]}"
    elsif @branch[1].extant
      outgroupString << "#{@branch[1].speciesId.split("/")[-1]}"
    else
      # Make sure that the sum of the extant progeny of @branch[0] and the extant progeny of @branch[1] equals the number of extant taxa.
      extantSum = 0
      extantSum += @branch[0].extantProgenyId.size
      extantSum += @branch[1].extantProgenyId.size
      raise "The sum of the sizes of extant progenies (#{@branch[0].extantProgenyId.size} + #{@branch[1].extantProgenyId.size} = #{extantSum}) does not match the number of extant taxa (#{@Np})!" unless extantSum == @Np

      # Find out which of the two root branches' extant progenies is the smaller one, then, write all progeny species ids to the outgroup string.
      if @branch[0].extantProgenyId.size <= @branch[1].extantProgenyId.size
        @branch.each {|b| outgroupString << "#{b.speciesId.split("/")[-1]}," if @branch[0].extantProgenyId.include?(b.id)}
      else
        @branch.each {|b| outgroupString << "#{b.speciesId.split("/")[-1]}," if @branch[1].extantProgenyId.include?(b.id)}
      end
      outgroupString.chomp!(",")
    end

    # Write outgroup string to stdout or to file 'fileName'
    if fileName == nil
      return outgroupString
    else
      file = File.open(fileName,"w")
      file.write(outgroupString)
      file.close
    end
  end

  def constraint(fileName = "constraint.txt", overwrite = false)
    unless overwrite
      # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
      if File.exists?(fileName)
        puts "Ok to replace '#{fileName}'? (y/N)"
        answer = gets
        unless answer[0].downcase == "y"
          puts "Please specify a new filename:"
          fileName = gets.strip
        end
        puts
      end
    end

    # Figure out which ids go into the two constraint groups.
    group0 = []
    group1 = []
    if @branch[0].extant
      group0 << @branch[0].speciesId.split("/")[-1]
      @branch[1..-1].each {|b| group1 << b.speciesId.split("/")[-1] if b.extant}
    elsif @branch[1].extant
      group1 << @branch[1].speciesId.split("/")[-1]
      @branch[2..-1].each {|b| group0 << b.speciesId.split("/")[-1] if b.extant}
    else
      @branch.each {|b| group0 << b.speciesId.split("/")[-1] if @branch[0].extantProgenyId.include?(b.id)}
      @branch.each {|b| group1 << b.speciesId.split("/")[-1] if @branch[1].extantProgenyId.include?(b.id)}
    end
    raise "The sum of species ids included in both constraint groups (#{group0.size} + #{group1.size} = #{group0.size+group1.size}) does not match the number of extant species #{@Np}!" unless group0.size+group1.size == @Np

    # Prepare the constraint string.
    constraintString = "(("
    group0.each {|g| constraintString << "#{g},"}
    constraintString.chomp!(",")
    constraintString << "),("
    group1.each {|g| constraintString << "#{g},"}
    constraintString.chomp!(",")
    constraintString << "));\n"

    # Write outgroup string to stdout or to file 'fileName'
    if fileName == nil
      return constraintString
    else
      file = File.open(fileName,"w")
      file.write(constraintString)
      file.close
    end
  end

  def monophyleticGroups(fileName = "monophyleticGroups.txt", overwrite = false)
    unless overwrite
      # If the file with name 'fileName' exists in the current directory, ask the user whether it's ok to overwrite it.
      if File.exists?(fileName)
        puts "Ok to replace '#{fileName}'? (y/N)"
        answer = gets
        unless answer[0].downcase == "y"
          puts "Please specify a new filename:"
          fileName = gets.strip
        end
        puts
      end
    end

    # For the id search below, prepare a sorted array of extant branches.
    extantBranch = []
    @branch.each {|b| extantBranch << b if b.extant}

    # For each branch, take its sampled extant progeny as one monophyletic group (after reconstruction, every extant species has been sampled).
    monophyleticGroups = []
    @branch.each do |b|
      unless b.extantProgenyId.size == 0
        monophyleticGroupOfThisBranch = []
        b.extantProgenyId.each do |i|
          idFound = false
          extantBranch.each do |eb|
            if eb.id == i
              monophyleticGroupOfThisBranch << eb.speciesId.split("/")[-1]
              idFound = true
            end
          end
          raise "Branch id #{i} was not found!" if idFound == false
        end
        monophyleticGroups << monophyleticGroupOfThisBranch.sort
      end
    end

    # Prepare monophyletic groups string.
    monophyleticGroupsString = ""
    monophyleticGroups.each do |m|
      m.each {|n| monophyleticGroupsString << "#{n},"}
      monophyleticGroupsString.chomp!(",")
      monophyleticGroupsString << "\n"
    end

    # Write monophyletic groups string to stdout or to file 'fileName'.
    if fileName == nil
      return monophyleticGroupsString
    else
      file = File.open(fileName,"w")
      file.write(monophyleticGroupsString)
      file.close
    end
  end

  def linkWithAlignment(fileName = "alignment.nex", fileType = "nexus", verbose = true)
    # Feedback.
    if verbose
      startTime = Time.now
      puts
      puts "--------------------------------------------------phylsim.rb | Alignment----------------------------------------------------"
      puts
    end

    # Make sure the file exists.
    unless File.exists?(fileName)
      raise "The alignment file \"#{fileName}\" could not be found!"
    end

    # Read the file.
    alignmentFile = File.open(fileName)
    alignmentLines = alignmentFile.readlines
    
    # Initiate variables ntax, nchar, and the arrays ids, and seqs.
    ntax = nil
    nchar = nil
    ids = []
    seqs = []    

    # If Nexus format is specified.
    if fileType == "nexus"
      # Make sure the file is in Nexus format.
      unless alignmentLines[0].strip.downcase[0..5] == "#nexus"
        raise "The alignment file #{fileName} does not seem to be in Nexus format!"
      end

      # Remove empty lines and strip each line.
      compressedAlignmentLines = []
      alignmentLines.each{|l| compressedAlignmentLines << l.strip unless l.strip == ""}

      # Read the specified nchar and ntax.
      compressedAlignmentLines.each do |l|
        if l.match(/dimensions.*ntax=(\d+)\s+nchar=(\d+)/)
          ntax = $1.to_i
          nchar = $2.to_i
        elsif l.match(/dimensions.*nchar=(\d+)\s+ntax=(\d+)/)
          ntax = $2.to_i
          nchar = $1.to_i
        end
      end
      raise "The number of taxa could not be read from the alignment file!" if ntax == nil
      raise "The number of characters could not be read from the alignment file!" if nchar == nil

      # Find the indeces of the lines before and after the sequence matrix.
      beforeIndex = nil
      afterIndex = nil
      compressedAlignmentLines.size.times do |x|
        if compressedAlignmentLines[x].downcase == "matrix"
          beforeIndex = x
        elsif compressedAlignmentLines[x] == ";" and beforeIndex != nil
          afterIndex = x
        end
      end
      raise "The alignment could not be read!" if afterIndex == nil
      matrixLines = compressedAlignmentLines[beforeIndex+1..afterIndex-1]

      # Store the alignment ids and sequences.
      matrixLines.each do |l|
        ids << l.split(/\s+/)[0]
        seqs << l.split(/\s+/)[1]
      end
      
    # If Phylip format has been specified.
    elsif fileType == "phylip"
      # Read the ntax and nchar from the first line of the file.
      ntax = alignmentLines[0].split(/\s+/)[0].to_i
      nchar = alignmentLines[0].split(/\s+/)[1].to_i

      # Store the alignment ids and sequences.
      matrixLines = []
      alignmentLines[1..-1].each{|l| matrixLines << l.strip}
      matrixLines.each do |l|
        ids << l.split(/\s+/)[0]
        seqs << l.split(/\s+/)[1]
      end

    # If another format has been specified.
    else
      raise "The alignment file format \"#{fileType}\" is not supported!"
    end

    # Make sure all sequences are the same length.
    seqs[1..-1].each{|s| raise "Not all sequences of alignment \"#{fileName}\" are of the same length!" if s.length != seqs[0].length}

    # Make sure that ntax and nchar matches the number of ids and the sequence length.      
    raise "The specified number of taxa doesn't match the number of sequences in alignment \"#{fileName}\"!" if ntax != seqs.size
    raise "The specified number of characters doesn't match the sequence length in alignment \"#{fileName}\"!" if nchar != seqs[0].length

    # Make sure that all sequence ids match species ids of branches.
    seqIdsWithoutSpecies = []
    ids.each do |i|
      speciesFound = false
      @branch.each do |b|
        if b.speciesId == i
          speciesFound = true
          break
        end
      end
      seqIdsWithoutSpecies << i if speciesFound == false
    end
    unless seqIdsWithoutSpecies == []
      raiseString = "#{seqIdsWithoutSpecies.size} sequence ids could not be assigned to species. These are "
      seqIdsWithoutSpecies.each{|si| raiseString << "#{si}, "}
      raiseString.chomp(", ")
      raiseString << "!"
      raise raiseString
    end

    # Make sure that sequences are available for all tip species.
    tipSpeciesWithoutSequences = []
    @branch.each do |b|
      unless b.endCause == "speciation"
        sequenceFound = false
        ids.each do |i|
          if b.speciesId == i
            sequenceFound = true
            break
          end
        end
        tipSpeciesWithoutSequences << b.speciesId if sequenceFound == false
      end
    end
    unless tipSpeciesWithoutSequences == []
      raiseString = "Sequences for #{tipSpeciesWithoutSequences.size} species are missing. These species are "
      tipSpeciesWithoutSequences.each{|ts| raiseString << "#{ts}, "}
      raiseString.chomp(", ")
      raiseString << "!"
      raise raiseString
    end

    # Add sequences to branches.
    @branch.each do |b|
      unless b.endCause == "speciation"
        ids.size.times do |x|
          if b.speciesId == ids[x]
            b.addEndSeq(seqs[x])
            break
          end
        end
      end
    end

  end

end

