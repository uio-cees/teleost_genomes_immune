class GTR
  attr_reader :pi, :s, :alpha
  def initialize(string)
    @bases = ["a","c","g","t"]
    @pi = [0,0,0,0]
    @s = [0,0,0,0,0,0]
    @alpha = 0
    if string.downcase =~ /mammal/
      # The following defines normal distribution with mean and standard variation derived from Table S1 (standard deviation is calculated as sqrt(var)) of Arbiza et al. 2011.
      @pi[0] = Rubystats::NormalDistribution.new(0.2538,0.0500).rng until @pi[0] > 0
      @pi[1] = Rubystats::NormalDistribution.new(0.2631,0.0480).rng until @pi[1] > 0
      @pi[2] = Rubystats::NormalDistribution.new(0.2653,0.0361).rng until @pi[2] > 0
      @pi[3] = Rubystats::NormalDistribution.new(0.2178,0.0400).rng until @pi[3] > 0
      @s[0]  = Rubystats::NormalDistribution.new(0.7212,0.2076).rng until @s[0] > 0
      @s[1]  = Rubystats::NormalDistribution.new(2.5220,0.4387).rng until @s[1] > 0
      @s[2]  = Rubystats::NormalDistribution.new(0.5245,0.2078).rng until @s[2] > 0
      @s[3]  = Rubystats::NormalDistribution.new(0.6198,0.2200).rng until @s[3] > 0
      @s[4]  = Rubystats::NormalDistribution.new(3.3510,0.7490).rng until @s[4] > 0
      @s[5]  = Rubystats::NormalDistribution.new(0.5517,0.1881).rng until @s[5] > 0
      if string.downcase =~ /meanalpha/
        @alpha = 0.5376
      else
        @alpha = Rubystats::NormalDistribution.new(0.5376,0.4159).rng until @alpha > 0
      end
      sumPi = @pi[0..3].sum
      @pi.collect! {|p| p = p/sumPi}
      @s.collect! {|s| s/@s[5]}
      raise "Something's wrong with rate G<->T: #{@s[5]}!" unless @s[5] == 1.0

    elsif string.downcase =~ /vertebrat/
      # The following defines normal distribution with mean and standard variation derived from Table S1 (standard deviation is calculated as sqrt(var)) of Arbiza et al. 2011.
      @pi[0] = Rubystats::NormalDistribution.new(0.2519,0.0346).rng until @pi[0] > 0
      @pi[1] = Rubystats::NormalDistribution.new(0.2737,0.0332).rng until @pi[1] > 0
      @pi[2] = Rubystats::NormalDistribution.new(0.2674,0.0224).rng until @pi[2] > 0
      @pi[3] = Rubystats::NormalDistribution.new(0.2070,0.0283).rng until @pi[3] > 0
      @s[0]  = Rubystats::NormalDistribution.new(0.8338,0.1414).rng until @s[0] > 0
      @s[1]  = Rubystats::NormalDistribution.new(2.2900,0.2766).rng until @s[1] > 0
      @s[2]  = Rubystats::NormalDistribution.new(0.6625,0.1637).rng until @s[2] > 0
      @s[3]  = Rubystats::NormalDistribution.new(0.6426,0.1670).rng until @s[3] > 0
      @s[4]  = Rubystats::NormalDistribution.new(3.1980,0.5306).rng until @s[4] > 0
      @s[5]  = Rubystats::NormalDistribution.new(0.5864,0.1249).rng until @s[5] > 0
      if string.downcase =~ /meanalpha/
        @alpha = 0.3658
      else
        @alpha = Rubystats::NormalDistribution.new(0.3658,0.2871).rng until @alpha > 0
      end
      sumPi = @pi[0..3].sum
      @pi.collect! {|p| p = p/sumPi}
      @s.collect! {|s| s/@s[5]}
      raise "Something's wrong with rate G<->T: #{@s[5]}!" unless @s[5] == 1.0

    elsif File.exists?(string)
      lines = File.open(string).readlines
      raise "File #{string} seems to be incorrectly formatted!" unless lines.size == 11 or lines.size == 13
      sLines = lines[0..5]
      sLines.each do |sl|
        @s << sl.split(":")[1].strip.to_f
      end
      piLines = lines[7..10]
      piLines.each do |pl|
        @pi << pl.split(":")[1].strip.to_f
      end
      @alpha = lines[12].split(":")[1].strip.to_f unless lines[12] == nil
    elsif string.split(",").size == 8 or string.split(",").size == 9
      ary = string.split(",")
      @s = [ary[0].to_f,ary[1].to_f,ary[2].to_f,ary[3].to_f,ary[4].to_f,1.0]
      @pi = [ary[5].to_f,ary[6].to_f,ary[7].to_f,(1.0-(ary[5].to_f+ary[6].to_f+ary[7].to_f))]
      @alpha = ary[8].to_f unless ary[8] == nil
    else
      raise "The string specified to initialize the GTR substitution model could not be understood. If it's supposed to be a file name, the file could not be read. If it's supposed to contain parameter values, it must have 8 comma-separated numbers, five rates first, then three base frequencies."
    end
    @numberOfRateCategories = 100
    @r = RSRuby.instance if @r == nil
    @siteSpecificRateMultiplier = @r.rgamma(@numberOfRateCategories,@alpha,@alpha) unless @alpha == 0

    # Check the number of values for s and pi.
    raise "There's a problem with the number of specified substitution rates" unless @s.size == 6
    raise "There's a problem with the number of specified base frequencies" unless @pi.size == 4

    # Make sure all values of pi are within 0 and 1.
    @pi.each do |p|
      raise "All base frequencies must be greater or equal to 0, but less or equal to 1!" unless p >= 0.0 and p <= 1.0
    end

    # Standardize all @pi (this only corrects for rounding errors).
    sum = @pi.sum
    @pi.size.times do |i|
      @pi[i] = @pi[i]/sum
    end
    if @pi.sum != 1.0 # This only corrects rounding errors, e.g. 0.99999999999999 -> 1.0
      @pi[0] = @pi[0] + 1.0 - @pi.sum
    end

    # The constant 'c' ensures that the average rate at equilibrium is 1.
    c = 1.0/(@s[0]*@pi[0]*@pi[1] + @s[1]*@pi[0]*@pi[2] + @s[2]*@pi[0]*@pi[3] + @s[3]*@pi[1]*@pi[2] + @s[4]*@pi[1]*@pi[3] + @s[5]*@pi[2]*@pi[3])

    # Complete the q matrix (see Yang 1994).
    @q = []
    @q[0] = [0              , c*@pi[1]*@s[0] , c*@pi[2]*@s[1] , c*@pi[3]*@s[2]]
    @q[1] = [c*@pi[0]*@s[0] , 0              , c*@pi[2]*@s[3] , c*@pi[3]*@s[4]]
    @q[2] = [c*@pi[0]*@s[1] , c*@pi[1]*@s[3] , 0              , c*@pi[3]*@s[5]]
    @q[3] = [c*@pi[0]*@s[2] , c*@pi[1]*@s[4] , c*@pi[2]*@s[5] , 0             ]
    @q.size.times do |i|
      @q[i][i] = -@q[i].sum
    end
    @qMatrix = Matrix.rows(@q)
    @u_infMatrix, @lMatrix, @uMatrix = @qMatrix.eigensystem
    
  end
  
  def randomSeq(seqLength)
    randomSeq = ""
    seqLength.times do |n|
      random = rand
      if random <= @pi[0]
        randomSeq << @bases[0]
      else
        (@pi.size-1).times do |i|
          if random > @pi[0..i].sum and random <= @pi[0..(i+1)].sum
            randomSeq << @bases[i+1]
            break
          end
        end
        puts "Something wrong about the base frequencies" if random > @pi.sum # should not be called at all
      end
    end
    if randomSeq.length != seqLength
      puts "The length of the randomly generated sequence is wrong"
    else
      randomSeq
    end
  end

  def evolveSeq(startSeq,expectedNumberOfSubstitutionsPerSite)
    # Calculate the substitution probability matrix pMatrix = exp(t*qMatrix), with t = expectedNumberOfSubstitutionsPerSite. See Lio & Goldman 1998 and http://en.wikipedia.org/wiki/Substitution_model.
    # Also check the Ruby documentation for class Matrix at http://www.ruby-doc.org/stdlib-1.9.3/libdoc/matrix/rdoc/Matrix.html.
    actualNumberOfSubstitutions = 0
    endSeq = ""
    
    # Unless the RSRuby instance had been created earlier, do so now.
    @r = RSRuby.instance if @r == nil

    # If there's no rate variation among sites, we just need one substitution matrix for all bases.
    if @alpha == 0
      t = expectedNumberOfSubstitutionsPerSite
      elt = []
      @lMatrix.row_size.times do |i|
        elt[i] = []
        @lMatrix.row_size.times do |j|
          if i == j
            elt[i][j] = Math::E**(@lMatrix.component(i,j)*t)
          else
            elt[i][j] = 0
          end
        end
      end
      eltMatrix = Matrix.rows(elt)
      pMatrix = @u_infMatrix * eltMatrix * @uMatrix
      seqBase = ""
      substitutionProbabilities = []
      startSeq.length.times do |n|
        seqBase = startSeq[n]
        substitutionProbabilities = pMatrix.row(@bases.index(seqBase)).to_a
        random = rand
        if random <= substitutionProbabilities[0]
          replaceBase = @bases[0]
          endSeq << replaceBase
          actualNumberOfSubstitutions += 1 unless replaceBase == seqBase
        else
          (substitutionProbabilities.size-1).times do |j|
            if random > substitutionProbabilities[0..j].sum and random <= substitutionProbabilities[0..(j+1)].sum
              replaceBase = @bases[j+1]
              endSeq << replaceBase
              actualNumberOfSubstitutions += 1 unless replaceBase == seqBase
              break
            end
          end
        end
      end

    # If there is rate variation among sites, we need to make new substitution matrices for each base.
    else
      startSeq.length.times do |n|
        t = @siteSpecificRateMultiplier[n.modulo(@numberOfRateCategories)] * expectedNumberOfSubstitutionsPerSite
        elt = []
        @lMatrix.row_size.times do |i|
          elt[i] = []
          @lMatrix.row_size.times do |j|
            if i == j
              elt[i][j] = Math::E**(@lMatrix.component(i,j)*t)
            else
              elt[i][j] = 0
            end
          end
        end
        eltMatrix = Matrix.rows(elt)
        pMatrix = @u_infMatrix * eltMatrix * @uMatrix
        seqBase = ""
        substitutionProbabilities = []
        seqBase = startSeq[n]
        substitutionProbabilities = pMatrix.row(@bases.index(seqBase)).to_a
        random = rand
        if random <= substitutionProbabilities[0]
          replaceBase = @bases[0]
          endSeq << replaceBase
          actualNumberOfSubstitutions += 1 unless replaceBase == seqBase
        else
          (substitutionProbabilities.size-1).times do |j|
            if random > substitutionProbabilities[0..j].sum and random <= substitutionProbabilities[0..(j+1)].sum
              replaceBase = @bases[j+1]
              endSeq << replaceBase
              actualNumberOfSubstitutions += 1 unless replaceBase == seqBase
              break
            end
          end
        end
      end # startSeq.length.times do |n|
    end
    return endSeq,actualNumberOfSubstitutions
  end

end
