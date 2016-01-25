class ECM
  def initialize(fileName)
    lines = IO.readlines(fileName)
    sLines = lines[0..59]
    piLines = lines[61]
    codonLines = lines[64..67]
    @s = []
    @s[0] = [nil]
    sLines.size.times do |i|
      @s[i+1] = sLines[i].split(" ")
      @s[i+1] << nil
    end
    @s.size.times do |i|
      @s.size.times do |j|
        @s[i][j] = @s[j][i] unless @s[j][i] == nil
      end
    end
    # convert all @s from string to float
    @s.size.times do |i|
      @s.size.times do |j|
        @s[i][j] = @s[i][j].to_f
      end
    end
    @sMatrix = Matrix.rows(@s)
    
    @pi = piLines.split(" ")
    # convert all @pi from string to float
    @pi.size.times do |i|
      @pi[i] = @pi[i].to_f
    end
    # standardize all @pi (this only corrects for rounding errors)
    sum = @pi.sum
    @pi.size.times do |i|
      @pi[i] = @pi[i]/sum
    end
    
    # calculate the Q matrix as q_ij = pi_j*s_ij
    @q = []
    @s.size.times do |i|
      @q[i] = []
      @s.size.times do |j|
        if j == i
          @q[i][j] = 0
        else
          @q[i][j] = @pi[j]*@s[i][j]
        end
      end
      @q[i][i] = -@q[i].sum
    end
    @qMatrix = Matrix.rows(@q)
    @u_infMatrix, @lMatrix, @uMatrix = @qMatrix.eigensystem

    @codon = []
    codonLines.size.times do |c|
      @codon << codonLines[c].split(" ")
    end
    @codon.flatten!
    @codon.size.times do |c|
      @codon[c] = @codon[c].downcase
    end
  end
  
  def randomSeq(seqLength)
    numberOfCodons = seqLength/3
    puts "WARNING: The chosen sequence length is not a multiple of 3. It will be shortened by #{seqLength-numberOfCodons*3} characters to a length of #{numberOfCodons*3} bp." if numberOfCodons*3 != seqLength
    randomSeq = ""
    numberOfCodons.times do |n|
      random = rand
      if random <= @pi[0]
        randomSeq << @codon[0]
      else
        (@pi.size-1).times do |i|
          if random > @pi[0..i].sum and random <= @pi[0..(i+1)].sum
            randomSeq << @codon[i+1]
            break
          end
        end
        puts "Something wrong about the codon frequencies" if random > @pi.sum # should not be called at all
      end
    end
    if randomSeq.length != numberOfCodons*3
      puts "The length of the randomly generated sequence is wrong"
    else
      randomSeq
    end
  end

  def evolveSeq(startSeq, expectedNumberOfSubstitutionsPerSite)
    # calculate the substitution probability matrix pMatrix = exp(t*qMatrix), with t = expectedNumberOfSubstitutionsPerSite. See Lio & Goldman 1998 and http://en.wikipedia.org/wiki/Substitution_model.
    # Also check the Ruby documentation for class Matrix at http://www.ruby-doc.org/stdlib-1.9.3/libdoc/matrix/rdoc/Matrix.html.
    numberOfCodons = startSeq.length/3
    puts "The start sequence length is not a multiple of 3. It will be shortened by #{startSeq.length-numberOfCodons*3} characters to a length of #{numberOfCodons*3} bp." if numberOfCodons*3 != startSeq.length
    endSeq = ""
    actualNumberOfSubstitutions = 0

    # If there's no rate variation among sites, we just need one substitution matrix for all codons.
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
    seqCodon = ""
    substitutionProbabilities = []
    numberOfCodons.times do |n|
      seqCodon = startSeq[(3*n)..(3*n+2)]
      substitutionProbabilities = pMatrix.row(@codon.index(seqCodon)).to_a
      random = rand
      if random <= substitutionProbabilities[0]
        replaceCodon = @codon[0]
        actualNumberOfSubstitutions += 1 unless replaceCodon[0] == seqCodon[0]
        actualNumberOfSubstitutions += 1 unless replaceCodon[1] == seqCodon[1]
        actualNumberOfSubstitutions += 1 unless replaceCodon[2] == seqCodon[2]
        endSeq << replaceCodon
      else
        (substitutionProbabilities.size-1).times do |j|
          if random > substitutionProbabilities[0..j].sum and random <= substitutionProbabilities[0..(j+1)].sum
            replaceCodon = @codon[j+1]
            endSeq << replaceCodon
            actualNumberOfSubstitutions += 1 unless replaceCodon[0] == seqCodon[0]
            actualNumberOfSubstitutions += 1 unless replaceCodon[1] == seqCodon[1]
            actualNumberOfSubstitutions += 1 unless replaceCodon[2] == seqCodon[2]
            break
          end
        end
      end
    end
    return endSeq,actualNumberOfSubstitutions
  end

end