class Branch
  attr_reader :id, :origin, :termination, :parentId, :daughterId, :progenyId, :speciesId, :endCause, :extant, :duration, :fossil, :progeniesComplete, :progenyPassedOn, :extantProgenyId, :originalExtantProgenyId, :extantDiversity, :originalExtantDiversity, :startRate, :endRate, :meanRate, :rate, :startSeq, :endSeq, :startTrait, :endTrait, :focusGroup, :actualNumberOfSubstitutions, :actualNumberOfSubstitutionsPerSite, :alignmentPosition
  attr_writer :progeniesComplete, :progenyPassedOn, :progenyId, :extantProgenyId
  def initialize(id, origin, termination, parentId, daughterId, endCause)
    @id = id
    @origin = origin
    @termination = termination
    @originBeforeErrorIntroduction = nil
    @terminationBeforeErrorIntroduction = nil
    @parentId = parentId
    @daughterId = daughterId # an array with two items
    @progenyId = []
    @extantProgenyId = []
    @originalExtantProgenyId = []
    @extantDiversity = 0
    @originalExtantDiversity = 0
    @endCause = endCause
    @extant = false
    @duration = @origin - @termination
    @progeniesComplete = 0
    @progenyPassedOn = false
    @speciesId = nil
    @rate = nil
    @fossil = []
  end
  def cutAtPresent(present)
    if @termination <= present
      @termination = present
      @duration = @origin-@termination
      @endCause = "present"
      @extant = true
      @daughterId = ["none","none"]
    else
      raise "Branch could not be cut at present, cause it terminates before the present!"
    end
  end
  def originBeforeErrorIntroduction
    if @originBeforeErrorIntroduction == nil
      @origin
    else
      @originBeforeErrorIntroduction
    end
  end
  def terminationBeforeErrorIntroduction
    if @terminationBeforeErrorIntroduction == nil
      @termination
    else
      @terminationBeforeErrorIntroduction
    end
  end
  def updateDaughterId(daughterId)
    @daughterId = daughterId
  end
  def updateParentId(parentId)
    @parentId = parentId
  end
  def addFossils(fossil)
    @fossil = fossil # an array with flexible number of items
  end
  def addStartRate(rate)
    @startRate = rate
  end
  def addEndRate(rate)
    @endRate = rate
  end
  def addMeanRate(rate)
    @meanRate = rate
    @rate = rate
  end
  def addStartSeq(seq)
    @startSeq = seq
  end
  def addEndSeq(seq)
    @endSeq = seq
  end
  def addStartTrait(trait)
    @startTrait = trait
  end
  def addEndTrait(trait)
    @endTrait = trait
  end
  def addSpeciesId(speciesId)
    @speciesId = speciesId
  end
  def updateSpeciesId(speciesId)
    @speciesId = speciesId
  end
  def updateEndCause(endCause)
    @endCause = endCause
  end
  def updateTermination(termination)
    @termination = termination
    @duration = @origin-@termination
  end
  def updateExtant(extant)
    @extant = extant
  end
  def updateRate(rate)
    @rate = rate
  end
  def updateDuration(duration)
    @duration = duration
  end
  def updateFossil(fossil)
    @fossil = fossil
  end
  def updateProgenyId(progenyId)
    @progenyId = progenyId
  end
  def updateExtantProgenyId(extantProgenyId)
    @originalExtantProgenyId = @extantProgenyId if @originalExtantProgenyId == []
    @extantProgenyId = extantProgenyId
  end
  def updateExtantDiversity(extantDiversity)
    @extantDiversity = extantDiversity
  end
  def updateOriginalExtantDiversity(originalExtantDiversity)
    @originalExtantDiversity = originalExtantDiversity
  end
  def addFocusGroup(focusGroup)
    @focusGroup = focusGroup
  end
  def addActualNumberOfSubstitutions(actualNumberOfSubstitutions)
    @actualNumberOfSubstitutions = actualNumberOfSubstitutions
  end
  def addActualNumberOfSubstitutionsPerSite(actualNumberOfSubstitutionsPerSite)
    @actualNumberOfSubstitutionsPerSite = actualNumberOfSubstitutionsPerSite
  end
  def addAlignmentPosition(alignmentPosition)
    @alignmentPosition = alignmentPosition
  end
  def updateEndSeq(seq)
    @endSeq = seq
  end
  def updateOriginAfterErrorIntroduction(origin)
    @originBeforeErrorIntroduction = @origin if @originBeforeErrorIntroduction == nil
    @origin = origin
    @duration = @origin-@termination
  end
  def updateTerminationAfterErrorIntroduction(termination)
    @terminationBeforeErrorIntroduction = @termination if @terminationBeforeErrorIntroduction == nil
    @termination = termination
    @duration = @origin-@termination
  end
  def ageError
    if @terminationBeforeErrorIntroduction == nil or @originBeforeErrorIntroduction == nil
      0
    else
      meanAge = (@origin+@termination)/2.0
      meanAgeBeforeErrorIntroduction = (@originBeforeErrorIntroduction+@terminationBeforeErrorIntroduction)/2.0
      if meanAge > meanAgeBeforeErrorIntroduction
        (meanAge-meanAgeBeforeErrorIntroduction)/meanAgeBeforeErrorIntroduction
      else
        (meanAgeBeforeErrorIntroduction-meanAge)/meanAgeBeforeErrorIntroduction
      end
    end
  end
  def expectedNumberOfSubstitutionsPerSite
    if @rate == nil
      nil
    else
      @duration * @rate
    end
  end
  def to_s
    string = ""
    string << "ID:                               #{@id}\n"
    string << "Origin:                           #{@origin}\n"
    string << "Termination:                      #{@termination}\n"
    string << "End cause:                        #{@endCause}\n"
    string << "Extant:                           #{@extant}\n"
    string << "Parent ID:                        #{@parentId}\n"
    string << "Daughter ID 1:                    #{@daughterId[0]}\n"
    string << "Daughter ID 2:                    #{@daughterId[1]}\n"
    string << "Species ID:                       #{@speciesId}\n"
    string << "Extant diversity:                 #{@extantDiversity}\n"
    string << "Original extant diversity:        #{@originalExtantDiversity}\n"
    string << "Focus group:                      #{@focusGroup}\n" unless @focusGroup == nil
    string << "Rate:                             #{@rate}\n" unless @rate == nil
    string << "Number of substitutions:          #{@numberOfSubstitutions}\n" unless @numberOfSubstitutions == nil
    string << "Number of substitutions per site: #{@numberOfSubstitutionsPerSite}\n" unless @numberOfSubstitutionsPerSite == nil
    string << "Alignment position:               #{@alignmentPosition}\n" unless @alignmentPosition == nil
    unless @fossil == []
      first_occurrence_age = 0
      @fossil.size.times {|f| first_occurrence_age = @fossil[f].age.round(3) if first_occurrence_age < @fossil[f].age.round(3)}
      string << "First occurrence age:             #{first_occurrence_age}\n"
    end
    if @traitsEvolved == true
      string << "Trait:                        #{@branch[b].endTrait.round(3)}\n"
    end
    string << "\n"
    string
  end
end
