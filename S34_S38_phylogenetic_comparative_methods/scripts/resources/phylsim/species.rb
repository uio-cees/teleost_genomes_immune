class Species
  attr_reader :id, :branch, :origin, :termination, :duration, :extant, :color
  def initialize(id, branch)
    @id = id
    @branch = [branch]
    @origin = branch.origin
    @termination = branch.termination
    @duration = @origin - @termination
    @color = rand
  end
  def addBranch(branch)
    raise "Branch to be added seems not to be connected to existing branch of species #{id}!" unless @termination == branch.origin
    @branch << branch
    @termination = branch.termination
    @duration = @origin - @termination
  end
  def addExtant(present)
    if @termination == present
      @extant = true
    else
      @extant = false
    end
  end
end
