class Fossil
  attr_reader :branchId, :age
  def initialize(branchId,age)
    @branchId = branchId
    @age = age
  end
end
