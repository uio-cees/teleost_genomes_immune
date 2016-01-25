module Enumerable
    def sum
      self.inject(0){|accum, i| accum + i }
    end
    def mean
      self.sum/self.length.to_f
    end
    def sample_variance
      m = self.mean
      sum = self.inject(0){|accum, i| accum +(i-m)**2 }
      sum/(self.length - 1).to_f
    end
    def standard_deviation
      return Math.sqrt(self.sample_variance)
    end
end

Fixnum.class_eval { def even?; self%2 == 0; end }
Fixnum.class_eval { def odd?; self%2 != 0; end }
