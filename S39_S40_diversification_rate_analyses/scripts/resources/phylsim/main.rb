require 'rubygems'
require 'rubystats'
require 'bio'
require 'matrix'
require 'rsruby'

require "#{$libPath}util.rb"
require "#{$libPath}ecm.rb"
require "#{$libPath}gtr.rb"
require "#{$libPath}branch.rb"
require "#{$libPath}species.rb"
require "#{$libPath}fossil.rb"
require "#{$libPath}tree.rb"

class RetryException < RuntimeError
  attr :okToRetry
  def initialize(okToRetry)
    @okToRetry = okToRetry
  end
end

