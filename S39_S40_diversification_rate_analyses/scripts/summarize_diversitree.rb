# Michael Matschiner, 2015-09-25

# Load script array_stats for array statistics.
require "./resources/array_stats.rb"

# Load FileUtils.
require 'fileutils'

class Analysis
	attr_reader :threshold, :replicate, :results_string, :lnLik_constrained, :lnLik_free
	attr_reader :parameters_constrained, :parameters_free
	def initialize(threshold,replicate,results_string)
		@threshold = threshold
		@replicate = replicate
		@results_string = results_string
		@results_lines = @results_string.split("\n")
		counts_lines = []
		@results_lines.size.times {|x| counts_lines << @results_lines[x+1] if @results_lines[x].strip == "$counts"}
		@parameters_constrained = @results_lines[0..5]
		@lnLik_constrained = @results_lines[6].to_f
		@parameters_free = @results_lines[7..12]
		@lnLik_free = @results_lines[13].to_f
	end
	def lambda_constrained
		@parameters_constrained[0].to_f
	end
	def lambda0_free
		@parameters_free[0].to_f
	end
	def lambda1_free
		@parameters_free[1].to_f
	end
	def mu_constrained
		@parameters_constrained[2].to_f
	end
	def mu_free
		@parameters_free[2].to_f
	end
	def net_div_constrained
		lambda_constrained - mu_constrained
	end
	def net_div0_free
		lambda0_free - mu_free
	end
	def net_div1_free
		lambda1_free - mu_free
	end
	def q01_constrained
		@parameters_constrained[4].to_f
	end
	def q10_constrained
		@parameters_constrained[5].to_f
	end
	def q01_free
		@parameters_free[4].to_f
	end
	def q10_free
		@parameters_free[5].to_f
	end
end

# Get command line arguments.
diversitree_dir = ARGV[0]
summary_dir = ARGV[1]

# Make the summary directory if it doesn't exist yet.
FileUtils.mkdir_p(summary_dir)

# Get analysis directory names.
analyses = []
thresholds = []
Dir.entries(diversitree_dir).each do |e|
	if e.match(/t\d\d_r\d\d\d/)
		if File.exists?("#{diversitree_dir}/#{e}/results.txt")
			threshold = e[1..2].to_i
			if threshold >= 10
				thresholds << threshold unless thresholds.include?(threshold)
				replicate = e[5..7].to_i
				results_file = File.open("#{diversitree_dir}/#{e}/results.txt")
				results_string = results_file.read
				analyses << Analysis.new(threshold,replicate,results_string)
			end
		end
	end
end

# Sort all thresholds.
sorted_thresholds = thresholds.sort

# Prepare a table for threshold, mean and stdev of constrained net_div, unconstrained lower net_div, unconstrained upper net_div, 
# mean and stdev of constrained likelihood and of unconstrained likelihood, and of the difference between the two.
outstring = "threshold\t"
outstring << "mean_constrained_netdiv\t"
outstring << "stdev_constrained_netdiv\t"
outstring << "mean_unconstrained_netdiv_low_copy_number\t"
outstring << "stdev_unconstrained_netdiv_low_copy_number\t"
outstring << "mean_unconstrained_netdiv_high_copy_number\t"
outstring << "stdev_unconstrained_netdiv_high_copy_number\t"
outstring << "median_constrained_log_likelihood\t"
outstring << "lower_constrained_log_likelihood\t"
outstring << "upper_constrained_log_likelihood\t"
outstring << "median_unconstrained_log_likelihood\t"
outstring << "lower_unconstrained_log_likelihood\t"
outstring << "upper_unconstrained_log_likelihood\t"
outstring << "median_log_likelihood_difference\t"
outstring << "lower_log_likelihood_difference\t"
outstring << "upper_log_likelihood_difference\n"
sorted_thresholds.each do |t|
	constrained_net_divs = []
	unconstrained_net_divs_low_copy_number = []
	unconstrained_net_divs_high_copy_number = []
	constrained_log_likelihoods = []
	unconstrained_log_likelihoods = []
	log_likelihood_differences = []
	analyses.each do |a|
		if a.threshold == t
			constrained_net_divs << a.net_div_constrained
			unconstrained_net_divs_low_copy_number << a.net_div0_free
			unconstrained_net_divs_high_copy_number << a.net_div1_free
			constrained_log_likelihoods << a.lnLik_constrained
			unconstrained_log_likelihoods << a.lnLik_free
			log_likelihood_differences << a.lnLik_free - a.lnLik_constrained
		end
	end
	outstring << "#{t}\t"
	outstring << "#{constrained_net_divs.mean}\t"
	outstring << "#{constrained_net_divs.standard_deviation}\t"
	outstring << "#{unconstrained_net_divs_low_copy_number.mean}\t"
	outstring << "#{unconstrained_net_divs_low_copy_number.standard_deviation}\t"
	outstring << "#{unconstrained_net_divs_high_copy_number.mean}\t"
	outstring << "#{unconstrained_net_divs_high_copy_number.standard_deviation}\t"
	outstring << "#{constrained_log_likelihoods.median}\t"
	outstring << "#{constrained_log_likelihoods.rough_percentile(0.05)}\t"
	outstring << "#{constrained_log_likelihoods.rough_percentile(0.95)}\t"
	outstring << "#{unconstrained_log_likelihoods.median}\t"
	outstring << "#{unconstrained_log_likelihoods.rough_percentile(0.05)}\t"
	outstring << "#{unconstrained_log_likelihoods.rough_percentile(0.95)}\t"
	outstring << "#{log_likelihood_differences.median}\t"
	outstring << "#{log_likelihood_differences.rough_percentile(0.05)}\t"
	outstring << "#{log_likelihood_differences.rough_percentile(0.95)}\n"
end

output_file_name = "summary.txt"
output_file = File.open("#{summary_dir}/#{output_file_name}","w")
output_file.write(outstring)
