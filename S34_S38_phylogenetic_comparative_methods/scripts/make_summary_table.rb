# Michael Matschiner, 2015-07-21

# Get command line argument.
slouch_output_dir = ARGV[0]
table_file_name = ARGV[1]

# Get the Slouch result file names.
result_file_names =[]
Dir.entries(slouch_output_dir).each {|e| result_file_names << e if e.match(/\.results\.txt/)}

# Prepare arrays.
regime_codes = []
alphas = []
sigsquares = []
thetas = []
lnliks = []
aics = []
aiccs = []
rsquares = []

# Read all result files.
result_file_names.each do |r|
	regime_code = r.chomp(".txt").chomp(".results")
	regime_codes << regime_code
	regime_ids = ["1"]
	regime_code.size.times do |x|
		regime_ids << x+2 if regime_code[x] == "1"
	end
	result_file = File.open("#{slouch_output_dir}/#{r}")
	result_lines = result_file.readlines
	in_theta_part = false
	theta_values = []
	result_lines.each do |l|
		alphas << l.split.last.to_f if l.include?("Rate of adaptation")
		sigsquares << l.split.last.to_f if l.include?("Stationary variance")
		lnliks << l.split.last.to_f if l.include?("Support")
		aics << l.split.last.to_f if l.include?("AIC ")
		aiccs << l.split.last.to_f if l.include?("AICc")
		rsquares << l.split.last.to_f if l.include?("r squared")
		if l.include?("PRIMARY OPTIMA")
			in_theta_part = true
		elsif l.include?("------------------")
			in_theta_part = false
		elsif in_theta_part and l.strip != "" and l.include?("mapped") == false and l.include?("Estimates") == false
			theta_values << l.split[-2].to_f
		end
	end
	if theta_values.size == regime_ids.size
		theta = [theta_values.shift]
		regime_code.size.times do |x|
			if regime_code[x] == "1"
				theta << theta_values.shift
			else
				theta << "NA"
			end
		end
		thetas << theta
	else
		warn "WARNING: Theta values could not be read for file #{r}!"
	end
end

# Calculate AIC weights.
best_aicc = aiccs.min
relative_likelihoods = []
aiccs.each do |a|
	relative_likelihoods << Math.exp(-0.5 * (a-best_aicc))
end
relative_likelihoods_sum = 0
relative_likelihoods.each {|a| relative_likelihoods_sum += a}
aicc_weights = []
relative_likelihoods.each {|a| aicc_weights << a/relative_likelihoods_sum}

# Prepare a table string.
table_string = "analysis\talpha\tsigsquare\ttheta0\ttheta1\ttheta2\ttheta3\ttheta4\ttheta5\ttheta6\ttheta7\ttheta8\trsquared\tlnlik\taic\taicc\taicc_weight\n"
alphas.size.times do |x|
	table_string << "#{regime_codes[x]}\t"
	table_string << "#{format("%.3f", alphas[x])}\t"
	table_string << "#{format("%.3f", sigsquares[x])}\t"
	if thetas[x][0] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][0])}\t"
	end
	if thetas[x][1] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][1])}\t"
	end
	if thetas[x][2] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][2])}\t"
	end
	if thetas[x][3] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][3])}\t"
	end
	if thetas[x][4] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][4])}\t"
	end
	if thetas[x][5] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][5])}\t"
	end
	if thetas[x][6] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][6])}\t"
	end
	if thetas[x][7] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][7])}\t"
	end
	if thetas[x][8] == "NA"
		table_string << "NA\t"
	else
		table_string << "#{format("%.3f", thetas[x][8])}\t"
	end
	table_string << "#{format("%.3f", rsquares[x])}\t"
	table_string << "#{format("%.3f", lnliks[x])}\t"
	table_string << "#{format("%.3f", aics[x])}\t"
	table_string << "#{format("%.3f", aiccs[x])}\t"
	table_string << "#{format("%.3f", aicc_weights[x])}\n"
end

# Write the table to file.
table_file = File.open(table_file_name,"w")
table_file.write(table_string)
table_file.close
