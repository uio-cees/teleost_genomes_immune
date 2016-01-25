# Michael Matschiner, 2015-06-27

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Define class Branch.
class Branch
	attr_reader :from, :to, :length
	def initialize(from, to, length)
		@from = from
		@to = to
		@length = length
		@consistency_indeces = []
	end
	def add_consistency_index(consistency_index)
		@consistency_indeces << consistency_index
	end
	def to_s
		"#{from}\t#{to}\t#{length}"
	end
	def number_of_consistent_indels
		number_of_consistent_indels = 0
		@consistency_indeces.each {|i| number_of_consistent_indels += 1 if i == 1.0}
		number_of_consistent_indels
	end
	def proportion_of_consistent_indels
		if @consistency_indeces.size == 0
			"NA"
		else
			self.number_of_consistent_indels/@consistency_indeces.size.to_f
		end
	end
end

# Get command line arguments.
paup_output_file_name = ARGV[0]
branch_table_file_name = ARGV[1]

# Make the branch table directory if it does not exist yet.
if branch_table_file_name.include?("/")
	branch_table_dir = branch_table_file_name.chomp(branch_table_file_name.split("/").last).chomp("/")
	# Create the output directories if they don't exist yet.
	FileUtils.mkdir_p(branch_table_dir)
end

# Parse the PAUP* output file and get the three tables.
paup_output_file = File.open(paup_output_file_name)
paup_output_lines = paup_output_file.readlines
just_before_first_table = false
in_first_table = false
after_first_table = false
first_table_lines = []
just_before_second_table = false
in_second_table = false
after_second_table = false
second_table_lines = []
just_before_third_table = false
in_third_table = false
after_third_table = false
third_table_lines = []
paup_output_lines.each do |l|
	if after_first_table == false and l.strip == "Branch lengths and linkages for tree #1"
		just_before_first_table = true
	elsif after_first_table and l.strip == "Branch lengths and linkages for tree #1"
		just_before_second_table = true
	elsif after_first_table == false and l.strip[0..4] == "Sum  "
		after_first_table = true
		just_before_first_table = false
	elsif after_first_table and l.strip[0..4] == "Sum  "
		after_second_table = true
		just_before_second_table = false
	elsif l.strip == "Apomorphy lists:"
		just_before_third_table = true
	end
	if just_before_first_table
		if in_first_table == false and l.strip.match(/^-+$/)
			in_first_table = true
		elsif in_first_table == true and l.strip.match(/^-+$/)
			in_first_table = false
		elsif in_first_table
			first_table_lines << l
		end
	end
	if just_before_second_table
		if in_second_table == false and l.strip.match(/^-+$/)
			in_second_table = true
		elsif in_second_table == true and l.strip.match(/^-+$/)
			in_second_table = false
		elsif in_second_table
			second_table_lines << l
		end
	end
	if just_before_third_table
		if in_third_table == false and l.strip.match(/^-+$/)
			in_third_table = true
		elsif in_third_table == true and l.strip.match(/^-+$/)
			in_third_table = false
		elsif in_third_table
			third_table_lines << l
		end
	end
end

# Extract branch info from the first two tables.
branch = []
first_table_lines.size.times do |x|
	from = second_table_lines[x+1][20..39].sub(/\(.+\)/,"").sub("*","").strip.sub(" ","_")
	from = "node_" + from if from.match(/^\d+$/)
	to = second_table_lines[x+1][0..19].sub(/\(.+\)/,"").sub("*","").strip.sub(" ","_")
	to = "node_" + to if to.match(/^\d+$/)
	length = first_table_lines[x][30..-1].strip.to_f
	branch << Branch.new(from, to, length)
end

# For each line in the third table, identify the respective branch.
current_branch = nil
third_table_lines.each do |l|
	if l[0..22].include?("-->")
		node_ary = l[0..22].split("-->")
		from = node_ary[0].strip.sub(" ","_")
		to = node_ary[1].strip.sub(" ","_")
		branch.each do |b|
			if b.from == from and b.to == to
				current_branch = b
				break
			end
		end
	end
	if current_branch == nil
		raise_string = "ERROR: The following table line could not be assigned to branch!\n"
		raise_string << l
		raise raise_string
	end
	consistency_index = l[23..-1].split[2].to_f
	current_branch.add_consistency_index(consistency_index)
end

# Prepare an output table.
branch_table_string = "from\tto\tlength\tlog_length\tproportion_consistent\n"
branch.each do |b|
	branch_table_string << "#{b.from}\t#{b.to}\t#{b.length}\t#{Math.log(b.length)}\t#{b.proportion_of_consistent_indels}\n"
end

# Write the output table.
branch_table_file_name = 
branch_table_file = File.open(branch_table_file_name,"w")
branch_table_file.write(branch_table_string)
branch_table_file.close
