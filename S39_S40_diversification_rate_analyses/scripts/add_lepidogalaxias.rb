# Michael Matschiner, 2015-07-01.

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Get command line arguments.
input_tree_file_name = ARGV[0]
output_tree_file_name = ARGV[1]

# Make the output tree directory if it doesn't exist yet.
if output_tree_file_name.include?("/")
	output_tree_directory = output_tree_file_name.chomp(output_tree_file_name.split("/")[-1])
	FileUtils.mkdir_p(output_tree_directory)
end

# Read the input tree file.
input_tree_file = File.open(input_tree_file_name)
input_tree_lines = input_tree_file.readlines

# Extract tree line.
newick_string = ""
input_tree_lines.each do |l|
	if l.include?("=") and l.downcase.include?("tree")
		newick_string = l.split("=")[1].gsub(/\[.+\]/,"").strip
	end
end

# Get the last branch length.
raise "ERROR: Last branch length could not be read!" unless newick_string.match(/:(\d+\.\d+)\);$/)
last_branch_length = $1.to_f

# Split this branch length randomly.
first_partial_branch_length = rand * (last_branch_length)
second_partial_branch_length = last_branch_length - first_partial_branch_length

# Get the branch length of Otomorpha (= the root height).
unless newick_string.match(/\(\(Astmex:(\d+\.\d+),Danrer:(\d+\.\d+)\):(\d+\.\d+),/)
	puts "ERROR: Unexpected tree string:"
	puts newick_string
	raise
end
raise "ERROR: Branch lengths differ for Astmex and Danrer!" unless $1 == $2
root_height = $1.to_f + $3.to_f

# Add Lepsal to the newick string.
newick_string.sub!(/\(\(Astmex:(\d+\.\d+),Danrer:(\d+\.\d+)\):(\d+\.\d+),/,"((Astmex:#{$1},Danrer:#{$2}):#{$3},(Lepsal:#{root_height-second_partial_branch_length},")
newick_string.sub!(/:\d+\.\d+\);$/,":#{first_partial_branch_length}):#{second_partial_branch_length});")

# Write the new newick string to file.
output_tree_file = File.open(output_tree_file_name,"w")
output_tree_file.write(newick_string)
output_tree_file.close
