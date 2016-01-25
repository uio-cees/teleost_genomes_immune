# Michael Matschiner, 2015-09-17

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Get command line arguments.
input_tree_file_name = ARGV[0]
branches_net_div_file_name = ARGV[1]
traits_file_name = ARGV[2]
output_tree_file_name = ARGV[3]
output_table_file_name = ARGV[4]

# Make the output directories if it doesn't exist yet.
FileUtils.mkdir_p(File.dirname(output_tree_file_name)) if output_tree_file_name.include?("/")
FileUtils.mkdir_p(File.dirname(output_table_file_name)) if output_table_file_name.include?("/")

# Read the input tree file and save it as the output tree.
input_tree_file = File.open(input_tree_file_name)
input_tree_lines = input_tree_file.read
output_tree_file = File.open(output_tree_file_name,"w")
output_tree_file.write(input_tree_lines)

# Read the file with trait data (the input file used for slouch).
traits_file = File.open(traits_file_name)
traits_lines = traits_file.readlines
trait_file_species_ids = []
trait_file_trait_values = []
traits_lines[1..-1].each do |l|
	trait_file_species_ids << l.split[1]
	trait_file_trait_values << l.split[5].to_f
end

# Read the file with net diversification estimates per branch.
branches_net_div_file = File.open(branches_net_div_file_name)
branches_net_div_lines = branches_net_div_file.readlines
branch_begins = []
branch_ends = []
branch_net_divs = []
branch_tip_labels = []
tip_labels = []
in_branch_begins = false
in_branch_ends = false
in_branch_net_divs = false
in_branch_tip_labels = false
branches_net_div_lines.each do |l|
	if l.strip == "Branch begin"
		in_branch_begins = true
	elsif l.strip == "Branch end"
		in_branch_ends = true
		in_branch_begins = false
	elsif l.strip == "Branch rate"
		in_branch_net_divs = true
		in_branch_ends = false
	elsif l.strip == "Branch tip label"
		in_branch_tip_labels = true
		in_branch_net_divs = false
	elsif in_branch_begins
		branch_begins << l.strip.to_f unless l.strip == ""
	elsif in_branch_ends
		branch_ends << l.strip.to_f unless l.strip == ""
	elsif in_branch_net_divs
		branch_net_divs << l.strip.to_f unless l.strip == ""
	elsif in_branch_tip_labels
		tip_labels << l.strip unless l.strip == ""
	end
end

# Determine maximum branch age.
branch_max_age = 0
branch_ends.each do |e|
	branch_max_age = e if e > branch_max_age
end

# Fill array branch_tip_labels according to branch_begins, branch_ends, and branch_net_divs.
branch_ends.size.times do |x|
	if (branch_max_age-branch_ends[x]).round(2) == 0
		branch_tip_labels << tip_labels.shift
	else
		branch_tip_labels << "unknown"
	end
end

# Prepare the output table string.
table_string = "Species\tlog_MHCI_copy_number\tnet_diversification\n"
trait_file_species_ids.size.times do |x|
	unless trait_file_species_ids[x] == "NA"
		species_id = trait_file_species_ids[x]
		trait_value = trait_file_trait_values[x]
		net_diversification = branch_net_divs[branch_tip_labels.index(species_id)]
		table_string << "#{species_id}\t#{trait_value}\t#{net_diversification}\n"
	end
end

# Write the output table.
output_table_file = File.open(output_table_file_name,"w")
output_table_file.write(table_string)
output_table_file.close
