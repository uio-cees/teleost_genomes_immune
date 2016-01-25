# Michael Matschiner, 2015-07-28

# Load FileUtils.
require 'fileutils'

# Get the command line arguments.
input_tree_file_name = ARGV[0]
clade_table_file_name = ARGV[1]
species_file_name = ARGV[2]
slouch_input_file_name = ARGV[3]
output_file_name_w_relative_path = ARGV[4]

# Make the output directory if it doesn't exist yet.
if output_file_name_w_relative_path.include?("/")
	output_dir = output_file_name_w_relative_path.chomp(output_file_name_w_relative_path.split("/").last)
	FileUtils.mkdir_p(output_dir)
end

# Read the tree file to get all clade names.
input_tree_file = File.open(input_tree_file_name,"r")
input_tree_string = input_tree_file.read
clade_names_raw = input_tree_string.scan(/[,\(].*?:/)
clade_names = []
clade_names_raw.each do |c|
	clade_name = c.gsub("(","").sub(",","").sub(":","")
	clade_names << clade_name
end

# Read the clade table file.
clade_table_file = File.open(clade_table_file_name)
clade_table_lines = clade_table_file.readlines
clade_table_file.close
clades = []
clade_species_numbers = []
clade_table_lines.each do |l|
	if l[0] != "#" and l.strip != ""
		if l[0] != "\t"
			clades << l.strip
			clade_species_numbers << 0
		else
			clade_species_numbers[-1] += l.split("\t")[-1].to_i
		end
	end
end

# Read the slouch input file.
slouch_input_file = File.open(slouch_input_file_name)
slouch_input_lines = slouch_input_file.readlines
slouch_species_names = []
slouch_trait_values = []
slouch_input_lines[1..-1].each do |l|
	line_ary = l.split
	unless line_ary[1] == "NA"
		slouch_species_names << line_ary[1]
		slouch_trait_values << Math.exp(line_ary[5].to_f)
	end
end

# Read the file with species-clade assignments.
species_file = File.open(species_file_name)
species_lines = species_file.readlines
species_file_species = []
species_file_clades = []
species_lines.each do |l|
	species_file_species << l.split[0]
	species_file_clades << l.split[1]
end

# Prepare the output string.
outstring = "tip.label\tnc\ttrait_values\n"
clade_names.each do |c|
	outstring << "#{c}\t#{clade_species_numbers[clades.index(c)]}"
	species_file_clades.size.times do |x|
		if species_file_clades[x] == c
			slouch_species_names.size.times do |y|
				if slouch_species_names[y] == species_file_species[x]
					outstring << "\t#{slouch_trait_values[y]}"
				end
			end
		end
	end
	outstring << "\n"
end

# Write the output string to file.
output_file = File.open(output_file_name_w_relative_path,"w")
output_file.write(outstring)
output_file.close
