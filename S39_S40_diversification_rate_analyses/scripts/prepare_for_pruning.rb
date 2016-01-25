# Michael Matschiner, 2015-07-28

# Get the command line arguments.
input_table_file_name = ARGV[0]
trait_table_file_name = ARGV[1]
output_table_file_name = ARGV[2]
input_tree_file_name = ARGV[3]
output_tree_file_name = ARGV[4]

# Read the table file.
input_table_file = File.open(input_table_file_name,"r")
input_table_lines = input_table_file.readlines
species = []
clades = []
input_table_lines[1..-1].each do |l|
	spc = l.split[0].strip
	species << spc
	clade = l.split[1].strip
	# clade = "Gadiformes1" if spc == "fish_97"
	# clade = "Gadiformes1" if spc == "fish_7"
	# clade = "Gadiformes1" if spc == "fish_2"
	# clade = "Gadiformes1" if spc == "fish_1"
	# clade = "Gadiformes1" if spc == "fish_6"
	# clade = "Gadiformes1" if spc == "fish_5"
	# clade = "Gadiformes1" if spc == "fish_4"
	# clade = "Gadiformes1" if spc == "fish_8"
	# clade = "Gadiformes1" if spc == "fish_3"
	# clade = "Gadiformes1" if spc == "fish_12"
	# clade = "Gadiformes1" if spc == "fish_10"
	# clade = "Gadiformes1" if spc == "fish_11"
	# clade = "Gadiformes1" if spc == "fish_9"
	# clade = "Gadiformes1" if spc == "fish_95"
	# clade = "Gadiformes2" if spc == "fish_18"
	# clade = "Gadiformes2" if spc == "fish_17"
	# clade = "Gadiformes2" if spc == "fish_19"
	# clade = "Gadiformes2" if spc == "fish_23"
	# clade = "Gadiformes2" if spc == "fish_22"
	# clade = "Gadiformes2" if spc == "fish_94"
	# clade = "Gadiformes2" if spc == "fish_36"
	# clade = "Gadiformes2" if spc == "fish_20"
	# clade = "Gadiformes2" if spc == "fish_16"
	# clade = "Gadiformes2" if spc == "fish_14"
	# clade = "Gadiformes2" if spc == "fish_13"
	# clade = "Gadiformes2" if spc == "fish_15"
	clades << clade
end

# Read the trait table file.
trait_table_file = File.open(trait_table_file_name,"r")
trait_table_lines = trait_table_file.readlines
trait_species = []
trait_table_lines.each do |l|
	trait_species << l.strip unless l.strip == ""
end

# Reduce the arrays for species and clades to include only species that also appear in the trait data table.
red_species = []
red_clades = []
number_of_duplicates = 0
species.size.times do |x|
	if trait_species.include?(species[x])
		red_species << species[x]
		red_clades << clades[x]
	else
		red_species << species[x]
		number_of_duplicates += 1
		red_clades << "prune_#{number_of_duplicates}"
	end
end

# Rename clades to remove clade duplicates.
unique_clades = []
red_clades.size.times do |x|
	unless red_clades[x].match(/prune_/)
		if unique_clades.include?(red_clades[x])
			number_of_duplicates += 1
			red_clades[x] = "prune_#{number_of_duplicates}"
		else
			unique_clades << red_clades[x]
		end
	end
end

# Read the input tree file.
input_tree_file = File.open(input_tree_file_name,"r")
input_tree_string = input_tree_file.read

# Replace species ids in the tree.
red_species.size.times do |x|
	input_tree_string.sub!(/#{red_species[x]}\n/,"#{red_clades[x]}\n")
	input_tree_string.sub!(/#{red_species[x]}:/,"#{red_clades[x]}:")
	input_tree_string.sub!(/#{red_species[x]}\)/,"#{red_clades[x]})")
end

# Write the edited tree string to file.
output_tree_file = File.open(output_tree_file_name,"w")
output_tree_file.write(input_tree_string)
output_tree_file.close

# Prepare the table of taxa to keep.
output_table_string = ""
unique_clades.each do |c|
	output_table_string << "#{c}\n"
end

# Write the table of taxa to keep to file.
output_table_file = File.open(output_table_file_name,"w")
output_table_file.write(output_table_string)
output_table_file.close
