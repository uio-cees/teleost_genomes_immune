# Michael Matschiner, 2015-07-01.

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Get the command line arguments.
species_table_file_name = ARGV[0]
clade_table_file_name = ARGV[1]
template_file_name = ARGV[2]
bamm_dir = ARGV[3]
number_of_replicates = ARGV[4].to_i
bamm_executable_name = ARGV[5]

# Read the species table file.
species_table_file = File.open(species_table_file_name)
species_table_lines = species_table_file.readlines
species = []
clades1 = []
species_table_lines[1..-1].each do |l|
	line_ary = l.split
	species << line_ary[0].strip
	clades1 << line_ary[1].strip
end

# Read the clade table file.
clade_table_file = File.open(clade_table_file_name)
clade_table_lines = clade_table_file.readlines
clade_table_file.close
clades2 = []
clade_species_numbers = []
clade_table_lines.each do |l|
	if l[0] != "#" and l.strip != ""
		if l[0] != "\t"
			clades2 << l.strip
			clade_species_numbers << 0
		else
			clade_species_numbers[-1] += l.split("\t")[-1].to_i
		end
	end
end

# Read the template file.
template_file = File.open(template_file_name)
template_string = template_file.read
template_file.close

# Repeat the next steps for all replicates.
number_of_replicates.times do |r|
	rep_id = "r#{(r+1).to_s.rjust(3).gsub(" ","0")}"

	# Add the species richness of Euclichthyidae and Macrouroidinae randomly to any of the other Gadiform clades.
	clades2.size.times do |x|
		if clades2[x] == "Euclichthyidae" or clades2[x] == "Macrouroidinae"
			clades2[x] = nil
			clade_species_numbers[x] = nil
		end
	end
	clades2.compact!
	clade_species_numbers.compact!
	euclichthyidae_sister = ["Bregmacerotidae","Muraenolepididae","Bathygadinae","Trachyrincinae","Macrourinae","Moridae","Melanonidae","Gadinae","Lotinae","Phycinae","Merlucciidae"].sample
	macrouroidinae_sister = ["Bregmacerotidae","Muraenolepididae","Bathygadinae","Trachyrincinae","Macrourinae","Moridae","Melanonidae","Gadinae","Lotinae","Phycinae","Merlucciidae"].sample
	clades2.size.times do |x|
		if clades2[x] == euclichthyidae_sister
			clade_species_numbers[x] += 1
		end
		if clades2[x] == macrouroidinae_sister
			clade_species_numbers[x] += 2
		end
	end

	# Prepare the BAMM sampling probability file. 
	bamm_sampling_probabilities_string = "1.0\n"
	species.size.times do |x|
		index2 = clades2.index(clades1[x])
		count1 = clades1.count(clades1[x])
		if index2 == nil
			raise "ERROR: #{clades1[x]} could not be found!"
		end
		proportion = count1/clade_species_numbers[index2].to_f
		bamm_sampling_probabilities_string << "#{species[x]}\t#{clades1[x]}\t#{proportion}\n"
	end

	# Write BAMM sampling probability files.
	bamm_sampling_probabilities_file_name = "sample_probs.txt"
	bamm_sampling_probabilities_file = File.open("#{bamm_dir}/#{rep_id}/#{bamm_sampling_probabilities_file_name}","w")
	bamm_sampling_probabilities_file.write(bamm_sampling_probabilities_string)
	bamm_sampling_probabilities_file.close

	# Write BAMM control files.
	bamm_control_file_name = "control.txt"
	bamm_control_file = File.open("#{bamm_dir}/#{rep_id}/#{bamm_control_file_name}","w")
	tree_file_name = ""
	Dir.entries("#{bamm_dir}/#{rep_id}").each {|e| tree_file_name = e if e.match(/\.tre/)}
	new_template_string = template_string.sub("XXX",tree_file_name)
	new_template_string.sub!("YYY",(rand(10000)).to_s)
	bamm_control_file.write(new_template_string)
	bamm_control_file.close

	# Copy the BAMM executable to the BAMM directory.
	FileUtils.cp(bamm_executable_name,"#{bamm_dir}/#{rep_id}/")

end
