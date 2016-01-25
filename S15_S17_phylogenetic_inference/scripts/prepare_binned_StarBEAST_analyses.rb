# Michael Matschiner, 2015-04-07.

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
bin_directory_in = ARGV[1].chomp("/")
tree_directory_out = ARGV[2].chomp("/")
constraint_file_name = ARGV[3]
starting_tree_directory = ARGV[4]

# If a starting tree directory has been specified, get the name of the starting tree in this directory.
starting_tree_file_name = nil
unless starting_tree_directory == nil
	starting_tree_directory_entries = Dir.entries(starting_tree_directory)
	starting_tree_directory_entries.each do |e|
		if e.match(/.tre$/)
			starting_tree_file_name = "#{starting_tree_directory}/#{e}"
		end
	end
end

# Create the output directories if they don't exist yet.
FileUtils.mkdir_p(tree_directory_out) unless Dir.exists?(tree_directory_out)
Dir.mkdir("#{tree_directory_out}") unless Dir.exists?("#{tree_directory_out}")
Dir.mkdir("#{tree_directory_out}/nexus") unless Dir.exists?("#{tree_directory_out}/nexus")
Dir.mkdir("#{tree_directory_out}/tree_linkage") unless Dir.exists?("#{tree_directory_out}/tree_linkage")
Dir.mkdir("#{tree_directory_out}/replicates") unless Dir.exists?("#{tree_directory_out}/replicates")
Dir.mkdir("#{tree_directory_out}/replicates/r01") unless Dir.exists?("#{tree_directory_out}/replicates/r01")
Dir.mkdir("#{tree_directory_out}/replicates/r02") unless Dir.exists?("#{tree_directory_out}/replicates/r02")
Dir.mkdir("#{tree_directory_out}/replicates/r03") unless Dir.exists?("#{tree_directory_out}/replicates/r03")
Dir.mkdir("#{tree_directory_out}/replicates/r04") unless Dir.exists?("#{tree_directory_out}/replicates/r04")
Dir.mkdir("#{tree_directory_out}/replicates/r05") unless Dir.exists?("#{tree_directory_out}/replicates/r05")

# Collect names of nexus files in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in)
filenames_in = []
dir_entries_in.each {|e| filenames_in << e if e.match(/.*.nex/)}

# Initiate arrays for the ids and sequences of each alignment.
nexus_ids_per_alignment = []
nexus_seqs_per_alignment = []
alignment_gene_names = []

# Do for each fasta file in the input directory.
filenames_in.each do |f|

	# Store the gene name for this alignment.
	alignment_gene_names << f.chomp(".nex")

	# Read the nexus file.
	nexus_file = File.open("#{alignment_directory_in}/#{f}")
	nexus_lines = nexus_file.readlines
	nexus_ids = []
	nexus_seqs = []
	in_matrix = false
	nexus_lines.each do |l|
		if l.strip == "matrix"
			in_matrix = true
		elsif l.strip == ";"
			in_matrix = false
		elsif in_matrix
			nexus_id = l.strip.split[0]
			nexus_seq = l.strip.split[1]
			unless nexus_id == "Gadmor"
				nexus_ids << nexus_id
				nexus_seqs << nexus_seq
			end
		end
	end
	nexus_ids_per_alignment << nexus_ids
	nexus_seqs_per_alignment << nexus_seqs
end

# Make sure the ids are identical in all alignments.
nexus_ids_per_alignment[1..-1].each do |ids|
	raise "Ids differ between alignments!" if ids != nexus_ids_per_alignment[0]
end

# Collect names of bin files from the bin directory.
dir_entries_in = Dir.entries(bin_directory_in)
file_names_in = []
dir_entries_in.each {|e| file_names_in << e if e.match(/bin\.\d+\.txt/)}

# Read and store bins and the genes they contain.
bins = []
bin_names = []
file_names_in.each do |f|
	bin_name = f.chomp(".txt")
	bin_file = File.open("#{bin_directory_in}/#{f}")
	bin_lines = bin_file.readlines
	bin = []
	bin_lines.each do |l|
		bin << l.strip unless l.strip == ""
	end
	bins << bin
	bin_names << bin_name.gsub(".","_")
end

# For each bin, prepare a concatenated alignment of all genes, split this concatenated
# alignment according to codon position, and prepare and write nexus files for this
# split alignments.
tree_linkage_string = ""
bins.size.times do |z|

	# Prepare a concatenated alignment.
	concatenated_seqs = []
	nexus_ids_per_alignment[0].size.times do |x|
		concatenated_seqs << ""
		nexus_seqs_per_alignment.size.times do |y|
			if bins[z].include?(alignment_gene_names[y])
				s = nexus_seqs_per_alignment[y]
				concatenated_seqs.last << s[x]
			end
		end
	end

	# Split the concatenated alignment according to codon position.
	concatenated_cp1_seqs = []
	concatenated_cp2_seqs = []
	concatenated_seqs.each do |s|
		concatenated_cp1_seqs << ""
		concatenated_cp2_seqs << ""
		s.size.times do |x|
			if (x/2)*2 == x
				concatenated_cp1_seqs.last << s[x]
			else
				concatenated_cp2_seqs.last << s[x]
			end
		end
	end

	# Prepare the string for two concatenated nexus files, one for each codon position.
	concatenated_cp1_string = "#nexus\n"
	concatenated_cp1_string << "\n"
	concatenated_cp1_string << "begin data;\n"
	concatenated_cp1_string << "dimensions  ntax=#{nexus_ids_per_alignment[0].size} nchar=#{concatenated_cp1_seqs[0].size};\n"
	concatenated_cp1_string << "format datatype=DNA gap=- missing=?;\n"
	concatenated_cp1_string << "matrix\n"
	nexus_ids_per_alignment[0].size.times do |x|
	    concatenated_cp1_string << "#{nexus_ids_per_alignment[0][x].ljust(12)}#{concatenated_cp1_seqs[x]}\n"
	end
	concatenated_cp1_string << ";\n"
	concatenated_cp1_string << "end;\n"
	concatenated_cp2_string = "#nexus\n"
	concatenated_cp2_string << "\n"
	concatenated_cp2_string << "begin data;\n"
	concatenated_cp2_string << "dimensions  ntax=#{nexus_ids_per_alignment[0].size} nchar=#{concatenated_cp2_seqs[0].size};\n"
	concatenated_cp2_string << "format datatype=DNA gap=- missing=?;\n"
	concatenated_cp2_string << "matrix\n"
	nexus_ids_per_alignment[0].size.times do |x|
	    concatenated_cp2_string << "#{nexus_ids_per_alignment[0][x].ljust(12)}#{concatenated_cp2_seqs[x]}\n"
	end
	concatenated_cp2_string << ";\n"
	concatenated_cp2_string << "end;\n"

	# Write the nexus files.
	concatenated_cp1_file_name = "#{tree_directory_out}/nexus/#{bin_names[z]}_cp1.nex"
	concatenated_cp1_file = File.open(concatenated_cp1_file_name,"w")
	concatenated_cp1_file.write(concatenated_cp1_string)
	puts "Wrote #{tree_directory_out}/nexus/#{bin_names[z]}_cp1.nex"
	concatenated_cp2_file_name = "#{tree_directory_out}/nexus/#{bin_names[z]}_cp2.nex"
	concatenated_cp2_file = File.open(concatenated_cp2_file_name,"w")
	concatenated_cp2_file.write(concatenated_cp2_string)
	puts "Wrote #{tree_directory_out}/nexus/#{bin_names[z]}_cp2.nex"

	# Prepare a tree linkage string.
	tree_linkage_string << "#{bin_names[z]}:#{bin_names[z]}_cp1,#{bin_names[z]}_cp2\n"

end
# Write the tree linkage file.
tree_linkage_file_name = "tree_linkage.txt"
tree_linkage_file = File.open("#{tree_directory_out}/tree_linkage/#{tree_linkage_file_name}","w")
tree_linkage_file.write(tree_linkage_string)
tree_linkage_file.close
puts "Wrote #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name}"

# Produce BEAST XML files with beauti.rb
if constraint_file_name == "../data/constraints/fossil_constraints_starbeast.xml"
	if starting_tree_file_name == nil
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r01 -c #{constraint_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r02 -c #{constraint_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r03 -c #{constraint_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r04 -c #{constraint_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r05 -c #{constraint_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")	
	else
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r01 -c #{constraint_file_name} -t #{starting_tree_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r02 -c #{constraint_file_name} -t #{starting_tree_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r03 -c #{constraint_file_name} -t #{starting_tree_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r04 -c #{constraint_file_name} -t #{starting_tree_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")
		system("ruby resources/beauti.rb -id 76g_nucl_binn_fossils -n #{tree_directory_out.gsub(" ","\\ ")}/nexus -o #{tree_directory_out.gsub(" ","\\ ")}/replicates/r05 -c #{constraint_file_name} -t #{starting_tree_file_name} -*tl #{tree_directory_out}/tree_linkage/#{tree_linkage_file_name} -bd -l 500000000 -m RB -g -e -u -s")	
	end
else
	raise_string = "ERROR: This script is designed for one of the following constraint files: \"../data/constraints/fossil_constraints_starbeast.xml\".\n"
	raise_string << "The specified constraint file \"#{constraint_file_name}\" can not be used."
	raise raise_string
end

# Remove the nexus and tree linkage directories required for beauti.rb.
system("rm -r #{tree_directory_out.gsub(" ","\\ ")}/nexus")
system("rm -r #{tree_directory_out.gsub(" ","\\ ")}/tree_linkage")


# Copy required resources to each directory.
system("cp resources/lib/beast.jar #{tree_directory_out.gsub(" ","\\ ")}/replicates/r01")
system("cp resources/lib/beast.jar #{tree_directory_out.gsub(" ","\\ ")}/replicates/r02")
system("cp resources/lib/beast.jar #{tree_directory_out.gsub(" ","\\ ")}/replicates/r03")
system("cp resources/lib/beast.jar #{tree_directory_out.gsub(" ","\\ ")}/replicates/r04")
system("cp resources/lib/beast.jar #{tree_directory_out.gsub(" ","\\ ")}/replicates/r05")
system("cp -r resources/RBS #{tree_directory_out.gsub(" ","\\ ")}/replicates/r01")
system("cp -r resources/RBS #{tree_directory_out.gsub(" ","\\ ")}/replicates/r02")
system("cp -r resources/RBS #{tree_directory_out.gsub(" ","\\ ")}/replicates/r03")
system("cp -r resources/RBS #{tree_directory_out.gsub(" ","\\ ")}/replicates/r04")
system("cp -r resources/RBS #{tree_directory_out.gsub(" ","\\ ")}/replicates/r05")
if constraint_file_name == "../data/constraints/fossil_constraints_starbeast.xml"
	system("cp -r resources/CA #{tree_directory_out.gsub(" ","\\ ")}/replicates/r01")
	system("cp -r resources/CA #{tree_directory_out.gsub(" ","\\ ")}/replicates/r02")
	system("cp -r resources/CA #{tree_directory_out.gsub(" ","\\ ")}/replicates/r03")
	system("cp -r resources/CA #{tree_directory_out.gsub(" ","\\ ")}/replicates/r04")
	system("cp -r resources/CA #{tree_directory_out.gsub(" ","\\ ")}/replicates/r05")
end

# Change the topology operator weights in each XML file.
["r01","r02","r03","r04","r05"].each do |rep|
	if constraint_file_name == "../data/constraints/fossil_constraints_starbeast.xml"
		xml_file = File.open("#{tree_directory_out}/replicates/#{rep}/76g_nucl_binn_fossils.xml")
	else
		raise_string = "ERROR: This script is designed for one of the following constraint files: \"../data/constraints/fossil_constraints_starbeast.xml\".\n"
		raise_string << "The specified constraint file \"#{constraint_file_name}\" can not be used."
		raise raise_string
	end
	xml_string = xml_file.read
	xml_file.close
	xml_string.sub!("<operator id=\"treeSubtreeSlide:Species\" spec=\"SubtreeSlide\" tree=\"@tree.t:Species\" weight=\"5.0\"/>","<operator id=\"treeSubtreeSlide:Species\" spec=\"SubtreeSlide\" tree=\"@tree.t:Species\" weight=\"50.0\"/>")
	xml_string.sub!("<operator id=\"treeExchange:Species\" spec=\"Exchange\" tree=\"@tree.t:Species\" weight=\"5\.0\"\/>","<operator id=\"treeExchange:Species\" spec=\"Exchange\" tree=\"@tree.t:Species\" weight=\"50\.0\"\/>")
	xml_string.sub!("<operator id=\"treeNarrowExchange:Species\" isNarrow=\"false\" spec=\"Exchange\" tree=\"@tree.t:Species\" weight=\"5\.0\"\/>","<operator id=\"treeNarrowExchange:Species\" isNarrow=\"false\" spec=\"Exchange\" tree=\"@tree.t:Species\" weight=\"50\.0\"\/>")
	xml_string.sub!("<operator id=\"treeWilsonBalding:Species\" spec=\"WilsonBalding\" tree=\"@tree.t:Species\" weight=\"3\.0\"\/>","<operator id=\"treeWilsonBalding:Species\" spec=\"WilsonBalding\" tree=\"@tree.t:Species\" weight=\"30\.0\"\/>")
	xml_string.gsub!("spec=\"ScaleOperator\" weight=\"0.1\"/>","spec=\"ScaleOperator\" weight=\"0.02\"/>")
	xml_string.gsub!("spec=\"DeltaExchangeOperator\" weight=\"0.1\"/>","spec=\"DeltaExchangeOperator\" weight=\"0.02\"/>")
	xml_string.gsub!("spec=\"RBScaleOperator\" weight=\"1.0\"/>","spec=\"RBScaleOperator\" weight=\"0.2\"/>")
	xml_string.gsub!("spec=\"RBOperator\" weight=\"1.0\"/>","spec=\"RBOperator\" weight=\"0.2\"/>")
	if constraint_file_name == "../data/constraints/fossil_constraints_starbeast.xml"
		xml_file = File.open("#{tree_directory_out}/replicates/#{rep}/76g_nucl_binn_fossils.xml","w")
	else
		raise_string = "ERROR: This script is designed for one of the following constraint files: \"../data/constraints/fossil_constraints_starbeast.xml\".\n"
		raise_string << "The specified constraint file \"#{constraint_file_name}\" can not be used."
		raise raise_string
	end
	xml_file.write(xml_string)
	xml_file.close
end

# Edit the start files.
["r01","r02","r03","r04","r05"].each do |rep|
	start_file = File.open("#{tree_directory_out}/replicates/#{rep}/start.sh")
	start_string = start_file.read
	start_file.close
	start_string.sub!("cp beast.jar $SCRATCH","cp beast.jar $SCRATCH\ncp -r RBS $SCRATCH\ncp -r CA $SCRATCH")
	start_file = File.open("#{tree_directory_out}/replicates/#{rep}/start.sh","w")
	start_file.write(start_string)
	start_file.close
end

# Produce start files for the accelerated abel partition.
["r01","r02","r03","r04","r05"].each do |rep|
	start_file = File.open("#{tree_directory_out}/replicates/#{rep}/start.sh")
	start_string = start_file.read
	start_file.close
	start_string.sub!("ntasks-per-node=4","ntasks-per-node=8")
	start_string.sub!("nodes=1","nodes=1\n#SBATCH --partition=accel\n#SBATCH --gres=gpu:2")
	start_string.sub!("module load beagle","module load beagle\nmodule load beagle_gpu")
	start_string.sub!("-threads 3","-threads 5")
	start_string.sub!("-beagle","-beagle -beagle_GPU -beagle_order 1,2")
	start_file = File.open("#{tree_directory_out}/replicates/#{rep}/start_accel.sh","w")
	start_file.write(start_string)
	start_file.close
end
