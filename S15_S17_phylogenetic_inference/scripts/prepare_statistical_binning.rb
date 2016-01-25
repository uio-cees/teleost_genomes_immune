# Michael Matschiner, 2015-04-13

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Get the input and output directories.
alignment_input_directory = ARGV[0]
statistical_binning_directory = ARGV[1]

# Get a list of all nexus file in the input directory.
dir_entries = Dir.entries(alignment_input_directory)
nexus_file_names = []
dir_entries.each {|e| nexus_file_names << e if e.match(/.nex$/)}

# Read each nexus file name.
nexus_file_names.each do |f|

	# Determine the marker name.
	marker_name = f.chomp(".nex")

	# Make the marker directory if it does not exist yet.
	unless Dir.exists?("#{statistical_binning_directory}/#{marker_name}")
		FileUtils.mkdir_p("#{statistical_binning_directory}/#{marker_name}")
	end

	# Read nexus file.
	nexus_file = File.open("#{alignment_input_directory}/#{f}")
	nexus_lines = nexus_file.readlines
	nexus_ids = []
	nexus_seqs = []
	in_matrix = false
	nexus_lines.each do |l|
		if l.downcase.strip == "matrix"
			in_matrix = true
		elsif l.strip == ";"
			in_matrix = false
		elsif l.strip != "" and in_matrix == true
			line_ary = l.split
			nexus_ids << line_ary[0]
			nexus_seqs << line_ary[1]
		end
	end

	# Prepare the fasta string.
	fasta_string = ""
	nexus_ids.size.times do |x|
		fasta_string << ">#{nexus_ids[x]}\n"
		fasta_string << "#{nexus_seqs[x]}\n"
	end

	# Write the fasta file.
	fasta_file_name = "#{marker_name}.fasta"
	fasta_file = File.open("#{statistical_binning_directory}/#{marker_name}/#{fasta_file_name}","w")
	fasta_file.write(fasta_string)
	fasta_file.close

	# Prepare a partitions string.
	partitions_string = "DNA, cp1 = 1-#{nexus_seqs[0].size}\\2\n"
	partitions_string << "DNA, cp1 = 2-#{nexus_seqs[0].size}\\2\n"

	# Write the partitions file.
	partitions_file_name = "partitions.txt"
	partitions_file = File.open("#{statistical_binning_directory}/#{marker_name}/#{partitions_file_name}","w")
	partitions_file.write(partitions_string)
	partitions_file.close

	# Prepare a constraint string.
	constraint_string = "(("
	nexus_ids.each do |i|
		if i == "Astmex" or i == "Danrer"
			constraint_string << "#{i},"
		end
	end
	constraint_string.chomp!(",")
	constraint_string << "),("
	nexus_ids.each do |i|
		if i != "Astmex" and i != "Danrer"
			constraint_string << "#{i},"
		end
	end
	constraint_string.chomp!(",")
	constraint_string << "));\n"

	# Write the constraint file.
	constraint_file_name = "constraint.txt"
	constraint_file = File.open("#{statistical_binning_directory}/#{marker_name}/#{constraint_file_name}","w")
	constraint_file.write(constraint_string)
	constraint_file.close

	# Run RAxML via run_raxml.rb with the fasta file.
	constraint_file_name_w_path = "#{statistical_binning_directory}/#{marker_name}/#{constraint_file_name}"
	partitions_file_name_w_path = "#{statistical_binning_directory}/#{marker_name}/#{partitions_file_name}"
	fasta_file_name_w_path = "#{statistical_binning_directory}/#{marker_name}/#{fasta_file_name}"
	tree_file_name_w_path = "#{statistical_binning_directory}/#{marker_name}/ML.tre"
	system("python3 resources/run_raxml.py -b auto -q #{partitions_file_name_w_path} -c #{constraint_file_name_w_path} --one-category #{fasta_file_name_w_path} #{tree_file_name_w_path}")
	system("rm #{partitions_file_name_w_path}")
	system("rm #{constraint_file_name_w_path}")

end