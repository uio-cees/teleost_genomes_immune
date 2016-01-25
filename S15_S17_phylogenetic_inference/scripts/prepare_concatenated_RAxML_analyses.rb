# Michael Matschiner, 2015-03-23.

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
tree_directory_out = ARGV[1].chomp("/")

# Create the output directory if it does not exist yet.
unless Dir.exists?(tree_directory_out)
	Dir.mkdir(tree_directory_out)
end

# Collect names of nexus files in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in)
filenames_in = []
dir_entries_in.each {|e| filenames_in << e if e.match(/.*.nex/)}

# Initiate arrays for the ids and sequences of each alignment.
nexus_ids_per_alignment = []
nexus_seqs_per_alignment = []

# Do for each fasta file in the input directory.
filenames_in.each do |f|

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
			nexus_ids << l.strip.split[0]
			nexus_seqs << l.strip.split[1]
		end
	end
	nexus_ids_per_alignment << nexus_ids
	nexus_seqs_per_alignment << nexus_seqs
end

# Make sure the ids are identical in all alignments.
nexus_ids_per_alignment[1..-1].each do |ids|
	raise "Ids differ between alignments!" if ids != nexus_ids_per_alignment[0]
end

# Prepare concatenated alignments.
cp1cp2_concatenated_seqs = []
cp1cp2cp3_concatenated_seqs = []
cp1cp2cp3ry_concatenated_seqs = []
nexus_ids_per_alignment[0].size.times do |x|
	cp1cp2_concatenated_seqs << ""
	cp1cp2cp3_concatenated_seqs << ""
	cp1cp2cp3ry_concatenated_seqs << ""
	marker_count = 0
	nexus_seqs_per_alignment.each do |s|
		codon = ""
		s[x].size.times do |pos|
			codon << s[x][pos].downcase
			if codon.size == 3
				if codon == "---"
					cp1cp2_concatenated_seqs.last << "--"
					cp1cp2cp3_concatenated_seqs.last << "---"
					cp1cp2cp3ry_concatenated_seqs.last << "---"
				else
					cp1cp2_concatenated_seqs.last << codon[0..1]
					cp1cp2cp3_concatenated_seqs.last << codon
					cp1cp2cp3ry_concatenated_seqs.last << codon[0..1]
					if ["a","g","r"].include?(codon[2])
						cp1cp2cp3ry_concatenated_seqs.last << "r"
					elsif ["c","t","y"].include?(codon[2])
						cp1cp2cp3ry_concatenated_seqs.last << "y"
					elsif ["n"].include?(codon[2])
						cp1cp2cp3ry_concatenated_seqs.last << "n"
					else
						raise "ERROR: Unexpected third codon position in codon \"#{codon}\" (pos = #{pos}) of #{nexus_ids_per_alignment[0][x]} in #{filenames_in[marker_count]}!"
					end
				end
				codon = ""
			end
		end
		marker_count += 1
	end
end

# Get the maximum name length.
max_name_length = 0
nexus_ids_per_alignment[0].each do |i|
	max_name_length = i.size if i.size > max_name_length
end

# Prepare the string for a concatenated phylip file.
cp1cp2_phylip_string = "#{nexus_ids_per_alignment[0].size} #{cp1cp2_concatenated_seqs[0].size}\n"
cp1cp2cp3_phylip_string = "#{nexus_ids_per_alignment[0].size} #{cp1cp2cp3_concatenated_seqs[0].size}\n"
cp1cp2cp3ry_phylip_string = "#{nexus_ids_per_alignment[0].size} #{cp1cp2cp3ry_concatenated_seqs[0].size}\n"
nexus_ids_per_alignment[0].size.times do |x|
	cp1cp2_phylip_string << "#{nexus_ids_per_alignment[0][x].ljust(max_name_length+2)} #{cp1cp2_concatenated_seqs[x]}\n"
	cp1cp2cp3_phylip_string << "#{nexus_ids_per_alignment[0][x].ljust(max_name_length+2)} #{cp1cp2cp3_concatenated_seqs[x]}\n"
	cp1cp2cp3ry_phylip_string << "#{nexus_ids_per_alignment[0][x].ljust(max_name_length+2)} #{cp1cp2cp3ry_concatenated_seqs[x]}\n"
end

# Write the concatenated phylip files.
cp1cp2_phylip_file_name = "#{tree_directory_out}/cp1cp2_concatenated.phy"
cp1cp2_phylip_file = File.open(cp1cp2_phylip_file_name,"w")
cp1cp2_phylip_file.write(cp1cp2_phylip_string)
cp1cp2_phylip_file.close
cp1cp2cp3_phylip_file_name = "#{tree_directory_out}/cp1cp2cp3_concatenated.phy"
cp1cp2cp3_phylip_file = File.open(cp1cp2cp3_phylip_file_name,"w")
cp1cp2cp3_phylip_file.write(cp1cp2cp3_phylip_string)
cp1cp2cp3_phylip_file.close
cp1cp2cp3ry_phylip_file_name = "#{tree_directory_out}/cp1cp2cp3ry_concatenated.phy"
cp1cp2cp3ry_phylip_file = File.open(cp1cp2cp3ry_phylip_file_name,"w")
cp1cp2cp3ry_phylip_file.write(cp1cp2cp3ry_phylip_string)
cp1cp2cp3ry_phylip_file.close

# Prepare partitions strings.
cp1cp2_partitions_string = ""
cp1cp2_partitions_string << "DNA, cp1 = 1-#{cp1cp2_concatenated_seqs[0].size}\\2\n"
cp1cp2_partitions_string << "DNA, cp2 = 2-#{cp1cp2_concatenated_seqs[0].size}\\2\n"
cp1cp2cp3_partitions_string = ""
cp1cp2cp3_partitions_string << "DNA, cp1 = 1-#{cp1cp2cp3_concatenated_seqs[0].size}\\3\n"
cp1cp2cp3_partitions_string << "DNA, cp2 = 2-#{cp1cp2cp3_concatenated_seqs[0].size}\\3\n"
cp1cp2cp3_partitions_string << "DNA, cp3 = 3-#{cp1cp2cp3_concatenated_seqs[0].size}\\3\n"

# Write the partitions files.
cp1cp2_partitions_file_name = "#{tree_directory_out}/cp1cp2_concatenated.parts"
cp1cp2_partitions_file = File.open(cp1cp2_partitions_file_name,"w")
cp1cp2_partitions_file.write(cp1cp2_partitions_string)
cp1cp2_partitions_file.close
cp1cp2cp3_partitions_file_name = "#{tree_directory_out}/cp1cp2cp3_concatenated.parts"
cp1cp2cp3_partitions_file = File.open(cp1cp2cp3_partitions_file_name,"w")
cp1cp2cp3_partitions_file.write(cp1cp2cp3_partitions_string)
cp1cp2cp3_partitions_file.close
cp1cp2cp3ry_partitions_file_name = "#{tree_directory_out}/cp1cp2cp3ry_concatenated.parts"
cp1cp2cp3ry_partitions_file = File.open(cp1cp2cp3ry_partitions_file_name,"w")
cp1cp2cp3ry_partitions_file.write(cp1cp2cp3_partitions_string)
cp1cp2cp3ry_partitions_file.close
