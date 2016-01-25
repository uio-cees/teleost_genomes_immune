# Michael Matschiner, 2015-06-23.

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
bin_directory_in = ARGV[1].chomp("/")
tree_directory_out = ARGV[2].chomp("/")

# Create the output directories if they don't exist yet.
FileUtils.mkdir_p(tree_directory_out) unless Dir.exists?(tree_directory_out)

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
bins.size.times do |z|

	# Create the output directories if they don't exist yet.
	Dir.mkdir("#{tree_directory_out}/#{bin_names[z]}") unless Dir.exists?("#{tree_directory_out}/#{bin_names[z]}")

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

	# Prepare the string for the concatenated phylip file.
	concatenated_string = "#{nexus_ids_per_alignment[0].size} #{concatenated_seqs[0].size}\n"
	nexus_ids_per_alignment[0].size.times do |x|
	    concatenated_string << "#{nexus_ids_per_alignment[0][x].ljust(12)}#{concatenated_seqs[x]}\n"
	end

	# Write the phylip file.
	concatenated_file_name = "#{tree_directory_out}/#{bin_names[z]}/concatenated.phy"
	concatenated_file = File.open(concatenated_file_name,"w")
	concatenated_file.write(concatenated_string)
	puts "Wrote #{tree_directory_out}/#{bin_names[z]}/concatenated.phy"

	# Prepare partitions strings.
	cp1cp2_partitions_string = ""
	cp1cp2_partitions_string << "DNA, cp1 = 1-#{concatenated_seqs[0].size}\\2\n"
	cp1cp2_partitions_string << "DNA, cp2 = 2-#{concatenated_seqs[0].size}\\2\n"

	# Write the partitions files.
	cp1cp2_partitions_file_name = "#{tree_directory_out}/#{bin_names[z]}/concatenated.parts"
	cp1cp2_partitions_file = File.open(cp1cp2_partitions_file_name,"w")
	cp1cp2_partitions_file.write(cp1cp2_partitions_string)
	cp1cp2_partitions_file.close

end
