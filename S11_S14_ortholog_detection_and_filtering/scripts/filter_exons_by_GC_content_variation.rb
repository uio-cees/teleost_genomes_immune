## Michael Matschiner 2015-02-28.

require "./resources/array_stats.rb"

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
maximum_gc_content_variation = ARGV[2].to_f

# Read the nuclear_queries_exons.txt file from the info directory.
nuclear_queries_exons_info_file_name = "../../info/nuclear_queries_exons.txt"
nuclear_queries_exons_info_file = File.open(nuclear_queries_exons_info_file_name)
nuclear_queries_exons_info_lines = nuclear_queries_exons_info_file.readlines
filtered_nuclear_exons_info_lines = []

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Collect names of nucleotide fasta files in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in)
filenames_in = []
dir_entries_in.each {|e| filenames_in << e if e.match(/.*_nucl.fasta/)}

# Do for each fasta file in the input directory.
filenames_in.each do |f|

	# Read the fasta file.
	fasta_file = File.open("#{alignment_directory_in}/#{f}")
	fasta_lines = fasta_file.readlines
	fasta_ids = []
	fasta_bitscores = []
	fasta_hits = []
	fasta_seqs = []
	fasta_lines.each do |l|
		if l[0] == ">"
			fasta_ids << l[1..-1].strip
			fasta_seqs << ""
		else
			fasta_seqs.last << l.strip
		end
	end

	# Get the variation in GC content in the 1st and 2nd codon position.
	cp1_gc_contents = []
	cp2_gc_contents = []
	fasta_ids.size.times do |x|
		unless fasta_seqs[x].match(/^-+$/)
			cp1_bases = []
			cp2_bases = []
			fasta_seqs[x].size.times do |pos|
				if (pos/2)*2 == pos
					cp1_bases << fasta_seqs[x][pos]
				elsif ((pos-1)/2)*2 == (pos-1)
					cp2_bases << fasta_seqs[x][pos]
				else
					raise "Unexpected codon position!"
				end
			end
			cp1_A_count = cp1_bases.count("A")
			cp1_C_count = cp1_bases.count("C")
			cp1_G_count = cp1_bases.count("G")
			cp1_T_count = cp1_bases.count("T")
			cp1_base_count = cp1_A_count + cp1_C_count + cp1_G_count + cp1_T_count
			cp2_A_count = cp2_bases.count("A")
			cp2_C_count = cp2_bases.count("C")
			cp2_G_count = cp2_bases.count("G")
			cp2_T_count = cp2_bases.count("T")
			cp2_base_count = cp2_A_count + cp2_C_count + cp2_G_count + cp2_T_count
			cp1_gc_contents << (cp1_C_count + cp1_G_count)/cp1_base_count.to_f
			cp2_gc_contents << (cp2_C_count + cp2_G_count)/cp2_base_count.to_f
		end
	end

	if cp1_gc_contents.standard_deviation <= maximum_gc_content_variation and cp2_gc_contents.standard_deviation <= maximum_gc_content_variation

		# Prepare the string for a new fasta file.
		new_fasta_string = ""
		fasta_ids.size.times do |x|
			new_fasta_string << ">#{fasta_ids[x]}\n"
			new_fasta_string << "#{fasta_seqs[x]}\n"
		end

		# Write the new fasta file.
		new_fasta_file = File.open("#{alignment_directory_out}/#{f}","w")
		new_fasta_file.write(new_fasta_string)

		# Add to the filtered nuclear exons info lines.
		exon_id = f.split("_")[0]
		nuclear_queries_exons_info_lines.each do |l|
			filtered_nuclear_exons_info_lines << l if l.include?(exon_id)
		end

	end

end

# Prepare the string for the filtered nuclear exons file.
filtered_nuclear_exons_info_string = "exon_id\tgene_id\ttranscript_id\ttranslation\tlength\tthreshold_bitscore\tcorrect_bitscores\tincorrect_bitscores\tmin(correct_bitscores)-max(incorrect_bitscores)\n"
filtered_nuclear_exons_info_lines.each {|l| filtered_nuclear_exons_info_string << l}

# Write the filtered nuclear exons file to the info directory.
filtered_nuclear_exons_info_file_name = "../../info/filtered_nuclear_exons.txt"
filtered_nuclear_exons_info_file = File.open(filtered_nuclear_exons_info_file_name,"w")
filtered_nuclear_exons_info_file.write(filtered_nuclear_exons_info_string)
filtered_nuclear_exons_info_file.close
