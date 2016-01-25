## Michael Matschiner 2015-03-29.

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Collect names of nucleotide fasta files in the input directory.
mitochondrial_marker_ids = ["ATP6","ATP8","COX1","COX2","COX3","CYTB","ND1","ND2","ND3","ND4","ND4L","ND5"]
dir_entries_in = Dir.entries(alignment_directory_in)
filenames_in = []
dir_entries_in.each do |e|
	if e.match(/.*_nucl.fasta/) or mitochondrial_marker_ids.include?(e.chomp(".fasta"))
		filenames_in << e
	end
end

# Do for each fasta file in the input directory.
filenames_in.each do |f|

	# Read the fasta file.
	fasta_file = File.open("#{alignment_directory_in}/#{f}")
	fasta_lines = fasta_file.readlines
	fasta_ids = []
	fasta_seqs = []
	fasta_lines.each do |l|
		if l[0] == ">"
			fasta_ids << l[1..-1].strip
			fasta_seqs << ""
		else
			fasta_seqs.last << l.strip
		end
	end

	# Prepare the string for a new fasta file.
	new_fasta_string = ""
	fasta_ids.size.times do |x|
		new_fasta_string << ">#{fasta_ids[x]}\n"
		fasta_seqs[x].size.times do |pos|
			unless ((pos-2)/3)*3 == (pos-2)
				new_fasta_string << fasta_seqs[x][pos]
			end
		end
		new_fasta_string << "\n"
	end

	# Write the filtered fasta file.
	new_fasta_file = File.open("#{alignment_directory_out}/#{f}","w")
	new_fasta_file.write(new_fasta_string)

end
