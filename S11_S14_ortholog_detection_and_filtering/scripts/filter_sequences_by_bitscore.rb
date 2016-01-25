## Michael Matschiner 2015-02-25.

require "./resources/array_stats.rb"

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
relative_bitscore_threshold = ARGV[2].to_f

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
			fasta_header = l[1..-1].strip
			fasta_ids << fasta_header.sub(/\[.*\]/,"")
			fasta_header.match(/bitscore=(.+?)[,\]]/)
			fasta_bitscores << $1
			fasta_header.match(/nhits=(.+?)[,\]]/)
			fasta_hits << $1.to_i
			fasta_seqs << ""
		else
			fasta_seqs.last << l.strip
		end
	end

	# Determine the minimum bitscore needed to qualify as ortholog.
	# It is calculates as relative_bitscore_threshold * maximum bitscore (excluding Danrer and Astmex).
	maximum_ingroup_bitscore = 0
	fasta_bitscores[2..-1].each do |bs|
		unless bs == "None"
			if maximum_ingroup_bitscore < bs.to_f
				maximum_ingroup_bitscore = bs.to_f
			end
		end
	end
	minimum_bitscore = relative_bitscore_threshold * maximum_ingroup_bitscore

	# Prepare the string for a new fasta file.
	new_fasta_string = ""
	fasta_ids.size.times do |x|
		new_fasta_string << ">#{fasta_ids[x]}\n"
		if fasta_bitscores[x] == "None" or fasta_bitscores[x].to_f < minimum_bitscore
			fasta_seqs[0].size.times {new_fasta_string << "-"}
			new_fasta_string << "\n"
		else
			new_fasta_string << "#{fasta_seqs[x]}\n"
		end
	end

	# Write the corrected fasta file.
	new_fasta_file = File.open("#{alignment_directory_out}/#{f}","w")
	new_fasta_file.write(new_fasta_string)

end
