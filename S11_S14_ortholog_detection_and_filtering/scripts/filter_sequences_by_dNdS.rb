## Michael Matschiner 2015-02-26.

require "./resources/array_stats.rb"

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
omega_threshold = ARGV[2].to_f

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Prepare a codeml.ctl file.
control_string = ""
control_string << "      seqfile = tmp.phy   * sequence data file name\n"
control_string << "      outfile = results.txt   * main result file name\n"
control_string << "\n"
control_string << "        noisy = 0      * 0,1,2,3,9: how much rubbish on the screen\n"
control_string << "      verbose = 0      * 1:detailed output\n"
control_string << "      runmode = -2     * -2:pairwise\n"
control_string << "\n"
control_string << "      seqtype = 1      * 1:codons\n"
control_string << "    CodonFreq = 1      * 0:equal, 1:F1X4, 2:F3X4, 3:F61\n"
control_string << "        model = 0      *\n"
control_string << "      NSsites = 0      *\n"
control_string << "        icode = 0      * 0:universal code\n"
control_string << "\n"
control_string << "    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated\n"
control_string << "        kappa = 1      * initial or fixed kappa\n"
control_string << "\n"
control_string << "    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated\n"
control_string << "        omega = 0.5    * initial omega value\n"

# Write the control file to the current directory.
control_file = File.open("codeml.ctl","w")
control_file.write(control_string)
control_file.close

# Collect names of nucleotide fasta files in the input directory.
dir_entries_in = Dir.entries("#{alignment_directory_in}")
filenames_in = []
dir_entries_in.each {|e| filenames_in << e if e.match(/.*_nucl.fasta/)}

# Open a file to which all omega estimates will be written.
estimates_file = File.open("omega_estimates.txt","a")

# Do for each fasta file in the input directory.
filenames_in.each do |f|

	# Feedback.
	print "Analysing file #{f}..."

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

	# For each non-empty sequence except Danrer and Astmex, write a reduced
	# fasta file with only that sequence and Danrer, then run codeml to estimate
	# omega for this 2-sequence comparison.
	danrer_empty = false
	species_ids_compared_with_danrer = []
	omega_estimates = []
	danrer_phylip_string = "2 #{fasta_seqs[0].size}\n"
	fasta_ids.size.times do |x|
		if fasta_ids[x] == "Danrer"
			danrer_empty = true if fasta_seqs[x].match(/^-+$/)
			danrer_phylip_string << "#{fasta_ids[x].ljust(12)}#{fasta_seqs[x]}\n"
		end
	end
	unless danrer_empty
		fasta_ids.size.times do |x|
			if fasta_ids[x] != "Danrer" and fasta_ids[x] != "Astmex"
				unless fasta_seqs[x].match(/^-+$/)

					# Memorize that this species has been used for dNdS estimation.
					species_ids_compared_with_danrer << fasta_ids[x]

					# Add the sequence of this taxon to the fasta string that already contains the Danrer
					# sequence.
					reduced_phylip_string = Marshal.load(Marshal.dump(danrer_phylip_string))
					reduced_phylip_string << "#{fasta_ids[x].ljust(12)}#{fasta_seqs[x]}\n"

					# Write the reduced fasta file with these two sequences.
					reduced_phylip_file = File.open("tmp.phy","w")
					reduced_phylip_file.write(reduced_phylip_string)
					reduced_phylip_file.close

					# Run codeml with codeml.ctl and file tmp.fasta
					system("codeml > /dev/null")
					
					# Read the codeml results file.
					result_file = File.open("results.txt")
					result_lines = result_file.readlines

					# Delete codeml output files.
					File.delete("results.txt")
					File.delete("2ML.dN")
					File.delete("2ML.dS")
					File.delete("2ML.t")
					File.delete("2NG.dN")
					File.delete("2NG.dS")
					File.delete("2NG.t")
					File.delete("rst")
					File.delete("rub")
					File.delete("rst1")

					# Delete the reduced phylip file.
					File.delete("tmp.phy")

					# Make sure codeml has run without errors.
					unless result_lines.size == 89
						raise "Unexpected length of codeml result file!"
					end

					# Parse the result file and store the omega estimate.
					result_lines[-1].match(/dN\/dS=\s+(\d+\.\d+)/)
					if $1 == nil
						puts "The omega estimate could not be parsed!"
						puts result_lines
						exit
					else
						omega_estimates << $1.to_f
					end

				end
			end
		end
	end

	# Prepare the string for a new fasta file.
	new_fasta_string = ""
	fasta_ids.size.times do |x|
		new_fasta_string << ">#{fasta_ids[x]}\n"
		if species_ids_compared_with_danrer.include?(fasta_ids[x]) and omega_estimates[species_ids_compared_with_danrer.index(fasta_ids[x])] > omega_threshold
			fasta_seqs[0].size.times {new_fasta_string << "-"}
			new_fasta_string << "\n"
		else
			new_fasta_string << "#{fasta_seqs[x]}\n"
		end
	end

	# Write the corrected fasta file.
	new_fasta_file = File.open("#{alignment_directory_out}/#{f}","w")
	new_fasta_file.write(new_fasta_string)

	puts " done."

end

# Delete the codeml.ctl file.
File.delete("codeml.ctl")
