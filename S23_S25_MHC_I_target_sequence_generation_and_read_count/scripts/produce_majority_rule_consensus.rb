# This script produces majority rule consensus sequences from aligned fasta files.
# Run this with like
# ruby produce_majority_rule_consensus.rb -f inputfilename -o outputfilename

class String
	def identical_with(str)
		return false unless str.class == String
		str1 = self.downcase
		str2 = str.downcase
		identical = true
		if str1.size != str2.size
			identical = false
		else
			str1.size.times do |x|
				if str1[x] != "n" and str1[x] != "-" and str1[x] != "?" and str2[x] != "n" and str2[x] != "-" and str2[x] != "?"
					identical = false unless str1[x] == str2[x]
				end
			end
		end
		identical
	end

	def informative_sites
		count = 0
		str = self.downcase
		str.size.times {|x| count += 1 unless str[x] == "n" or str[x] == "-" or str[x] == "?"}
		count
	end

	def more_informative_than(str)
		return true unless str.class == String
		str1 = self.downcase
		str2 = str.downcase
		more_informative = false
		more_informative = true if str1.informative_sites > str.informative_sites
		more_informative
	end
end

# Read the command line arguments.
if ARGV.include?("-f")
	fasta_file_name = ARGV[ARGV.index("-f")+1]
else
	raise "The fasta file name should be specified with the '-f' option!"
end
if ARGV.include?("-o")
	output_file_name = ARGV[ARGV.index("-o")+1]
else
	raise "The output file name should be specified with the '-o' option!"
end

# Open the fasta file.
fasta_file = File.open(fasta_file_name)
fasta_lines = fasta_file.readlines
fasta_file.close

# Read the fasta file.
seqs = []
fasta_lines.each do |l|
	line = l.strip
	if line[0..0] == ">"
		seqs << ""
	elsif line != ""
		seqs.last << line
	end
end

# Make sure all sequences have the same length.
1.upto(seqs.size-1) {|x| raise "Sequences in file #{fasta_file_name} do not have the same length!" unless seqs[x].size == seqs[0].size}

# Remove identical sequences.
seqs.uniq!

# Remove sequences that are identical with others, taking into account gaps and ns.
0.upto(seqs.size-2) do |x|
	(x+1).upto(seqs.size-1) do |y|
		unless seqs[x] == nil
			if seqs[x].identical_with(seqs[y])
				if seqs[x].more_informative_than(seqs[y])
					seqs[y] = nil
				else
					seqs[x] = nil
				end
			end
		end
	end
end
seqs.compact!

consensus = ""
seqs[0].size.times do |pos|
	characters_at_this_site = []
	frequencies_of_characters_at_this_site = []
	only_gaps_at_this_pos = true
	seqs.size.times {|x| only_gaps_at_this_pos = false if seqs[x][pos..pos] != "-"}
	if only_gaps_at_this_pos
		consensus << "-"
	else
		seqs.size.times do |x|
			character = seqs[x][pos..pos].downcase
			unless ["n","-","?"].include?(character)
				if characters_at_this_site.include?(character)
					frequencies_of_characters_at_this_site[characters_at_this_site.index(character)] += 1
				else
					characters_at_this_site << character
					frequencies_of_characters_at_this_site << 1
				end
			end
		end
		if characters_at_this_site == []
			consensus << "n"
		else
			consensus << characters_at_this_site[frequencies_of_characters_at_this_site.index(frequencies_of_characters_at_this_site.max)]
		end
	end
end

output_string = ">Consensus\n"
while consensus.size > 0 do 
	output_string << "#{consensus.slice!(0..59)}\n"
end
output_string << "\n"


output_file = File.new(output_file_name,"w")
output_file.write(output_string)
output_file.close
