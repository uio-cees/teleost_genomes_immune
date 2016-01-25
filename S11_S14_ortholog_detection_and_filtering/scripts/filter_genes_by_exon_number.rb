## Michael Matschiner 2015-02-27.

require "./resources/array_stats.rb"

class String
	def reverse_complement
		reverse_complement = ""
		(self.size-1).downto(0) do |x|
			if self[x].upcase == "A"
				reverse_complement << "T"
			elsif self[x].upcase == "C"
				reverse_complement << "G"
			elsif self[x].upcase == "G"
				reverse_complement << "C"
			elsif self[x].upcase == "T"
				reverse_complement << "A"
			else
				reverse_complement << "N"
			end
		end
		reverse_complement
	end
	def translate
		translation = ""
		(self.size/3).times do |x|
			codon = self[3*x..3*(x+1)-1].upcase
			if codon == "AAA"
				translation << "K"
			elsif codon == "AAC"
				translation << "N"
			elsif codon == "AAG"
				translation << "K"
			elsif codon == "AAT"
				translation << "N"
			elsif codon == "ACA"
				translation << "T"
			elsif codon == "ACC"
				translation << "T"
			elsif codon == "ACG"
				translation << "T"
			elsif codon == "ACT"
				translation << "T"
			elsif codon == "AGA"
				translation << "R"
			elsif codon == "AGC"
				translation << "S"
			elsif codon == "AGG"
				translation << "R"
			elsif codon == "AGT"
				translation << "S"
			elsif codon == "ATA"
				translation << "I"
			elsif codon == "ATC"
				translation << "I"
			elsif codon == "ATG"
				translation << "M"
			elsif codon == "ATT"
				translation << "I"
			elsif codon == "CAA"
				translation << "Q"
			elsif codon == "CAC"
				translation << "H"
			elsif codon == "CAG"
				translation << "Q"
			elsif codon == "CAT"
				translation << "H"
			elsif codon == "CCA"
				translation << "P"
			elsif codon == "CCC"
				translation << "P"
			elsif codon == "CCG"
				translation << "P"
			elsif codon == "CCT"
				translation << "P"
			elsif codon == "CGA"
				translation << "R"
			elsif codon == "CGC"
				translation << "R"
			elsif codon == "CGG"
				translation << "R"
			elsif codon == "CGT"
				translation << "R"
			elsif codon == "CTA"
				translation << "L"
			elsif codon == "CTC"
				translation << "L"
			elsif codon == "CTG"
				translation << "L"
			elsif codon == "CTT"
				translation << "L"
			elsif codon == "GAA"
				translation << "E"
			elsif codon == "GAC"
				translation << "D"
			elsif codon == "GAG"
				translation << "E"
			elsif codon == "GAT"
				translation << "D"
			elsif codon == "GCA"
				translation << "A"
			elsif codon == "GCC"
				translation << "A"
			elsif codon == "GCG"
				translation << "A"
			elsif codon == "GCT"
				translation << "A"
			elsif codon == "GGA"
				translation << "G"
			elsif codon == "GGC"
				translation << "G"
			elsif codon == "GGG"
				translation << "G"
			elsif codon == "GGT"
				translation << "G"
			elsif codon == "GTA"
				translation << "V"
			elsif codon == "GTC"
				translation << "V"
			elsif codon == "GTG"
				translation << "V"
			elsif codon == "GTT"
				translation << "V"
			elsif codon == "TAA"
				translation << "*"
			elsif codon == "TAC"
				translation << "Y"
			elsif codon == "TAG"
				translation << "*"
			elsif codon == "TAT"
				translation << "Y"
			elsif codon == "TCA"
				translation << "S"
			elsif codon == "TCC"
				translation << "S"
			elsif codon == "TCG"
				translation << "S"
			elsif codon == "TCT"
				translation << "S"
			elsif codon == "TGA"
				translation << "*"
			elsif codon == "TGC"
				translation << "C"
			elsif codon == "TGG"
				translation << "W"
			elsif codon == "TGT"
				translation << "C"
			elsif codon == "TTA"
				translation << "L"
			elsif codon == "TTC"
				translation << "F"
			elsif codon == "TTG"
				translation << "L"
			elsif codon == "TTT"
				translation << "F"
			else
				translation << "X"
			end
		end
		translation
	end
end

class Gene
	attr_reader :id, :transcripts, :chromosome, :gene_tree, :ortholog_ids
	def initialize(id)
		@id = id
		@transcripts = []
		@ortholog_ids = []
	end
	def add_chromosome(chromosome)
		@chromosome = chromosome
	end
	def add_transcript(transcript)
		@transcripts << transcript
	end
	def longest_transcript
		longest_transcript = @transcripts[0]
		@transcripts[1..-1].each do |t|
			if t.length > longest_transcript.length
				longest_transcript = t
			end
		end
		longest_transcript
	end
	def selected
		self.longest_transcript.selected
	end
	def add_gene_tree(gene_tree)
		@gene_tree = gene_tree
	end
	def add_ortholog_id(ortholog_id)
		@ortholog_ids << ortholog_id
	end
end

class Transcript
	attr_reader :id, :exons, :selected, :strand
	def initialize(id)
		@id = id
		@exons = []
	end
	def add_start(transcript_start)
		@start = transcript_start
	end
	def add_end(transcript_end)
		@end = transcript_end
	end
	def add_strand(strand)
		@strand = strand
	end
	def add_exon(exon)
		@exons << exon
	end
	def length
		length = 0
		@exons.each do |e|
			length += e.length
		end
		length
	end
	def long_exons
		long_exons = []
		@exons.each do |e|
			long_exons << e if e.length >= 150
		end
		long_exons
	end
	def selected
		number_of_selected_exons = 0
		@exons.each do |e|
			number_of_selected_exons += 1 if e.selected
		end
		if number_of_selected_exons >= 5
			true
		else
			false
		end
	end
	def sort_exons
		exon_starts = []
		@exons.each do |e|
			exon_starts << e.start
		end
		exon_starts.sort!
		if strand == -1
			exon_starts.reverse!
		end
		tmp_exons = []
		if exon_starts.size != exon_starts.uniq.size
			raise "Some exons seem to have the same start!"
		end
		exon_starts.each do |s|
			@exons.each do |e|
				if e.start == s
					tmp_exons << e
				end
			end
		end
		@exons = tmp_exons
	end
end

class Exon
	attr_reader :id, :selected, :seq, :start, :end, :phase, :translation, :correct_bitscores, :incorrect_bitscores
	def initialize(id)
		@id = id
		@selected = false
	end
	def add_start(exon_start)
		@start = exon_start
	end
	def add_end(exon_end)
		@end = exon_end
	end
	def add_phase(phase)
		@phase = phase
	end
	def length
		(@end - @start).abs + 1
	end
	def select
		@selected = true
	end
	def deselect
		@selected = false
	end
	def add_seq(exon_seq)
		@seq = exon_seq
	end
	def translate
		@translation = @seq.translate
	end
	def add_correct_bitscores(correct_bitscores)
		@correct_bitscores = correct_bitscores
	end
	def add_incorrect_bitscores(incorrect_bitscores)
		@incorrect_bitscores = incorrect_bitscores
	end
end


# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
genes_dump_file = ARGV[2]
minimum_number_of_exons_per_gene = ARGV[3].to_i

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Collect names of nucleotide fasta files in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in)
filenames_in = []
dir_entries_in.each {|e| filenames_in << e if e.match(/.*_nucl.fasta/)}

# Get exon names of alignments from the filenames.
exon_alignment_names = []
filenames_in.each do |f|
	exon_alignment_names << f.chomp("_nucl.fasta")
end

# Load the genes dump file.
genes = Marshal.load(File.read(genes_dump_file))

# Create an array to store the number of missing sequences per species.
exons_alignments_per_gene = []

# For each gene, test whether at least three exons alignments are still present.
genes.each do |g|
	exon_alignments_for_this_gene = []
	g.transcripts.each do |t|
		t.exons.each do |e|
			exon_alignments_for_this_gene << e.id if exon_alignment_names.include?(e.id)
		end
	end
	# If the minimum number of exon alignments is still present, convert the fasta
	# alignment files to phylip format and save these in a directory inside the
	# alignments output directory.
	if exon_alignments_for_this_gene.uniq.size >= minimum_number_of_exons_per_gene
		Dir.mkdir("#{alignment_directory_out}/#{g.id}")
		exon_alignments_for_this_gene.uniq.each do |e|
			fasta_file_name = "#{e}_nucl.fasta"
			fasta_file = File.open("#{alignment_directory_in}/#{fasta_file_name}")
			fasta_lines = fasta_file.readlines
			fasta_ids = []
			fasta_seqs = []
			fasta_lines.each do |l|
				if l[0] == ">"
					fasta_ids << l[1..-1].strip
					fasta_seqs << ""
				elsif l.strip != ""
					fasta_seqs.last << l.strip
				end
			end
			nexus_string = "#nexus\n"
			nexus_string << "\n"
			nexus_string << "begin data;\n"
			nexus_string << "dimensions  ntax=#{fasta_ids.size} nchar=#{fasta_seqs[0].size};\n"
			nexus_string << "format datatype=DNA gap=- missing=?;\n"
			nexus_string << "matrix\n"
			fasta_ids.size.times do |x|
				nexus_string << "#{fasta_ids[x].ljust(12)}#{fasta_seqs[x]}\n"
			end
			nexus_string << ";\n"
			nexus_string << "end;\n"
			nexus_file_name = "#{e}.seq"
			nexus_file = File.open("#{alignment_directory_out}/#{g.id}/#{nexus_file_name}","w")
			nexus_file.write(nexus_string)
		end
	end
end
