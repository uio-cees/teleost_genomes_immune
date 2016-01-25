# This script reads the table of zebrafish exons from Biomart
# and filters them as follows:
#
# 1. If multiple transcripts are known for a gene, only exons of
#    the transcript with the longest overall length are used.
# 2. Exons shorter than 300 bp are discarded.
# 3. Exons are discarded if their gene has less than three exons
#    that passed filters 1 and 2.
#
# Exons that pass these filters are translated and written to a
# fasta file.
# Note that this file will work only with ENSEMBL v.78 and with
# the ten teleost genomes of this ENSEMBL version. If more or less
# genomes should be used, multiple parts of this script need to be
# changes. Search for 'ENS' to find those parts.

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

# Definitions.
subject_database_path = "../data"
subject_databases = []
subject_databases << "Astmex.cdna.fasta"
subject_databases << "Danrer.cdna.fasta"
subject_databases << "Gadmor.cdna.fasta"
subject_databases << "Gasacu.cdna.fasta"
subject_databases << "Orenil.cdna.fasta"
subject_databases << "Orylat.cdna.fasta"
subject_databases << "Poefor.cdna.fasta"
subject_databases << "Takrub.cdna.fasta"
subject_databases << "Tetnig.cdna.fasta"
subject_databases << "Xipmac.cdna.fasta"
query_file_name = "tmp_query.txt"
subject_tags = []
subject_tags << "ENSAMX"
subject_tags << "ENSDAR"
subject_tags << "ENSGMO"
subject_tags << "ENSGAC"
subject_tags << "ENSONI"
subject_tags << "ENSORL"
subject_tags << "ENSPFO"
subject_tags << "ENSTRU"
subject_tags << "ENSTNI"
subject_tags << "ENSXMA"
# The next two parameters specify the factor by which a bitscore must be
# better than the worst incorrect bitscore in order to count as an ortholog,
# and the number of taxa for which orthologs must be identified in this way,
# otherwise, the exon will be removed.
ortholog_bitscore_ratio_threshold = 1.5
number_of_required_orthologs = 8
ortholog_minimum_bitscore = 50

# See if a dump file exists already.
dump_file_name = "../data/genes.dmp"
if File.exists?(dump_file_name)

	# Read the dump file.
	print "Reading dump file..."
	dump_file = File.open(dump_file_name)
	genes = Marshal.load(dump_file.read)
	puts " done."

else

	# Read the biomart table.
	biomart_table_file_name = "../data/exons_from_biomart.txt"
	print "Reading file #{biomart_table_file_name}..."
	biomart_table_file = File.open(biomart_table_file_name)
	biomart_table = biomart_table_file.readlines
	biomart_table = biomart_table[1..-1].sort!
	puts " done."

	# Read the genome sequence fasta file.
	fasta_file_name = "../data/Danrer.fasta"
	print "Reading file #{fasta_file_name}..."
	fasta_file = File.open(fasta_file_name)
	fasta_lines = fasta_file.readlines
	fasta_ids = []
	fasta_seqs = []
	read_this_seq = false
	fasta_lines.each do |l|
		if l[0] == ">"
			if l.match(/>(\d+) dna:chromosome/)
				read_this_seq = true
				fasta_ids << $1.to_i
				fasta_seqs << ""
			else
				read_this_seq = false
			end
		elsif l.strip != "" and read_this_seq
			fasta_seqs.last << l.strip
		end
	end
	puts " done."

	# Analyse the biomart table.
	print "Analysing the biomart table..."
	genes = []
	gene_id = ""
	transcript_id = ""
	biomart_table.each do |l|
		last_gene_id = gene_id
		last_transcript_id = transcript_id
		row_ary = l.split
		gene_id = row_ary[0]
		transcript_id = row_ary[1]
		chromosome = row_ary[2].to_i
		if chromosome > 0
			transcript_start = row_ary[3].to_i
			transcript_end = row_ary[4].to_i
			strand = row_ary[5].to_i
			exon_id = row_ary[6]
			exon_start = row_ary[7].to_i
			exon_end = row_ary[8].to_i
			phase = row_ary[9].to_i
			exon = Exon.new(exon_id)
			exon.add_start(exon_start)
			exon.add_end(exon_end)
			exon.add_phase(phase)
			if transcript_id != last_transcript_id
				transcript = Transcript.new(transcript_id)
				transcript.add_start(transcript_start)
				transcript.add_end(transcript_end)
				transcript.add_strand(strand)
				transcript.add_exon(exon)
				if gene_id != last_gene_id
					gene = Gene.new(gene_id)
					gene.add_chromosome(chromosome)
					genes << gene
				end
				genes.last.add_transcript(transcript)
			else
				genes.last.transcripts.last.add_exon(exon)
			end
		end
	end
	puts " done."

	# Sort all exons per transcript.
	print "Sorting exons for each transcript..."
	genes.each do |g|
		g.transcripts.each do |t|
			t.sort_exons
		end
	end
	puts " done."

	# Select exons that pass filters 1, 2, and 3, as well as their transcripts and genes.
	print "Marking filtered genes, transcripts, and exons..."
	genes.each do |g|
		long_exons = g.longest_transcript.long_exons
		if long_exons.size > 2
			long_exons.each do |e|
				e.select if e.phase != -1
			end
		end
	end
	puts " done."

	# Add the sequence and translation to each exon, and deselect it if the translation contains
	# stop codons before the last codon.
	print "Adding the sequence and sequence tranlation of each marked exon..."
	genes.each do |g|
		if g.selected
			chromosome_index = fasta_ids.index(g.chromosome)
			g.transcripts.each do |t|
				if t.selected
					t.exons.each do |e|
						if e.selected
							if t.strand == 1
								exon_seq = fasta_seqs[chromosome_index][e.start-1..e.end-1]
							elsif t.strand == -1
								exon_seq = fasta_seqs[chromosome_index][e.start-1..e.end-1].reverse_complement
							else
								raise "Unexpected transcript strand!"
							end
							if e.phase == 1
								exon_seq = exon_seq[2..-1]
							elsif e.phase == 2
								exon_seq = exon_seq[1..-1]
							elsif e.phase != 0
								raise "Unexpected phase!"
							end
							exon_seq = exon_seq[0..((exon_seq.size/3)*3)-1]
							e.add_seq(exon_seq)
							e.translate
							if e.translation[0..-2].include?("*") or e.translation.size < 50
								e.deselect
							end
						end
					end
				end
			end
		end
	end
	puts " done."

	# Read the ENSEMBL gene trees file to memory.
	gene_tree_file_name = "../data/Compara.78.protein.nhx.emf"
	print "Reading the ENSEMBL gene tree file..."
	gene_tree_file = File.open(gene_tree_file_name)
	gene_trees = []
	previous_line = ""
	while line = gene_tree_file.gets do
		if previous_line[0..3] == "DATA"
			gene_trees << line
		end
		previous_line = line
	end
	puts " done."

	# For each selected gene, check the gene tree.
	print "Checking for gene id presence in gene trees..."
	genes.each do |g|
		if g.selected
			gene_trees.each do |gt|
				if gt.include?(g.id)
					g.add_gene_tree(gt)
				end
			end
		end
	end
	puts " done."

	# Deselect exons of genes without gene trees.
	print "Deselecting exons of genes without gene tres..."
	genes.each do |g|
		g.transcripts.each do |t|
			t.exons.each do |e|
				e.deselect if g.gene_tree == nil
			end
		end
	end
	puts " done."

	# Check whether the gene tree for each gene conforms to expectations.
	print "Checking gene trees and collecting ortholog ids..."
	genes.each do |g|
		use_gene = false
		use_gene = true if g.selected
		if use_gene
			gene_id_length = g.id.length
			index = gene_id_length
			g.gene_tree.size.times do |x|
				if g.id == g.gene_tree[x..x+gene_id_length-1]
					index = x
					break
				end
			end
			num_closing_brackets = 0
			forward_index = index + gene_id_length - 1
			until num_closing_brackets == 1
				forward_index += 1
				if g.gene_tree[forward_index] == ")"
					num_closing_brackets += 1
				elsif g.gene_tree[forward_index] == "("
					num_closing_brackets -= 1
				end
			end
			num_opening_brackets = 0
			backwards_index = index
			until num_opening_brackets == 1
				backwards_index -= 1
				if g.gene_tree[backwards_index] == "("
					num_opening_brackets += 1
				elsif g.gene_tree[backwards_index] == ")"
					num_opening_brackets -= 1
				end
			end
			subject_tags.each do |st|
				if st == "ENSAMX" or st == "ENSDAR"
					use_gene = false if g.gene_tree[backwards_index..forward_index].scan("#{st}P").size != 1
				else
					use_gene = false if g.gene_tree[backwards_index..forward_index].scan("#{st}P").size != 0
				end
			end
		end
		if use_gene
			until num_closing_brackets == 2
				forward_index += 1
				if g.gene_tree[forward_index] == ")"
					num_closing_brackets += 1
				elsif g.gene_tree[forward_index] == "("
					num_closing_brackets -= 1
				end
			end
			until num_opening_brackets == 2
				backwards_index -= 1
				if g.gene_tree[backwards_index] == "("
					num_opening_brackets += 1
				elsif g.gene_tree[backwards_index] == ")"
					num_opening_brackets -= 1
				end
			end
			subject_tags.each do |st|
				use_gene = false if g.gene_tree[backwards_index..forward_index].scan("#{st}P").size != 1
			end
			# Filter those trees, where a non-teleost species is nested within the ten teleosts in the
			# ENSEMBL gene tree.
			# The number of occurrences of 'ENS' should be 20 cause each species, and with ENSEMBL v. 78
			# we have 10 teleost species, is listed with gene and protein id.
			use_gene = false if use_gene and g.gene_tree[backwards_index..forward_index].scan("ENS").size != 2*subject_tags.size
			# Filter trees that include an annotated duplication node.
			use_gene = false if use_gene and g.gene_tree[backwards_index..forward_index].include?(":D=Y")
		end
		if use_gene
			ext_forward_index = forward_index
			ext_backwards_index = backwards_index
			until num_closing_brackets == 3
				ext_forward_index += 1
				if g.gene_tree[ext_forward_index] == ")"
					num_closing_brackets += 1
				elsif g.gene_tree[ext_forward_index] == "("
					num_closing_brackets -= 1
				end
			end
			until num_opening_brackets == 3
				ext_backwards_index -= 1
				if g.gene_tree[ext_backwards_index] == "("
					num_opening_brackets += 1
				elsif g.gene_tree[ext_backwards_index] == ")"
					num_opening_brackets -= 1
				end
			end
			subject_tags.each do |st|
				use_gene = false if g.gene_tree[ext_backwards_index..ext_forward_index].scan("#{st}P").size != 1
			end
		end
		if use_gene
			ortholog_ids = []
			subject_tags.each do |st|
				ortholog_ids << g.gene_tree[backwards_index..forward_index].match(/#{st}G\d+/).to_s
			end
			# Make sure exactly 10 orthologs are included in the string (this only works with the ten
			# teleosts of ENSEMBL v.78).
			if ortholog_ids.size == subject_tags.size
				ortholog_ids.each do |o|
					g.add_ortholog_id(o)
				end
			else
				raise "Wrong number of ortholog ids found: #{ortholog_ids.size}!"
			end
		end
		unless use_gene
			g.transcripts.each do |t|
				t.exons.each do |e|
					e.deselect
				end
			end
		end
	end
	puts " done."

	# Get the worst tblastn score for each selected exon.
	print "Running tblastn to determine the worst score for each exon..."
	genes.each do |g|
		if g.selected
			g.transcripts.each do |t|
				if t.selected
					t.exons.each do |e|
						if e.selected
							# Write a query file with the amino acid sequence of
							# the current exon.
							query_file = File.new(query_file_name,"w")
							query_file.write(">query\n#{e.translation}\n")
							query_file.close
							orthologs = []
							correct_bitscores = []
							incorrect_bitscores = []
							# For each database, run tblastn searches.
							subject_databases.size.times do |x|
								db = subject_databases[x]
								hits_string = `tblastn -db #{subject_database_path}/#{db} -max_target_seqs 20 -query #{query_file_name} -outfmt '6 stitle bitscore'`
								hits = hits_string.split("\n")
								correct_bitscores_for_this_subject = []
								incorrect_bitscores_for_this_subject = []
								hits.each do |h|
									stitle = h.split("\t")[0]
									bitscore = h.split("\t")[1].to_f
									if stitle.include?(g.ortholog_ids[x])
										correct_bitscores_for_this_subject << bitscore
									else
										incorrect_bitscores_for_this_subject << bitscore
									end
								end
								if correct_bitscores_for_this_subject == []
									correct_bitscores << 0
								else
									correct_bitscores << correct_bitscores_for_this_subject.max
								end
								if incorrect_bitscores_for_this_subject == []
									incorrect_bitscores << 0
								else
									incorrect_bitscores << incorrect_bitscores_for_this_subject.max
								end
								orthologs << g.ortholog_ids[x]
							end
							e.add_correct_bitscores(correct_bitscores)
							e.add_incorrect_bitscores(incorrect_bitscores)
						end
					end
				end
			end
		end
	end

	# Write the genes to a dump file.
	print "Writing dump file..."
	dump_file = File.open(dump_file_name,"w")
	dump_file.write(Marshal.dump(genes))
	puts " done."

end

# Deselect exons according to blast bitscores.
print "Deselecting exons according to blast bitscores..."
genes.each do |g|
	if g.selected
		g.transcripts.each do |t|
			if t.selected
				t.exons.each do |e|
					if e.selected
						number_of_orthologs = 0
						e.correct_bitscores.each do |bs|
							if bs > ortholog_bitscore_ratio_threshold * e.incorrect_bitscores.max
								if bs > ortholog_minimum_bitscore
									number_of_orthologs += 1
								end
							end
						end
						if number_of_orthologs < number_of_required_orthologs
							e.deselect
						end
					end
				end
			end
		end
	end
end
puts " done."

# Prepare the output table and fasta string.
print "Preparing output strings..."
gene_count = 0
transcript_count = 0
exon_count = 0
gene_table_string = "danrer_gene_id\tdanrer_transcript_id\tnumber_of_exons\tastmex_gene_id\tdanrer_gene_id\tgadmor_gene_id\tgasacu_gene_id\torenil_gene_id\torylat_gene_id\tpoefor\ttakrub_gene_id\ttetnig_gene_id\txipmac_gene_id\n"
exon_table_string = "exon_id\tgene_id\ttranscript_id\ttranslation\tlength\tthreshold_bitscore\tcorrect_bitscores\tincorrect_bitscores\tmin(correct_bitscores)-max(incorrect_bitscores)\n"
exon_fasta_string = ""
genes.each do |g|
	if g.selected
		exon_count_for_this_gene = 0
		gene_count += 1
		g.transcripts.each do |t|
			if t.selected
				transcript_count += 1
				t.exons.each do |e|
					if e.selected
						exon_count_for_this_gene += 1
						ortholog_threshold_bitscore = [(ortholog_bitscore_ratio_threshold * e.incorrect_bitscores.max),ortholog_minimum_bitscore].max.round(2)
						exon_table_string << "#{e.id}\t#{g.id}\t#{t.id}\t#{e.translation[0..9]}..#{e.translation[-10..-1]}\t#{e.length}\t#{ortholog_threshold_bitscore}\t["
						e.correct_bitscores.each do |bs|
							exon_table_string << "#{bs},"
						end
						exon_table_string.chomp!(",")
						exon_table_string << "]\t["
						e.incorrect_bitscores.each do |bs|
							exon_table_string << "#{bs},"
						end
						exon_table_string.chomp!(",")
						exon_table_string << "]\t#{(e.correct_bitscores.min-e.incorrect_bitscores.max).round(1)}\n"
						exon_fasta_string << ">#{e.id}[&bitscore=#{ortholog_threshold_bitscore}]\n#{e.translation}\n\n"
						exon_count += 1
					end
				end
			end
		end
		gene_table_string << "#{g.id}\t#{g.transcripts[0].id}\t#{exon_count_for_this_gene}\t"
		g.ortholog_ids.each do |o|
			gene_table_string << "#{o}\t"
		end
		gene_table_string.chomp("\t")
		gene_table_string << "\n"
	end
end
puts " done."
puts "#{gene_count} genes, #{transcript_count} transcripts, and #{exon_count} exons remain after filtering."

# Write the gene output table.
gene_table_file_name = "../analysis/nuclear_queries_genes.txt"
gene_table_file = File.open(gene_table_file_name,"w")
gene_table_file.write(gene_table_string)
puts "Wrote file #{gene_table_file_name}."

# Write the exon output table.
exon_table_file_name = "../analysis/nuclear_queries_exons.txt"
exon_table_file = File.open(exon_table_file_name,"w")
exon_table_file.write(exon_table_string)
puts "Wrote file #{exon_table_file_name}."

# Write the exon output fasta file.
exon_fasta_file_name = "../analysis/exons.fasta"
exon_fasta_file = File.open(exon_fasta_file_name,"w")
exon_fasta_file.write(exon_fasta_string)
puts "Wrote file #{exon_fasta_file_name}."
