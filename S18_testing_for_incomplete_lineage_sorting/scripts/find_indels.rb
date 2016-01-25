# Michael Matschiner, 2015-06-15

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

class Indel
	attr_reader :gene_id, :species, :from, :to, :overlaps, :species_with_this_indel, :species_without_this_indel
	def initialize(gene_id,species,from,to)
		@gene_id = gene_id
		@species = [species]
		@from = from
		@to = to
		@overlaps = false
		@species_with_this_indel = nil
		@species_without_this_indel = nil
	end
	def length
		@to-@from+1
	end
	def set_overlaps(overlaps)
		@overlaps = overlaps
	end
	def add_species(spc)
		@species << spc
	end
	def add_species_with_this_indel(species_with_this_indel)
		@species_with_this_indel = species_with_this_indel
	end
	def add_species_without_this_indel(species_without_this_indel)
		@species_without_this_indel = species_without_this_indel
	end
	def to_s
		outstring = ""
		outstring << "#{@gene_id}#{@from.to_s.rjust(8)}#{@to.to_s.rjust(8)}    "
		@species_with_this_indel.each {|s| outstring << "#{s},"}
		outstring.chomp!(",")
		outstring << "   "
		@species_without_this_indel.each {|s| outstring << "#{s},"}
		outstring.chomp!(",")
		outstring << "\n"
		outstring
	end
end

# Get the command line arguments.
fasta_dir = ARGV[0]
output_file_name = ARGV[1]
minimum_number_of_ambiguous_species = ARGV[2].to_i
tree_directory = ARGV[3] # The tree directory is only needed to write PAUP* files including a gettrees command.
paup_directory = ARGV[4] # The PAUP* directory is only needed to be included in the PAUP* block written.

# Get all tree file names in the tree directory.
tree_file_names = []
Dir.entries(tree_directory).each {|e| tree_file_names << e if e.match(/.tre/)}

# Make the paup directory if it does not exist yet.
FileUtils.mkdir_p(paup_directory)

# Make the output directory if it does not exist yet.
if output_file_name.include?("/")
	output_dir = output_file_name.chomp(output_file_name.split("/").last).chomp("/")
	# Create the output directories if they don't exist yet.
	FileUtils.mkdir_p(output_dir)
end

# Read the fasta file names.
fasta_file_names = []
Dir.entries(fasta_dir).each {|e| fasta_file_names << e if e.match(/_nucl.fasta/)}

# Read the first fasta file to get the taxon ids.
fasta_file = File.open("#{fasta_dir}/#{fasta_file_names[0]}")
fasta_lines = fasta_file.readlines
fasta_file.close
taxon_ids = []
fasta_lines.each do |l|
	if l[0] == ">"
		taxon_ids << l[1..-1].strip
	end
end

# Initiate the output string.
output_string = ""

# For each fasta file, find indels.
indels = []
fasta_file_names.each do |n|

	# Read the nexus file.
	gene_id = n.chomp("_nucl.fasta")
	fasta_file = File.open("#{fasta_dir}/#{n}")
	fasta_lines = fasta_file.readlines
	fasta_file.close
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

	# Make sure the fasta ids are identical to the taxon ids.
	raise "ERROR: Taxon ids in fasta files differ!" unless taxon_ids == fasta_ids

	# Find indels in sequences.
	indels_of_this_gene = []
	fasta_seqs.size.times do |x|
		1.upto(fasta_seqs[x].size-1) do |pos|
			if fasta_seqs[x][pos] == "-" and fasta_seqs[x][pos-1] != "-"
				continue = true unless pos == fasta_seqs[x].size-1
				add_to_pos = 0
				while continue
					add_to_pos += 1
					if pos+add_to_pos == fasta_seqs[x].size-1
						continue = false
					else
						if fasta_seqs[x][pos+add_to_pos] != "-" and fasta_seqs[x][pos-1+add_to_pos] == "-"
							continue = false
							indels_of_this_gene << Indel.new(gene_id,fasta_ids[x],pos,pos-1+add_to_pos)
						end
					end
				end
			end
		end
	end

	# Remove indels with lengths that are not multiples of three.
	indels_of_this_gene.size.times do |x|
		indels_of_this_gene[x] = nil unless (indels_of_this_gene[x].length/3)*3 == indels_of_this_gene[x].length
	end
	indels_of_this_gene.compact!

	# Remove overlapping indels.
	0.upto(indels_of_this_gene.size-2) do |x|
		(x+1).upto(indels_of_this_gene.size-1) do |y|
			overlaps = false
			overlaps = true if indels_of_this_gene[x].from < indels_of_this_gene[y].from and indels_of_this_gene[x].to >= indels_of_this_gene[y].from-1
			overlaps = true if indels_of_this_gene[x].to > indels_of_this_gene[y].to and indels_of_this_gene[x].from <= indels_of_this_gene[y].to+1
			overlaps = true if indels_of_this_gene[y].from < indels_of_this_gene[x].from and indels_of_this_gene[y].to >= indels_of_this_gene[x].from-1
			overlaps = true if indels_of_this_gene[y].to > indels_of_this_gene[x].to and indels_of_this_gene[y].from <= indels_of_this_gene[x].to+1
			if overlaps
				indels_of_this_gene[x].set_overlaps(overlaps)
				indels_of_this_gene[y].set_overlaps(overlaps)
			end
		end
	end
	indels_of_this_gene.size.times do |x|
		indels_of_this_gene[x] = nil if indels_of_this_gene[x].overlaps
	end
	indels_of_this_gene.compact!

	# If indels have identical start and end positions in multiple species, combine them to a single indel with multiple species.
	0.upto(indels_of_this_gene.size-2) do |x|
		unless indels_of_this_gene[x] == nil
			(x+1).upto(indels_of_this_gene.size-1) do |y|
				unless indels_of_this_gene[y] == nil
					if indels_of_this_gene[x].from == indels_of_this_gene[y].from and indels_of_this_gene[x].to == indels_of_this_gene[y].to
						indels_of_this_gene[y].species.each {|s| indels_of_this_gene[x].add_species(s)}
						indels_of_this_gene[y] = nil
					end
				end
			end
		end
	end
	indels_of_this_gene.compact!

	# For each indel of this gene, find out which species share it, which ones don't, and which are ambiguous due to missing data.
	indels_of_this_gene.each do |i|
		species_with_indel = []
		species_without_indel = []
		ambiguous_species = []
		fasta_ids.size.times do |x|
			if fasta_seqs[x][i.from] == "-" and fasta_seqs[x][i.from-1] != "-" and fasta_seqs[x][i.to] == "-" and fasta_seqs[x][i.to+1] != "-"
				species_with_indel << fasta_ids[x]
			elsif fasta_seqs[x][i.from] != "-" and fasta_seqs[x][i.from-1] != "-" and fasta_seqs[x][i.to] != "-" and fasta_seqs[x][i.to+1] != "-"
				species_without_indel << fasta_ids[x]
			elsif fasta_seqs[x].match("-+$") or fasta_seqs[x].match("^-+")
				ambiguous_species << fasta_ids[x]
			else
				raise "ERROR: Unexpected indel pattern in gene #{gene_id}, pos #{i.from} - #{i.to}, for species #{fasta_ids[x]}!"
			end
		end
		i.add_species_with_this_indel(species_with_indel)
		i.add_species_without_this_indel(species_without_indel)
	end

	# Remove all indels that don't have at least one species with and one species without this indel.
	indels_of_this_gene.size.times do |x|
		if indels_of_this_gene[x].species_with_this_indel.size == 0 or indels_of_this_gene[x].species_without_this_indel.size == 0
			indels_of_this_gene[x] = nil
		end
	end
	indels_of_this_gene.compact!

	# Remove all indels with more than the minimum number of ambiguous species.
	indels_of_this_gene.size.times do |x|
		if fasta_ids.size - (indels_of_this_gene[x].species_with_this_indel.size + indels_of_this_gene[x].species_without_this_indel.size) > minimum_number_of_ambiguous_species
			indels_of_this_gene[x] = nil
		end
	end
	indels_of_this_gene.compact!

	# Add to the array of all indels.
	indels_of_this_gene.each do |i|
		indels << i
	end

	# Add to output string.
	indels_of_this_gene.each do |i|
		output_string << i.to_s
	end

end

# Write the indel list.
output_file = File.open(output_file_name,"w")
output_file.write(output_string)
output_file.close

# Prepare a nexus file with indel information.
nexus_string = "#nexus\n"
nexus_string << "begin data;\n"
nexus_string << "  dimensions  ntax=#{taxon_ids.size} nchar=#{indels.size};\n"
nexus_string << "  format datatype=Standard gap=- missing=?;\n"
nexus_string << "  matrix\n"
taxon_ids.each do |t|
	nexus_string << "  #{t.ljust(10)}"
	indels.each do |i|
		if i.species_with_this_indel.include?(t)
			nexus_string << "1"
		elsif i.species_without_this_indel.include?(t)
			nexus_string << "0"
		else
			nexus_string << "?"
		end
	end
	nexus_string << "\n"
end
nexus_string << "  ;\n"
nexus_string << "end;\n"

# Write the nexus file.
nexus_file_name = output_file_name.chomp(".txt") + ".nex"
nexus_file = File.open(nexus_file_name,"w")
nexus_file.write(nexus_string)
nexus_file.close

# Prepare a PAUP* block.
paup_block_string = "#nexus\n"
paup_block_string << "begin paup;\n"
tree_file_names.each do |t|
	paup_block_string << "  log start file='../#{paup_directory}/#{t.chomp(".tre")}.log' replace=yes;\n"
	paup_block_string << "  gettrees file='../#{tree_directory}/#{t}';\n"
	# paup_block_string << "  outgroup Astmex Danrer;\n"
	# paup_block_string << "  root;\n"
	paup_block_string << "  describetrees / brlens=yes;\n"
	paup_block_string << "  execute ../#{nexus_file_name};\n"
	paup_block_string << "  gettrees file='../#{tree_directory}/#{t}';\n"
	paup_block_string << "  outgroup Astmex Danrer;\n"
	paup_block_string << "  root;\n"
	paup_block_string << "  describetrees / brlens=yes apolist;\n"
	paup_block_string << "  log stop;\n"
end
paup_block_string << "  quit;\n"
paup_block_string << "end;\n"

# Write the PAUP* block.
paup_block_file_name = output_file_name.chomp(".txt") + ".paup"
paup_block_file = File.open(paup_block_file_name,"w")
paup_block_file.write(paup_block_string)
paup_block_file.close
