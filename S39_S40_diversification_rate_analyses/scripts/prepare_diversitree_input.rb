# Michael Matschiner, 2015-09-17

class Clade
	attr_reader :name, :species, :ntax, :age, :traits
	def initialize(name,ntax,traits)
		@name = name
		@ntax = ntax
		@traits = traits
		@species = []
		@speciation_rates = []
		@extinction_rates = []
	end
	def add_age(age)
		@age = age
	end
	def large
		if @ntax >= 100
			true
		else
			false
		end
	end
	def mean_speciation
		speciation_rates_sum = 0
		@speciation_rates.each do |r|
			speciation_rates_sum += r
		end
		speciation_rates_sum/@speciation_rates.size.to_f
	end
	def mean_extinction
		extinction_rates_sum = 0
		@extinction_rates.each do |r|
			extinction_rates_sum += r
		end
		extinction_rates_sum/@extinction_rates.size.to_f
	end
	def add_species(species)
		@species << species
	end
	def add_speciation_rate(speciation_rate)
		@speciation_rates << speciation_rate
	end
	def add_extinction_rate(extinction_rate)
		@extinction_rates << extinction_rate
	end
	def info
		info_string = "#{@name.to_s.ljust(20)}#{@age.to_s.rjust(10)}#{@ntax.to_s.rjust(5)}#{self.mean_speciation.round(3).to_s.rjust(10)}#{self.mean_extinction.round(3).to_s.rjust(10)}#{@traits.size.to_s.rjust(10)}#{@species.size.to_s.rjust(10)}"
		info_string
	end
end

# Load phylsim.
$libPath = "./resources/phylsim/"
require "./resources/phylsim/main.rb"

# Load FileUtils.
require 'fileutils'

# Get the command line arguments.
input_tree_file_name = ARGV[0]
input_table_file_name = ARGV[1]
analysis_dir_prefix = ARGV[2]
analysis_dir_prefix = analysis_dir_prefix.chomp("/")
species_table_file_name = ARGV[3]
bamm_branch_file_name = ARGV[4]
start_file_name = ARGV[5]

# Initiate an array to hold all clades.
clades = []

# Read the input table file.
input_table_file = File.open(input_table_file_name)
input_table_lines = input_table_file.readlines
input_table_lines[1..-1].each do |l|
	line_ary = l.split
	name = line_ary[0]
	ntax = line_ary[1].to_i
	traits = line_ary[2..-1]
	clades << Clade.new(name,ntax,traits)
end

# Read the tree as a string.
input_tree_file = File.open(input_tree_file_name)
tree_string = input_tree_file.read

# Read the age of each clade from the string and add it to the clade objects.
clades.each do |c|
	tree_string_fragment = tree_string.match(/#{c.name}:\d+\.\d+/)
	clade_age = tree_string_fragment.to_s.split(":")[1].to_f
	c.add_age(clade_age)
end

# Read the species table.
species_table_file = File.open(species_table_file_name)
species_table_lines = species_table_file.readlines
species_table_clades = []
species_table_species = []
species_table_lines[1..-1].each do |l|
	line_ary = l.split
	species_table_species << line_ary[0]
	species_table_clades << line_ary[1]
end

# Assign species to all clades.
clades.each do |c|
	species_table_clades.size.times do |x|
		if c.name == species_table_clades[x]
			c.add_species(species_table_species[x])
		end
	end
end

# Read the file with net diversification estimates per branch.
bamm_branch_file = File.open(bamm_branch_file_name)
bamm_branch_lines = bamm_branch_file.readlines
branch_begins = []
branch_ends = []
branch_speciation = []
branch_extinction = []
branch_tip_labels = []
tip_labels = []
in_branch_begins = false
in_branch_ends = false
in_branch_speciation = false
in_branch_extinction = false
in_branch_tip_labels = false
bamm_branch_lines.each do |l|
	if l.strip == "Branch begin"
		in_branch_begins = true
	elsif l.strip == "Branch end"
		in_branch_ends = true
		in_branch_begins = false
	elsif l.strip == "Branch speciation"
		in_branch_speciation = true
		in_branch_ends = false
	elsif l.strip == "Branch extinction"
		in_branch_extinction = true
		in_branch_speciation = false
	elsif l.strip == "Branch tip label"
		in_branch_tip_labels = true
		in_branch_extinction = false
	elsif in_branch_begins
		branch_begins << l.strip.to_f unless l.strip == ""
	elsif in_branch_ends
		branch_ends << l.strip.to_f unless l.strip == ""
	elsif in_branch_speciation
		branch_speciation << l.strip.to_f unless l.strip == ""
	elsif in_branch_extinction
		branch_extinction << l.strip.to_f unless l.strip == ""
	elsif in_branch_tip_labels
		tip_labels << l.strip unless l.strip == ""
	end
end

# Determine maximum branch age.
branch_max_age = 0
branch_ends.each do |e|
	branch_max_age = e if e > branch_max_age
end

# Fill array branch_tip_labels according to branch_begins, branch_ends, and branch_net_divs.
branch_ends.size.times do |x|
	if (branch_max_age-branch_ends[x]).round(2) == 0
		branch_tip_labels << tip_labels.shift
	else
		branch_tip_labels << "unknown"
	end
end

# Add speciation and extinction rates to all clades.
clades.each do |c|
	c.species.each do |s|
		branch_tip_labels.size.times do |x|
			if branch_tip_labels[x] == s
				c.add_speciation_rate(branch_speciation[x])
				c.add_extinction_rate(branch_extinction[x])
			end
		end
	end
end

# Repeat multiple times to produce replicates.
number_of_replicates = 25
number_of_replicates.times do |rep|

	# Randomly resolve large clades into subclades.
	clades_to_be_removed = []
	new_clades = []
	clades.each do |c|
		if c.large
			min_lambda = 0.8 * c.mean_speciation
			max_lambda = 1.2 * c.mean_speciation
			min_mu = 0.8 * c.mean_extinction
			max_mu = 1.2 * c.mean_extinction
			treeOrigin = c.age
			min_np = (0.99 * c.ntax).to_i
			max_np = (1.01 * c.ntax).to_i
			tree = Tree.generate([min_lambda,max_lambda], [min_mu,max_mu], treeOrigin, 0, 0, false, np = [min_np,max_np], npEach = [0,'inf'], checkProbabilities = false, algorithm = "forward", verbose = true, threads = 1)
			tree.reconstruct("diversified",c.ntax/20,nil,false)
			clade_tree_age = tree.treeOrigin
			clade_tree_string = tree.to_newick(nil,"duration",false,false,true,false,false)
			clade_tree_taxon_names = []
			clade_tree_taxon_diversities = []
			tree.branch.each do |b|
				if b.extant
					if b.speciesId.include?("/")
						clade_tree_taxon_names << b.speciesId.split("/").last
					else
						clade_tree_taxon_names << b.speciesId
					end
					clade_tree_taxon_diversities << b.originalExtantDiversity
				end
			end
			clade_tree_taxon_names.size.times do |x|
				new_taxon_name = "#{c.name}_#{(x+1).to_s.rjust(3).gsub(" ","0")}"
				clade_tree_string.sub!(/#{clade_tree_taxon_names[x]}:/,"#{new_taxon_name}:")
				new_clades << Clade.new(new_taxon_name,clade_tree_taxon_diversities[x],c.traits)
			end
			tree_string.sub!(/#{c.name}:\d+\.\d+/,"#{clade_tree_string}:#{c.age-clade_tree_age}")
			clades_to_be_removed << c.name
		else
			new_clades << c
		end
	end

	# Make directories for replicates with different thresholds for high/low trait values.
	26.times do |t|
		threshold = 2 * t + 10
		analysis_dir_name = "#{analysis_dir_prefix}/t#{threshold.to_s.rjust(2).sub(" ","0")}_r#{(rep+1).to_s.rjust(3).gsub(" ","0")}"
		FileUtils.mkdir_p(analysis_dir_name)
		output_tree_file = File.open("#{analysis_dir_name}/unresolved.tre","w")
		output_tree_file.write(tree_string)
		output_tree_file.close
		resolved_vector_string = "tip.label\tstate\n"
		unresolved_table_string = "tip.label\tNc\tn0\tn1\n"
		new_clades.each do |c|
			if c.ntax == 1
				if c.traits[0].to_f < threshold
					resolved_vector_string << "#{c.name}\t0\n"
				else
					resolved_vector_string << "#{c.name}\t1\n"
				end
			else
				resolved_vector_string << "#{c.name}\tNA\n"
				unresolved_table_string << "#{c.name}\t#{c.ntax}"
				n0 = 0
				n1 = 0
				c.traits.each do |i|
					if i.to_f < threshold
						n0 += 1
					else
						n1 += 1
					end
				end
				n0_all = 0
				n1_all = 0
				c.ntax.times do
					if rand < n0/(n0+n1).to_f
						n0_all += 1
					else
						n1_all += 1
					end
				end
				unresolved_table_string << "\t#{n0_all}\t#{n1_all}\n"
			end
		end
		unresolved_table_file = File.open("#{analysis_dir_name}/unresolved.txt","w")
		unresolved_table_file.write(unresolved_table_string)
		unresolved_table_file.close
		resolved_vector_file = File.open("#{analysis_dir_name}/resolved.txt","w")
		resolved_vector_file.write(resolved_vector_string)
		resolved_vector_file.close

		# Copy the start file to the analysis directory.
		FileUtils.cp(start_file_name,"#{analysis_dir_name}/start.r")
	end

# end
