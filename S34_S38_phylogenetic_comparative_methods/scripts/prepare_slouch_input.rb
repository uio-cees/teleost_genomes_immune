# Michael Matschiner, 2015-07-16

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Load the phylsim library.
$libPath = "./resources/phylsim/"
require "#{$libPath}/main.rb"

class TableLine
	attr_reader :branch_id, :number, :regime, :species_id
	def initialize(number, branch_id, species_id, time)
		@number = number
		@branch_id = branch_id
		@species_id = species_id
		@time = time
		@ancestor = nil
		@regime = nil
		@trait = nil
		@trait_stderr = nil
	end
	def set_ancestor(ancestor)
		@ancestor = ancestor
	end
	def set_regime(regime)
		@regime = regime
	end
	def set_trait(trait)
		@trait = trait
	end
	def set_trait_stderr(trait_stderr)
		@trait_stderr = trait_stderr
	end
	def to_s
		str = "#{@number}\t#{@species_id}\t#{@ancestor}\t#{format("%.4f", @time)}\t#{@regime}\t"
		if @trait == "NA"
			str << "NA\t"
		else
			str << "#{format("%.4f", @trait)}\t"
		end
		if @trait_stderr == "NA"
			str << "NA\n"
		else
			str << "#{format("%.4f", @trait_stderr)}\n"
		end
		str
	end
end

# Get the command line arguments.
tree_file_name = ARGV[0]
trait_file_name = ARGV[1]
slouch_input_directory = ARGV[2]

# Make the slouch directory if it doesn't exist yet.
FileUtils.mkdir_p(slouch_input_directory)

# Read the trait data file.
trait_file = File.open(trait_file_name)
trait_lines = trait_file.readlines
species_ids = []
traits = []
trait_stderrs = []
trait_lines.each do |l|
	line_ary = l.split("\t")
	species_ids << line_ary[0].strip
	traits << line_ary[1].strip.to_f
	trait_stderrs << line_ary[2].strip.to_f
end

# Parse the tree.
tree = Tree.parse(tree_file_name, fileType = "nexus", diversityFileName = nil, treeNumber = 0, verbose = false)

# Find the tree root age.
root_age = 0
tree.branch.each do |b|
	root_age = b.origin if b.origin > root_age
end

# Sort branches by their termination age.
branches = tree.branch
sorted = false
until sorted
	sorted = true
	0.upto(branches.size-2) do |x|
		if branches[x].termination < branches[x+1].termination
			branches[x], branches[x+1] = branches[x+1], branches[x]
			sorted = false
		end
	end
end

# Prepare all table lines for the slouch input table.
table_lines = []
table_lines << TableLine.new(1, "root", "NA", 0)
table_lines.last.set_ancestor("0")
table_lines.last.set_regime(1)
table_lines.last.set_trait("NA")
table_lines.last.set_trait_stderr("NA")
branch_count = 1
branches.each do |b|
	branch_count += 1
	if b.extant
		table_lines << TableLine.new(branch_count, b.id, b.speciesId, root_age - b.termination)
	else
		table_lines << TableLine.new(branch_count, b.id, "NA", root_age - b.termination)
	end
end

# Add the ancestor of each branch.
table_lines[1..-1].each do |t|
	parent_id = nil
	branches.each do |b|
		if b.id == t.branch_id
			parent_id = b.parentId
			break
		end
	end
	raise "ERROR: Parent id could not be found!" if parent_id == nil
	if parent_id == "treeOrigin"
		ancestor = "1"
	else
		ancestor = nil
		table_lines.each do |tt|
			if tt.branch_id == parent_id
				ancestor = tt.number
				break
			end
		end
		raise "ERROR: Ancestor could not be found!" if ancestor == nil
	end
	t.set_ancestor(ancestor)
end

# The regime code specifies which of the seven non-background regimes should be used.
# regime_code = "11111111" means that all regimes are used, a total of 9.
# regime_code = "00000000" means that the entire phylogeny has the background regime.
0.upto(255) do |rc|
	regime_code = rc.to_s(2).rjust(8).gsub(" ","0")
	table_lines_copy = Marshal.load(Marshal.dump(table_lines))

	# Add the regime of each branch.
	table_lines_copy[1..-1].each do |t|
		# The background regime is '1'.
		t.set_regime(1) if t.regime == nil
		branches.each do |b|
			if b.id == t.branch_id
				extant_species_of_this_branch = []
				if b.extant
					extant_species_of_this_branch << b.speciesId
				else
					b.extantProgenyId.each do |p|
						branches.each do |bb|
							extant_species_of_this_branch << bb.speciesId if bb.id == p and bb.extant
						end
					end
				end
				# The regime of Percomorphaceae excluding Ophidiiformes is '2'.
				if regime_code[0] == "1"
					if extant_species_of_this_branch.include?("fish_40") and extant_species_of_this_branch.include?("fish_91")
						unless extant_species_of_this_branch.include?("fish_31")
							regime = 2
							t.set_regime(regime)
							table_lines_copy[1..-1].each do |tt|
								tt.set_regime(regime) if b.progenyId.include?(tt.branch_id)
							end
						end
					end
				end
				# The regime of Berycoidei is '3'.
				if regime_code[1] == "1"
					if extant_species_of_this_branch.include?("fish_69")
						unless extant_species_of_this_branch.include?("fish_79")
							regime = 3
							t.set_regime(regime)
							table_lines_copy[1..-1].each do |tt|
								tt.set_regime(regime) if b.progenyId.include?(tt.branch_id)
							end
						end
					end
				end
				# The regime of Myctophiformes is '4'.
				if regime_code[2] == "1"
					if extant_species_of_this_branch.include?("fish_66")
						unless extant_species_of_this_branch.include?("fish_24")
							regime = 4
							t.set_regime(regime)
							table_lines_copy[1..-1].each do |tt|
								tt.set_regime(regime) if b.progenyId.include?(tt.branch_id)
							end
						end
					end
				end
				# The regime of Paracanthomorphacea is '5'.
				if regime_code[3] == "1"
					if extant_species_of_this_branch.include?("fish_24") and extant_species_of_this_branch.include?("fish_97")
						unless extant_species_of_this_branch.include?("fish_47")
							regime = 5
							t.set_regime(regime)
							table_lines_copy[1..-1].each do |tt|
								tt.set_regime(regime) if b.progenyId.include?(tt.branch_id)
							end
						end
					end
				end
				# The regime of Polymixiiformes+Percopsiformes is '6'.
				if regime_code[4] == "1"
					if extant_species_of_this_branch.include?("fish_24") and extant_species_of_this_branch.include?("fish_26")
						unless extant_species_of_this_branch.include?("fish_28")
							regime = 6
							t.set_regime(regime)
							table_lines_copy[1..-1].each do |tt|
								tt.set_regime(regime) if b.progenyId.include?(tt.branch_id)
							end
						end
					end
				end
				# The regime of Gadiformes including Bregmacerotidae is '7'.
				if regime_code[5] == "1"
					if extant_species_of_this_branch.include?("fish_21") and extant_species_of_this_branch.include?("fish_15")
						unless extant_species_of_this_branch.include?("fish_80")
							regime = 7
							t.set_regime(regime)
							table_lines_copy[1..-1].each do |tt|
								tt.set_regime(regime) if b.progenyId.include?(tt.branch_id)
							end
						end
					end
				end
				# The regime of Gadiformes without Bregmacerotidae is '8'.
				if regime_code[6] == "1"
					if extant_species_of_this_branch.include?("fish_97") and extant_species_of_this_branch.include?("fish_15")
						unless extant_species_of_this_branch.include?("fish_21")
							regime = 8
							t.set_regime(regime)
							table_lines_copy[1..-1].each do |tt|
								tt.set_regime(regime) if b.progenyId.include?(tt.branch_id)
							end
						end
					end
				end
				# The regime of Southern Gadiformes is '9'.
				if regime_code[7] == "1"
					if extant_species_of_this_branch.include?("fish_18") and extant_species_of_this_branch.include?("fish_15")
						unless extant_species_of_this_branch.include?("fish_97")
							regime = 9
							t.set_regime(regime)
							table_lines_copy[1..-1].each do |tt|
								tt.set_regime(regime) if b.progenyId.include?(tt.branch_id)
							end
						end
					end
				end
				break
			end
		end
	end

	# Add the trait and trait standard errors to each branch.
	table_lines_copy[1..-1].each do |t|
		if t.species_id == "NA"
			t.set_trait("NA")
			t.set_trait_stderr("NA")
		else
			t.set_trait(traits[species_ids.index(t.species_id)])
			t.set_trait_stderr(trait_stderrs[species_ids.index(t.species_id)])
		end
	end		

	# Prepare a string for the slouch table file.
	slouch_table_string = ""
	slouch_table_string << "	species	ancestor	time	regime	trait	trait.me\n"
	table_lines_copy.each {|t| slouch_table_string << t.to_s}

	# Write the slouch table to file.
	slouch_table_file_name = "#{slouch_input_directory}/#{regime_code}.txt"
	slouch_table_file = File.open(slouch_table_file_name,"w")
	slouch_table_file.write(slouch_table_string)
	slouch_table_file.close

end # 0.upto(255) do |rc|
