$libPath = "./resources/phylsim/"
require "./resources/phylsim/main.rb"

class Line
	attr_reader :x_start, :x_end, :y_start, :y_end, :z_start, :z_end, :color
	def initialize(x_start,x_end,y_start,y_end,z_start,z_end,color,stroke)
		@x_start = x_start
		@x_end = x_end
		@y_start = y_start
		@y_end = y_end
		@z_start = z_start
		@z_end = z_end
		@original_x_start = x_start
		@original_x_end = x_end
		@original_y_start = y_start
		@original_y_end = y_end
		@original_z_start = z_start
		@original_z_end = z_end
		@color = color
		@stroke = stroke
	end
	def to_s
		string = ""
		string << "Start:  #{@x_start.round(3)},#{@y_start.round(3)},#{@z_start.round(3)}\n"
		string << "End:    #{@x_end.round(3)},#{@y_end.round(3)},#{@z_end.round(3)}\n"
		string << "Color:  #{@color}\n"
		string << "Stroke: #{@stroke}\n"
		string << "\n"
		string
	end
	def to_svg
		svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)-@z_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{(@y_end-@z_end).round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\"  />"
		svg
	end
	def rotate(angle_in_degrees,x_center,y_center)
		# See http://stackoverflow.com/questions/13695317/rotate-a-point-around-another-point
		angle_in_radians = angle_in_degrees * (Math::PI/180.0)
		cos_theta = Math.cos(angle_in_radians)
		sin_theta = Math.sin(angle_in_radians)
		new_x_start = cos_theta * (@x_start - x_center) - sin_theta * (@y_start - y_center) + x_center
		new_y_start = sin_theta * (@x_start - x_center) + cos_theta * (@y_start - y_center) + y_center
		new_x_end = cos_theta * (@x_end - x_center) - sin_theta * (@y_end - y_center) + x_center
		new_y_end = sin_theta * (@x_end - x_center) + cos_theta * (@y_end - y_center) + y_center
		@x_start = new_x_start
		@x_end = new_x_end
		@y_start = new_y_start
		@y_end = new_y_end
	end
	def tilt(angle_in_degrees,y_center)
		angle_in_radians = angle_in_degrees * (Math::PI/180.0)
		tan_theta = Math.sin(angle_in_radians)
		new_y_start = y_center - tan_theta * (y_center-y_start)
		new_y_end = y_center - tan_theta * (y_center-y_end)
		@y_start = new_y_start
		@y_end = new_y_end
	end
	def max_y
		[@y_start,@y_end].max
	end
	def min_original_x
		[@original_x_start,@original_x_end].min
	end
	def max_original_x
		[@original_x_start,@original_x_end].max
	end
	def avg_original_x
		(@original_x_start+@original_x_end)/2.0
	end
	def min_original_y
		[@original_y_start,@original_y_end].min
	end
	def max_original_y
		[@original_y_start,@original_y_end].max
	end
end

class Polygon
	attr_reader :x1, :y1, :z1, :x2, :y2, :z2, :x3, :y3, :z3, :x4 ,:y4, :z4, :fill, :stroke, :stroke_width, :opacity
	def initialize(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, fill, stroke, stroke_width, opacity)
		@x1 = x1
		@y1 = y1
		@z1 = z1
		@x2 = x2
		@y2 = y2
		@z2 = z2
		@x3 = x3
		@y3 = y3
		@z3 = z3
		@x4 = x4
		@y4 = y4
		@z4 = z4
		@original_x1 = x1
		@original_y1 = y1
		@original_z1 = z1
		@original_x2 = x2
		@original_y2 = y2
		@original_z2 = z2
		@original_x3 = x3
		@original_y3 = y3
		@original_z3 = z3
		@original_x4 = x4
		@original_y4 = y4
		@original_z4 = z4
		@fill = fill
		@stroke = stroke
		@stroke_width = stroke_width
		@opacity = opacity
	end
	def to_s
		string = ""
		string << "Point 1:      #{@x1.round(3)},#{@y1.round(3)},#{@z1.round(3)}\n"
		string << "Point 2:      #{@x2.round(3)},#{@y2.round(3)},#{@z2.round(3)}\n"
		string << "Point 3:      #{@x3.round(3)},#{@y3.round(3)},#{@z3.round(3)}\n"
		string << "Point 4:      #{@x4.round(3)},#{@y4.round(3)},#{@z4.round(3)}\n"
		string << "Fill:         #{@fill}\n"
		string << "Stroke:       #{@stroke}\n"
		string << "Stroke-width: #{@stroke_width}\n"
		string << "Opacity:      #{@opacity}\n"
		string << "\n"
		string
	end
	def to_svg
		svg = "<polygon points=\"#{@x1.round(3)},#{(@y1-@z1).round(3)} #{@x2.round(3)},#{(@y2-@z2).round(3)} #{@x3.round(3)},#{(@y3-@z3).round(3)} #{@x4.round(3)},#{(@y4-@z4).round(3)}\" fill=\"#{@fill}\" stroke=\"#{@stroke}\" stroke-width=\"#{@stroke_width}\" opacity=\"#{@opacity}\" />"
		svg
	end
	def rotate(angle_in_degrees,x_center,y_center)
		# See http://stackoverflow.com/questions/13695317/rotate-a-point-around-another-point
		angle_in_radians = angle_in_degrees * (Math::PI/180.0)
		cos_theta = Math.cos(angle_in_radians)
		sin_theta = Math.sin(angle_in_radians)
		new_x1 = cos_theta * (@x1 - x_center) - sin_theta * (@y1 - y_center) + x_center
		new_y1 = sin_theta * (@x1 - x_center) + cos_theta * (@y1 - y_center) + y_center
		new_x2 = cos_theta * (@x2 - x_center) - sin_theta * (@y2 - y_center) + x_center
		new_y2 = sin_theta * (@x2 - x_center) + cos_theta * (@y2 - y_center) + y_center
		new_x3 = cos_theta * (@x3 - x_center) - sin_theta * (@y3 - y_center) + x_center
		new_y3 = sin_theta * (@x3 - x_center) + cos_theta * (@y3 - y_center) + y_center
		new_x4 = cos_theta * (@x4 - x_center) - sin_theta * (@y4 - y_center) + x_center
		new_y4 = sin_theta * (@x4 - x_center) + cos_theta * (@y4 - y_center) + y_center
		@x1 = new_x1
		@x2 = new_x2
		@x3 = new_x3
		@x4 = new_x4
		@y1 = new_y1
		@y2 = new_y2
		@y3 = new_y3
		@y4 = new_y4
	end
	def tilt(angle_in_degrees,y_center)
		angle_in_radians = angle_in_degrees * (Math::PI/180.0)
		tan_theta = Math.sin(angle_in_radians)
		new_y1 = y_center - tan_theta * (y_center-y1)
		new_y2 = y_center - tan_theta * (y_center-y2)
		new_y3 = y_center - tan_theta * (y_center-y3)
		new_y4 = y_center - tan_theta * (y_center-y4)
		@y1 = new_y1
		@y2 = new_y2
		@y3 = new_y3
		@y4 = new_y4
	end
	def max_y
		[@y1,@y2,@y3,@y4].max
	end
	def min_original_x
		[@original_x1,@original_x2,@original_x3,@original_x4].min
	end
	def max_original_x
		[@original_x1,@original_x2,@original_x3,@original_x4].max
	end
	def avg_original_x
		(@original_x1+@original_x2+@original_x3+@original_x4)/4.0
	end
	def min_original_y
		[@original_y1,@original_y2,@original_y3,@original_y4].min
	end
	def max_original_y
		[@original_y1,@original_y2,@original_y3,@original_y4].max
	end
end

# Read the tree.
tree_file_name = ARGV[0]
tree = Tree.parse(tree_file_name,"newick",nil,0,false)

# Read the table with regime specifications and trait values.
table_file_name = ARGV[1]
table_file = File.open(table_file_name)
table_lines = table_file.readlines
table_species = []
table_regimes = []
table_traits = []
table_lines[1..-1].each do |l|
	line_ary = l.split
	unless line_ary[1] == "NA"
		table_species << line_ary[1]
		table_regimes << line_ary[4].to_i
		table_traits << line_ary[5].to_f
	end
end

# Read the file that specifies the order in which tips should appear in the figure.
order_file_name = ARGV[2]
order_file = File.open(order_file_name)
order_lines = order_file.readlines
order = []
order_lines.each do |l|
	order << l.strip
end

# Read the Slouch results file.
res_file_name = ARGV[3]
res_file = File.open(res_file_name)
res_lines = res_file.readlines
alpha = nil
sigma_square = nil
regimes_in_res = []
optima_in_res = []
root_optimum = []
in_optima_part = false
in_estimates_part = false
res_lines.each do |l|
	# Read sigma square and alpha.
	line_ary = l.split(/\s\s+/)
	sigma_square = line_ary[1].to_f if line_ary[0] == "Stationary variance"
	alpha = line_ary[1].to_f if line_ary[0] == "Rate of adaptation"
	# Read the optima.
	in_optima_part = true if l.strip == "PRIMARY OPTIMA"
	in_estimates_part = true if l.strip == "Estimates Std.error" and in_optima_part
	in_optima_part = false if l.include?("-------------") or l.include?("=============")
	in_estimates_part = false if l.include?("-------------") or l.include?("=============")
	if in_estimates_part
		unless l.strip == "Estimates Std.error" or l.strip == ""
			line_ary = l.split
			raise "ERROR: Slouch result file #{res_file_name} could not be parsed!" if line_ary.size == 2
			regimes_in_res << line_ary[0].to_i
			optima_in_res << line_ary[1].to_f
		end
	end
end
root_optimum = optima_in_res[regimes_in_res.index(regimes_in_res.min)]

# Assign y axis positions to all tips.
branch_pos_on_y_axis = []
tree.branch.size.times {|x| branch_pos_on_y_axis[x] = "unknown"}

tree.branch.size.times do |x|
	if order.include?(tree.branch[x].speciesId)
		branch_pos_on_y_axis[x] = order.index(tree.branch[x].speciesId)
	end
end

# Assign y axis positions to all internal branches.
all_y_axis_pos_assigned = false
while all_y_axis_pos_assigned == false
	# Initiate a variable to check whether the loop has finished.
	still_unknown_positions = false
	# For all branches that have no y axis position yet.
	tree.branch.size.times do |x|
		if branch_pos_on_y_axis[x] == "unknown"
			# If unknown positions are still present, set the check variable so that the loop is repeated once more.
			still_unknown_positions = true
			# Check whether the two daughters already have y axis positions.
			daughter1_y_pos = nil
			daughter2_y_pos = nil
			tree.branch.size.times do |y|
				if tree.branch[x].daughterId[0] == tree.branch[y].id
					daughter1_y_pos = branch_pos_on_y_axis[y]
					break
				end
			end
			tree.branch.size.times do |y|
				if tree.branch[x].daughterId[1] == tree.branch[y].id
					daughter2_y_pos = branch_pos_on_y_axis[y]
					break
				end
			end
			raise "ERROR: Daughter 1 could not be found!" if daughter1_y_pos == nil
			raise "ERROR: Daughter 2 could not be found!" if daughter2_y_pos == nil
			# If both daughters already have y axis positions, this branch should receive a position exactly in the middle between them.
			if daughter1_y_pos != "unknown" and daughter2_y_pos != "unknown"
				branch_pos_on_y_axis[x] = (daughter1_y_pos+daughter2_y_pos)/2.0
			end
		end
	end
	# Finish the loop if no changes were made in the last round.
	all_y_axis_pos_assigned = true unless still_unknown_positions
end

# Initiate an array to hold regime assignments and traits for each branch (for now only extant branches).
branch_regimes = []
branch_traits_end = []
tree.branch.size.times do |x|
	if tree.branch[x].endCause == "present"
		branch_regimes << table_regimes[table_species.index(tree.branch[x].speciesId)]
		branch_traits_end << table_traits[table_species.index(tree.branch[x].speciesId)]
	else
		branch_regimes << "unknown"
		branch_traits_end << "unknown"
	end
end

# Find the highest regime number.
max_regime = 0
branch_regimes.each do |r|
	max_regime = r if r > max_regime unless r == "unknown"
end

# For each internal branch, assign the smallest regime number found in the progeny.
tree.branch.size.times do |x|
	unless tree.branch[x].endCause == "present"
		lowest_regime_in_progeny = max_regime
		tree.branch.size.times do |y|
			if tree.branch[x].extantProgenyId.include?(tree.branch[y].id)
				lowest_regime_in_progeny = branch_regimes[y] if branch_regimes[y] < lowest_regime_in_progeny
			end
		end
		branch_regimes[x] = lowest_regime_in_progeny
	end
end

# Make sure, no branch regimes are unknown anymore.
tree.branch.size.times do |x|
	raise "Regime of branch #{tree.branch[x].id} is unknown!" if branch_regimes[x] == "unknown"
end

# For each branch, assign the optimum based on regime.
branch_optima = []
tree.branch.size.times do |x|
	branch_optima[x] = optima_in_res[regimes_in_res.index(branch_regimes[x])]
end

# For each branch, calculate the predicted mean at its start and end.
branch_predicted_mean_start = []
branch_predicted_mean_end = []
# First set the predicted mean of the start of the two root branches.
branch_predicted_mean_start[0] = root_optimum
branch_predicted_mean_start[1] = root_optimum
# Prepare a loop that goes through all branches until no more branches have predicted means at their start, but not at their end.
continue = true
while continue
	continue = false
	tree.branch.size.times do |x|
		if branch_predicted_mean_start[x] != nil and branch_predicted_mean_end[x] == nil
			continue = true
			# Calculate the predicted mean at branch end, given its start value, the regime's optimum, and alpha.
			branch_predicted_mean_end[x] = branch_predicted_mean_start[x] + (branch_optima[x]-branch_predicted_mean_start[x]) - Math.exp(-alpha*tree.branch[x].duration)*(branch_optima[x]-branch_predicted_mean_start[x])
			# Find the two daughters and transfer the predicted mean end value to them as their start value
			tree.branch.size.times do |y|
				branch_predicted_mean_start[y] = branch_predicted_mean_end[x] if tree.branch[x].daughterId.include?(tree.branch[y].id)
			end
		end
	end
end


################################
# The matrix calculation part
################################

# Prepare an array for the extant residuals (that's called "y-m" in Thomas' equation).
extant_trait_values = []
extant_predicted_means = []
tree.branch.size.times do |x|
	if tree.branch[x].extant
		extant_trait_values << branch_traits_end[x]
		extant_predicted_means << branch_predicted_mean_end[x]
	end
end
extant_residuals_row = []
extant_trait_values.size.times {|x| extant_residuals_row[x] = extant_trait_values[x]-extant_predicted_means[x]}
extant_residuals_rows = []
extant_residuals_rows << extant_residuals_row
extant_residuals = Matrix.rows(extant_residuals_rows).transpose

# Prepare the matrix of extant covariance (called "Var[Y]" in Thomas' equation).
extant_covariance_rows = []
tree.branch.size.times do |x|
	if tree.branch[x].extant
		extant_covariance_row = []
		tree.branch.size.times do |y|
			if tree.branch[y].extant
				# Find the sum of branch lengths txy that separates the two extant species tree.branch[x] and tree.branch[y].
				txy = nil
				if x == y
					txy = 0
				else
					# Unless x = y, txy is twice the termination age of the youngest branch that still has both species in its extant progeny.
					termination_ages_of_branches_ancestral_to_both = [tree.branch[0].origin] # So that the root age is always included.
					tree.branch.size.times do |z|
						unless tree.branch[z].extant
							if tree.branch[z].extantProgenyId.include?(tree.branch[x].id) and tree.branch[z].extantProgenyId.include?(tree.branch[y].id)
								termination_ages_of_branches_ancestral_to_both << tree.branch[z].termination
							end
						end
					end
					# As tree.branch[x].termination and tree.branch[y].termination should be both 0 (both species are extant), this is equivalent to twice the minimum termination age.
					tx = termination_ages_of_branches_ancestral_to_both.min - tree.branch[x].termination
					ty = termination_ages_of_branches_ancestral_to_both.min - tree.branch[y].termination
					txy = tx + ty
				end
				extant_covariance_cell = sigma_square/(2*alpha) * Math.exp(-alpha*txy)
				extant_covariance_row << extant_covariance_cell
			end
		end
		extant_covariance_rows << extant_covariance_row
	end
end
extant_covariance = Matrix.rows(extant_covariance_rows)
extant_covariance_inverse = extant_covariance.inverse

# Prepare the matrix of ancestral covariance (called "Var[A,Y]" in Thomas' equation).
ancestral_covariance_rows = []
tree.branch.size.times do |x|
	unless tree.branch[x].extant
		ancestral_covariance_row = []
		tree.branch.size.times do |y|
			if tree.branch[y].extant
				# Find the sum of branch lengths txy that separates the internal node tree.branch[x] (its termination) and the extant tree.branch[y].
				txy = nil
				if tree.branch[x].extantProgenyId.include?(tree.branch[y].id)
					txy = tree.branch[x].termination
				else
					# If tree.branch[y] is no descendant of tree.branch[x], a species is picked from the extant progeny of tree.branch[x],
					# and the youngest branch that still has both this picked species and tree.branch[y] in its extant progeny is selected.
					termination_ages_of_branches_ancestral_to_both = [tree.branch[0].origin] # So that the root age is always included.
					tree.branch.size.times do |z|
						unless tree.branch[z].extant
							if tree.branch[z].extantProgenyId.include?(tree.branch[x].extantProgenyId[0]) and tree.branch[z].extantProgenyId.include?(tree.branch[y].id)
								termination_ages_of_branches_ancestral_to_both << tree.branch[z].termination
							end
						end
					end
					# This time, only tree.branch[y].termination should be 0.
					tx = termination_ages_of_branches_ancestral_to_both.min - tree.branch[x].termination
					ty = termination_ages_of_branches_ancestral_to_both.min - tree.branch[y].termination
					txy = tx + ty
				end
				ancestral_covariance_cell = sigma_square/(2*alpha) * Math.exp(-alpha*txy)
				ancestral_covariance_row << ancestral_covariance_cell
			end
		end
		ancestral_covariance_rows << ancestral_covariance_row
	end
end
ancestral_covariance = Matrix.rows(ancestral_covariance_rows)

# Calculate the ancestral residuals matrix according to Thomas' equation, and transform it into an array.
ancestral_residuals = (ancestral_covariance * extant_covariance_inverse * extant_residuals).transpose.to_a[0]

# Prepare an array for the ancestral predicted means.
ancestral_predicted_means = []
tree.branch.size.times do |x|
	ancestral_predicted_means << branch_predicted_mean_end[x] unless tree.branch[x].extant
end

# Prepare the array for the ancestral estimates.
ancestral_estimates = []
ancestral_residuals.size.times do |x|
	ancestral_estimates[x] = ancestral_residuals[x] + ancestral_predicted_means[x]
end

# Assign ancestral estimates as end traits to non-extant branches.
# For extant branches, the actual trait values have been added already above.
tree.branch.size.times do |x|
	branch_traits_end[x] = ancestral_estimates.shift unless tree.branch[x].extant
end

# For all branches except the root branches, use the end trait of the parent as the start trait.
branch_traits_start = []
tree.branch.size.times do |x|
	unless tree.branch[x].parentId == "treeOrigin"
		tree.branch.size.times do |y|
			if tree.branch[y].id == tree.branch[x].parentId
				branch_traits_start[x] = branch_traits_end[y]
			end
		end
	end
end

# Add the optimum of the most ancestral regime as the start trait of the two root branches.
branch_traits_start[0] = root_optimum
branch_traits_start[1] = root_optimum

# Transform all traits back to the original scale instead of the lognormal scale.
logscale = false
if logscale == false
	branch_traits_start.size.times do |x|
		branch_traits_start[x] = Math.exp(branch_traits_start[x])
		branch_traits_end[x] = Math.exp(branch_traits_end[x])
	end
end

################################
# The graphics part
################################

# Specifation of some graphic parameters before calculating lines.
x_dim = 600
y_dim = 600
z_dim = 300 # Adjust!
z_add = 0 # Adjust!
margin = 100
colors = ["#859900","#cb4b16","#2aa198","#268bd2","#6c71c4","#d33682", "#dc322f", "#b58900"]
polygon_fill_color = "#d6d6d6"
polygon_stroke = "white"
polygon_stroke_width = 0.5
polygon_opacity = 0.5

# Initiate an array for all lines.
elements = []

# Find the maximum age, the maximum trait value, and the maximum y position.
max_age = 0
max_trait = 0
max_y = 0
tree.branch.size.times do |x|
	max_age = tree.branch[x].origin if tree.branch[x].origin > max_age
	max_y = branch_pos_on_y_axis[x] if branch_pos_on_y_axis[x] > max_y
	max_trait = branch_traits_start[x] if branch_traits_start[x] > max_trait
	max_trait = branch_traits_end[x] if branch_traits_end[x] > max_trait
end

# Add a line to connect the two root branches.
x_start = margin+((max_age-tree.branch[0].origin)/max_age)*(x_dim-2*margin)
x_end = margin+((max_age-tree.branch[1].origin)/max_age)*(x_dim-2*margin)
y_start = margin+(branch_pos_on_y_axis[0]/max_y.to_f)*(y_dim-2*margin)
y_end = margin+(branch_pos_on_y_axis[1]/max_y.to_f)*(y_dim-2*margin)
z_start = (branch_traits_start[0]/max_trait)*z_dim + z_add
z_end = (branch_traits_start[1]/max_trait)*z_dim + z_add
# lines << Line.new(x_start,x_end,y_start,y_end,0,0,colors[0],1)
elements << Line.new(x_start,x_end,y_start,y_end,z_start,z_end,colors[0],1)
x1 = x_start
y1 = y_start
z1 = 0
x2 = x_end
y2 = y_end
z2 = 0
x3 = x_end
y3 = y_end
z3 = z_end
x4 = x_start
y4 = y_start
z4 = z_start
elements << Polygon.new(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, polygon_fill_color, polygon_stroke, polygon_stroke_width, polygon_opacity)

# Walk down the tree to create SVG elements in a sorted manner.
branches_drawn = []
all_descendants_drawn = []
tree.branch.size.times do |x|
	branches_drawn << false
	all_descendants_drawn << false
end
this_branch_index = 0
while branches_drawn.include?(false)

	# If this branch has not been drawn already, do so.
	unless branches_drawn[this_branch_index]
		# Add a line for this branch.
		x_start = margin+((max_age-tree.branch[this_branch_index].origin)/max_age)*(x_dim-2*margin)
		x_end = margin+((max_age-tree.branch[this_branch_index].termination)/max_age)*(x_dim-2*margin)
		y_start = margin+(branch_pos_on_y_axis[this_branch_index]/max_y.to_f)*(y_dim-2*margin)
		y_end = margin+(branch_pos_on_y_axis[this_branch_index]/max_y.to_f)*(y_dim-2*margin)
		z_start = (branch_traits_start[this_branch_index]/max_trait)*z_dim + z_add
		z_end = (branch_traits_end[this_branch_index]/max_trait)*z_dim + z_add
		elements << Line.new(x_start,x_end,y_start,y_end,z_start,z_end,colors[branch_regimes.uniq.sort.index(branch_regimes[this_branch_index])],1)
		x1 = x_start
		y1 = y_start
		z1 = 0
		x2 = x_end
		y2 = y_end
		z2 = 0
		x3 = x_end
		y3 = y_end
		z3 = z_end
		x4 = x_start
		y4 = y_start
		z4 = z_start
		elements << Polygon.new(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, polygon_fill_color, polygon_stroke, polygon_stroke_width, polygon_opacity)
		# Add a line to connect daughters with each other.
		daughter1_y_pos = nil
		daughter2_y_pos = nil
		unless tree.branch[this_branch_index].daughterId == ["none","none"]
			connecting_branch_regime = nil
			connecting_branch_regime1 = nil
			connecting_branch_regime2 = nil
			tree.branch.size.times do |y|
				if tree.branch[this_branch_index].daughterId[0] == tree.branch[y].id
					daughter1_y_pos = branch_pos_on_y_axis[y]
					connecting_branch_regime1 = branch_regimes[y]
					break
				end
			end
			tree.branch.size.times do |y|
				if tree.branch[this_branch_index].daughterId[1] == tree.branch[y].id
					daughter2_y_pos = branch_pos_on_y_axis[y]
					connecting_branch_regime2 = branch_regimes[y]
					break
				end
			end
			connecting_branch_regime = [connecting_branch_regime1,connecting_branch_regime2].min
			x_start = margin+((max_age-tree.branch[this_branch_index].termination)/max_age)*(x_dim-2*margin)
			x_end = margin+((max_age-tree.branch[this_branch_index].termination)/max_age)*(x_dim-2*margin)
			y_start = margin+(daughter1_y_pos/max_y.to_f)*(y_dim-2*margin)
			y_end = margin+(daughter2_y_pos/max_y.to_f)*(y_dim-2*margin)
			z_start = (branch_traits_end[this_branch_index]/max_trait)*z_dim + z_add
			z_end = (branch_traits_end[this_branch_index]/max_trait)*z_dim + z_add
			elements << Line.new(x_start,x_end,y_start,y_end,z_start,z_end,colors[branch_regimes.uniq.sort.index(connecting_branch_regime)],1)
			x1 = x_start
			y1 = y_start
			z1 = 0
			x2 = x_end
			y2 = y_end
			z2 = 0
			x3 = x_end
			y3 = y_end
			z3 = z_end
			x4 = x_start
			y4 = y_start
			z4 = z_start
			elements << Polygon.new(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, polygon_fill_color, polygon_stroke, polygon_stroke_width, polygon_opacity)
		end
		branches_drawn[this_branch_index] = true
	end

	# Find out if all descendants have been drawn.
	all_descendants_drawn_of_this_branch = false
	if tree.branch[this_branch_index].daughterId == ["none","none"]
		all_descendants_drawn_of_this_branch = true
	else
		# Find the two daughters and check if the descendants of both of them have been drawn.
		this_branch_daughter0_index = nil
		this_branch_daughter1_index = nil
		tree.branch.size.times do |xx|
			this_branch_daughter0_index = xx if tree.branch[xx].id == tree.branch[this_branch_index].daughterId[0]
			this_branch_daughter1_index = xx if tree.branch[xx].id == tree.branch[this_branch_index].daughterId[1]
		end
		if branch_pos_on_y_axis[this_branch_daughter0_index] > branch_pos_on_y_axis[this_branch_daughter1_index]
			this_branch_daughter0_index,this_branch_daughter1_index = this_branch_daughter1_index,this_branch_daughter0_index
		end
		all_descendants_drawn_of_this_branch = true if all_descendants_drawn[this_branch_daughter0_index] and all_descendants_drawn[this_branch_daughter1_index]
	end

	# If all branches descending from this branch have been drawn already, mark this branch accordingly, and then change the current branch again to one of the root branches.
	if all_descendants_drawn_of_this_branch
		all_descendants_drawn[this_branch_index] = true
		if all_descendants_drawn[0] == false and all_descendants_drawn[1] == false
			this_branch_index = 0
		elsif all_descendants_drawn[0] == true and all_descendants_drawn[1] == false
			this_branch_index = 1
		else
			raise "ERROR!"
		end
	else
		# If not, change the current branch.
		if all_descendants_drawn[this_branch_daughter0_index] == false
			this_branch_index = this_branch_daughter0_index
		elsif all_descendants_drawn[this_branch_daughter1_index] == false
			this_branch_index = this_branch_daughter1_index
		end
	end

end

# Rotate all lines around the center.
elements.each {|e| e.rotate(-105,x_dim/2.0,y_dim/2.0)}

# Tilt all lines with respect to the lower margin.
elements.each {|e| e.tilt(15,y_dim-margin)}

# Change the order of all elements, so that the first ones are printed last, and thus appear on top.
elements.reverse!

# Prepare the svg string.
svg_string = ""
svg_string << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
svg_string << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
svg_string << "<svg width=\"#{x_dim}\" height=\"#{y_dim}\" viewBox=\"0 0 #{x_dim} #{y_dim}\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n"
elements.each {|e| svg_string << "  #{e.to_svg}\n"}
svg_string << "</svg>\n"

# Write the svg string to file.
svg_file_name = ARGV[4]
svg_file = File.new(svg_file_name,"w")
svg_file.write(svg_string)

# Prepare a report of ancestral states for each branch.
report_string = "#{"branch_id".rjust(14)}\t#{"origin".rjust(14)}\t#{"termination".rjust(14)}\t#{"species_id".rjust(14)}\t#{"regime".rjust(14)}\t#{"start_trait".rjust(14)}\t#{"end_trait".rjust(14)}\n"
tree.branch.size.times do |x|
	report_string << "#{tree.branch[x].id.rjust(14)}\t#{(format("%.2f", tree.branch[x].origin).rjust(14))}\t#{(format("%.2f", tree.branch[x].termination).rjust(14))}\t#{tree.branch[x].speciesId.rjust(14)}\t#{branch_regimes[x].to_s.rjust(14)}\t#{(format("%.2f", branch_traits_start[x]).rjust(14))}\t#{(format("%.2f", branch_traits_end[x]).rjust(14))}\n"
end

# Write the report to file.
report_file_name = ARGV[5]
report_file = File.new(report_file_name,"w")
report_file.write(report_string)
