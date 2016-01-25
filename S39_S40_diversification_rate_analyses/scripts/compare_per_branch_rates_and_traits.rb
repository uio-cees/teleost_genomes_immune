# Michael Matschiner, 2015-07-29

# Get command line arguments.
branches_net_div_file_name = ARGV[0]
trait_reconstruction_file_name = ARGV[1]
output_table_file_name = ARGV[2]
output_svg_file_name = ARGV[3]

# Read the file with net diversification estimates per branch.
branches_net_div_file = File.open(branches_net_div_file_name)
branches_net_div_lines = branches_net_div_file.readlines
branch_begins = []
branch_ends = []
branch_net_divs = []
branch_tip_labels = []
tip_labels = []
in_branch_begins = false
in_branch_ends = false
in_branch_net_divs = false
in_branch_tip_labels = false
branches_net_div_lines.each do |l|
	if l.strip == "Branch begin"
		in_branch_begins = true
	elsif l.strip == "Branch end"
		in_branch_ends = true
		in_branch_begins = false
	elsif l.strip == "Branch rate"
		in_branch_net_divs = true
		in_branch_ends = false
	elsif l.strip == "Branch tip label"
		in_branch_tip_labels = true
		in_branch_net_divs = false
	elsif in_branch_begins
		branch_begins << l.strip.to_f unless l.strip == ""
	elsif in_branch_ends
		branch_ends << l.strip.to_f unless l.strip == ""
	elsif in_branch_net_divs
		branch_net_divs << l.strip.to_f unless l.strip == ""
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

# Read the trait reconstruction file.
trait_reconstruction_file = File.open(trait_reconstruction_file_name)
trait_reconstruction_lines = trait_reconstruction_file.readlines

# For each branch, identify regime, mean trait, branch length, and net diversification.
regime_per_branch = []
mean_trait_per_branch = []
branch_length_per_branch = []
extant_per_branch = []
net_diversification_per_branch = []
trait_reconstruction_lines[1..-1].each do |l|
	line_ary = l.split
	origin = line_ary[1].to_f
	termination = line_ary[2].to_f
	if termination == 0.0
		extant_per_branch << true
	else
		extant_per_branch << false
	end
	species = line_ary[3]
	net_diversification_candidates = []
	branch_begins.size.times do |x|
		if (branch_max_age-branch_begins[x]).round(2) == origin.round(2) and (branch_max_age-branch_ends[x]).round(2) == termination.round(2) and species == branch_tip_labels[x]
			net_diversification_candidates << branch_net_divs[x]
		end
	end
	if net_diversification_candidates.size != 1
		raise "ERROR: Branches could not unambiguously be assigned to net diversification rates!"
	else
		net_diversification = net_diversification_candidates[0]
	end
	regime = line_ary[4]
	mean_trait = (Math.log(line_ary[5].to_f) + Math.log(line_ary[6].to_f))/2.0
	regime_per_branch << regime
	mean_trait_per_branch << mean_trait
	branch_length_per_branch << origin - termination
	net_diversification_per_branch << net_diversification
end

# Prepare the output table string.
output_table = ""
output_table << "branch_length".rjust(20)
output_table << "regime".rjust(20)
output_table << "mean_trait".rjust(20)
output_table << "net_diversification".rjust(20)
output_table << "\n"
branch_length_per_branch.size.times do |x|
	output_table << "#{format("%.2f", branch_length_per_branch[x])}".rjust(20)
	output_table << "#{regime_per_branch[x]}".rjust(20)
	output_table << "#{format("%.4f", mean_trait_per_branch[x])}".rjust(20)
	output_table << "#{format("%.4f", net_diversification_per_branch[x])}".rjust(20)
	output_table << "\n"
end

# Write the output table to file.
output_table_file = File.open(output_table_file_name,"w")
output_table_file.write(output_table)
output_table_file.close

# Prepare the svg string.
x_dim = 600
y_dim = 580
margin_top = 20
margin_bottom = 60
margin_left = 80
margin_right = 20
min_x_value = 1
max_x_value = 5
min_y_value = 0.01
max_y_value = 0.09
max_radius = 10
min_radius = 2
frame_stroke_width = 1
circle_stroke_width = 1
stroke_color = "black"
max_sqrt_branch_length = 0
tick_length = 5
text_shift_vertical_x_ticks = 20
text_shift_vertical_y_ticks = 4
text_shift_horizontal = 10
font_size = 14
label_font_size = 16
text_shift_vertical_x_label = 40
text_shift_horizontal_y_label = 40
branch_length_per_branch.size.times {|x| max_sqrt_branch_length = Math.sqrt(branch_length_per_branch[x]) if Math.sqrt(branch_length_per_branch[x]) > max_sqrt_branch_length}
min_sqrt_branch_length = 1000
branch_length_per_branch.size.times {|x| min_sqrt_branch_length = Math.sqrt(branch_length_per_branch[x]) if Math.sqrt(branch_length_per_branch[x]) < min_sqrt_branch_length}
svg_string = ""
svg_string << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
svg_string << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
svg_string << "<svg width=\"#{x_dim}\" height=\"#{y_dim}\" viewBox=\"0 0 #{x_dim} #{y_dim}\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n"
svg_string << "  <style><![CDATA[\n"
svg_string << "    text{\n"
svg_string << "      font: #{font_size}px Helvetica, Arial, sans-serif;\n"
svg_string << "      }\n"
svg_string << "    ]]>\n"
svg_string << "  </style>\n"
svg_string << "  <rect x=\"#{margin_left}\" y=\"#{margin_top}\" width=\"#{x_dim-margin_left-margin_right}\" height=\"#{y_dim-margin_top-margin_bottom}\" fill=\"none\" stroke=\"#{stroke_color}\" stroke-width=\"#{frame_stroke_width}\"  />\n"
branch_length_per_branch.size.times do |x|
	circle_x = ((mean_trait_per_branch[x]-min_x_value)/(max_x_value-min_x_value)) * (x_dim - margin_left - margin_right) + margin_left
	circle_y = y_dim - (((net_diversification_per_branch[x]-min_y_value)/(max_y_value-min_y_value)) * (y_dim - margin_top - margin_bottom) + margin_bottom)
	# if circles should be according to branch length:
	#  circle_radius = (Math.sqrt(branch_length_per_branch[x])-min_sqrt_branch_length)/(max_sqrt_branch_length-min_sqrt_branch_length) * (max_radius-min_radius) + min_radius
	# else:
	if extant_per_branch[x] == true
		circle_radius = 10
	else
		circle_radius = 5
	end
	if regime_per_branch[x] == "1"
		circle_color = "#859900"
	elsif regime_per_branch[x] == "2"
		circle_color = "#cb4b16"
	elsif regime_per_branch[x] == "3"
		circle_color = "#2aa198"
	elsif regime_per_branch[x] == "6"
		circle_color = "#268bd2"
	elsif regime_per_branch[x] == "8"
		circle_color = "#6c71c4"
	elsif regime_per_branch[x] == "9"
		circle_color = "#d33682"
	else
		raise "ERROR: Unexpected regime : #{regime_per_branch[x]}!"
	end
	svg_string << "  <circle cx=\"#{circle_x}\" cy=\"#{circle_y}\" r=\"#{circle_radius}\" fill=\"#{circle_color}\" stroke=\"#{stroke_color}\" stroke-width=\"#{circle_stroke_width}\"  />\n"
end
svg_string << "  <line x1=\"#{margin_left}\" y1=\"#{y_dim-margin_bottom}\" x2=\"#{margin_left}\" y2=\"#{y_dim-margin_bottom+tick_length}\" stroke=\"#{stroke_color}\" stroke-width=\"#{frame_stroke_width}\"  />\n"
svg_string << "  <line x1=\"#{x_dim-margin_right}\" y1=\"#{y_dim-margin_bottom}\" x2=\"#{x_dim-margin_right}\" y2=\"#{y_dim-margin_bottom+tick_length}\" stroke=\"#{stroke_color}\" stroke-width=\"#{frame_stroke_width}\"  />\n"
svg_string << "  <line x1=\"#{margin_left-tick_length}\" y1=\"#{y_dim-margin_bottom}\" x2=\"#{margin_left}\" y2=\"#{y_dim-margin_bottom}\" stroke=\"#{stroke_color}\" stroke-width=\"#{frame_stroke_width}\"  />\n"
svg_string << "  <line x1=\"#{margin_left-tick_length}\" y1=\"#{margin_top}\" x2=\"#{margin_left}\" y2=\"#{margin_top}\" stroke=\"#{stroke_color}\" stroke-width=\"#{frame_stroke_width}\"  />\n"
svg_string << "  <text text-anchor=\"middle\" x=\"#{margin_left}\" y=\"#{y_dim-margin_bottom+text_shift_vertical_x_ticks}\" fill=\"#{stroke_color}\">#{min_x_value}</text>\n"
svg_string << "  <text text-anchor=\"middle\" x=\"#{x_dim-margin_right}\" y=\"#{y_dim-margin_bottom+text_shift_vertical_x_ticks}\" fill=\"#{stroke_color}\">#{max_x_value}</text>\n"
svg_string << "  <text text-anchor=\"end\" x=\"#{margin_left-text_shift_horizontal}\" y=\"#{y_dim-margin_bottom+text_shift_vertical_y_ticks}\" fill=\"#{stroke_color}\">#{min_y_value}</text>\n"
svg_string << "  <text text-anchor=\"end\" x=\"#{margin_left-text_shift_horizontal}\" y=\"#{margin_top+text_shift_vertical_y_ticks}\" fill=\"#{stroke_color}\">#{max_y_value}</text>\n"
svg_string << "  <text text-anchor=\"middle\" x=\"#{margin_left+(x_dim-margin_left-margin_right)/2.0}\" y=\"#{y_dim-margin_bottom+text_shift_vertical_x_label}\" style=\"font-size:#{label_font_size}px;\" fill=\"#{stroke_color}\">log(MHCI copy number)</text>\n"
svg_string << "  <text text-anchor=\"middle\" x=\"#{margin_left-text_shift_horizontal_y_label}\" y=\"#{margin_top+(y_dim-margin_bottom-margin_top)/2.0}\" style=\"font-size:#{label_font_size}px;\" transform=\"rotate(-90,#{margin_left-text_shift_horizontal_y_label},#{margin_top+(y_dim-margin_bottom-margin_top)/2.0})\" fill=\"#{stroke_color}\">Net diversification rate</text>\n"
svg_string << "</svg>\n"

# Write the svg string to file.
output_svg_file = File.open(output_svg_file_name,"w")
output_svg_file.write(svg_string)
