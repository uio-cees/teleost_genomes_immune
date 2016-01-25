# Michael Matschiner, 2015-07-09

# Load FileUtils (needed to recursively make directories).
require 'fileutils'

# Define a function to calculate the mean.
module Enumerable
	def sum
		self.inject(0){|accum, i| accum + i }
	end
	def mean
		if self.length == 0
			nil
		else
			self.sum/self.length.to_f
		end
	end
	def sample_variance
		if self.length == 0
			nil
		else
			m = self.mean
			sum = self.inject(0){|accum, i| accum +(i-m)**2 }
			sum/(self.length - 1).to_f
		end
	end
	def standard_deviation
		if self.length == 0
			nil
		else
			return Math.sqrt(self.sample_variance)
		end
	end
	def standard_error
		if self.length == 0
			nil
		else
			return self.standard_deviation/Math.sqrt(self.length)
		end
	end
end

# Get the command line arguments.
u_lineage_table_file_name = ARGV[0]
z_lineage_table_file_name = ARGV[1]
species_table_file_name = ARGV[2]
log_copy_number_table_file_name = ARGV[3]
taxa_list_file_name = ARGV[4]

# Make the output directory if it does not exist yet.
if log_copy_number_table_file_name.include?("/")
	output_tree_directory = log_copy_number_table_file_name.chomp(log_copy_number_table_file_name.split("/")[-1])
	FileUtils.mkdir_p(output_tree_directory)
end

# Read the species table.
species_table_file = File.open(species_table_file_name)
species_table_lines = species_table_file.readlines
species_names = []
species_ids = []
species_table_lines.each do |l|
	unless l.strip == ""
		line_ary = l.split("\t")
		species_names << line_ary[0].strip
		species_ids << line_ary[1].strip
	end
end

# Read the U-lineage and Z-lineage tables.
u_lineage_table_file = File.open(u_lineage_table_file_name)
u_lineage_table_lines = u_lineage_table_file.readlines
z_lineage_table_file = File.open(z_lineage_table_file_name)
z_lineage_table_lines = z_lineage_table_file.readlines

# Get the mean from both tables, for each species id.
log_copy_number_table_string = ""
species_ids.size.times do |x|

	# Get the column index for this species in both tables.
	u_lineage_table_header = u_lineage_table_lines[0].split("\t")
	u_lineage_table_column_names = []
	u_lineage_table_header.each {|i| u_lineage_table_column_names << i.strip}
	u_lineage_index = u_lineage_table_column_names.index(species_names[x])
	z_lineage_table_header = z_lineage_table_lines[0].split("\t")
	z_lineage_table_column_names = []
	z_lineage_table_header.each {|i| z_lineage_table_column_names << i.strip}
	z_lineage_index = z_lineage_table_column_names.index(species_names[x])

	# Get all values for this species.
	u_lineage_values_for_this_species = []
	z_lineage_values_for_this_species = []
	u_lineage_table_lines[1..-1].each do |l|
		line_ary = l.split("\t")
		u_lineage_values_for_this_species << line_ary[u_lineage_index].to_f
	end
	z_lineage_table_lines[1..-1].each do |l|
		line_ary = l.split("\t")
		z_lineage_values_for_this_species << line_ary[z_lineage_index].to_f
	end
	log_both_values_for_this_species = []
	u_lineage_values_for_this_species.size.times do |z|
		log_both_values_for_this_species << Math.log(u_lineage_values_for_this_species[z] + z_lineage_values_for_this_species[z])
	end
	log_copy_number_for_this_species_mean = log_both_values_for_this_species.mean
	log_copy_number_for_this_species_stderr = log_both_values_for_this_species.standard_error
	log_copy_number_table_string << "#{species_ids[x]}\t#{log_copy_number_for_this_species_mean}\t#{log_copy_number_for_this_species_stderr}\n"
end

# Write a table for the copy numbers.
log_copy_number_table_file = File.open(log_copy_number_table_file_name,"w")
log_copy_number_table_file.write(log_copy_number_table_string)
log_copy_number_table_file.close

# Prepare a list of all taxa.
taxa_list_string = ""
species_ids.each {|i| taxa_list_string << "#{i}\n"}

# Write the taxa list to file.
taxa_list_file = File.open(taxa_list_file_name,"w")
taxa_list_file.write(taxa_list_string)
taxa_list_file.close
