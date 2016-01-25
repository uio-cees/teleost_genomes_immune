## Michael Matschiner 2015-03-06.

require "./resources/array_stats.rb"

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
    Dir.mkdir(alignment_directory_out)
end

# Get the names of the directories with exon alignments for each gene.
alignment_directory_entries = Dir.entries(alignment_directory_in)
gene_alignment_directories = []
alignment_directory_entries.each {|e| gene_alignment_directories << e if e[0..6] == "ENSDARG"}

# Read the nuclear_queries_exons.txt file from the info directory.
filtered_nuclear_exons_info_file_name = "../../info/filtered_nuclear_exons.txt"
filtered_nuclear_exons_info_file = File.open(filtered_nuclear_exons_info_file_name)
filtered_nuclear_exons_info_lines = filtered_nuclear_exons_info_file.readlines
selected_nuclear_exons_info_lines = []

# Prepare strings for concatenated nexus files.
cp1_concatenated_seqs = []
cp2_concatenated_seqs = []
concatenated_seqs = []
concatenated_loci_count = 0

# For each of the alignment directories, do the following.
is_first_dir = true
nexus_ids_in_first_alignment = []
gene_alignment_directories.each do |dir|

    # Read the BEAST log file.
    log_file = File.open("#{alignment_directory_in}/#{dir}/#{dir}.log")
    log_lines = log_file.readlines
    coefficients_of_variation = []
    log_lines.each do |l|
        coefficients_of_variation << l.split[12].to_f if l[0] != "#"
    end
    coefficients_of_variation = coefficients_of_variation[1..-1]
    number_of_burnin_samples = (coefficients_of_variation.size)/5
    coefficients_of_variation_wo_burning = coefficients_of_variation[number_of_burnin_samples+1..-1]

    if coefficients_of_variation_wo_burning.mean < 1 and coefficients_of_variation_wo_burning.hpd_upper(0.95) < 1.2

        # Get the ids and seqs from the Nexus format alignment.
        nexus_file = File.open("#{alignment_directory_in}/#{dir}/#{dir}.nex")
        nexus_lines = nexus_file.readlines
        nexus_ids_in_this_alignment = []
        nexus_seqs_in_this_alignment = []
        in_matrix = false
        nexus_lines.each do |l|
            if l.strip.downcase == "matrix"
                in_matrix = true
            elsif l.strip == ";"
                in_matrix = false
            elsif l.strip != "" and in_matrix
                nexus_ids_in_this_alignment << l.split[0].strip
                nexus_seqs_in_this_alignment << l.split[1].strip
            end
        end

        # Store the nexus ids of the first alignment and prepare the concatenated alignments.
        if is_first_dir
            nexus_ids_in_first_alignment = nexus_ids_in_this_alignment
            nexus_ids_in_first_alignment.size.times do
                cp1_concatenated_seqs << ""
                cp2_concatenated_seqs << ""
                concatenated_seqs << ""
            end
        else
            raise "Nexus is differ between alignments!" if nexus_ids_in_first_alignment != nexus_ids_in_this_alignment
        end

        # Add sequences of this alignment to the concatenated alignment.
        nexus_seqs_in_this_alignment.size.times do |x|
            cp1_seq = ""
            cp2_seq = ""
            nexus_seqs_in_this_alignment[x].size.times do |y|
                if (y/2)*2 == y
                    cp1_seq << nexus_seqs_in_this_alignment[x][y]
                else
                    cp2_seq << nexus_seqs_in_this_alignment[x][y]
                end
            end
            cp1_concatenated_seqs[x] << cp1_seq
            cp2_concatenated_seqs[x] << cp2_seq
            concatenated_seqs[x] << nexus_seqs_in_this_alignment[x]
        end
        concatenated_loci_count += 1

        # Prepare the string for a new nexus file.
        new_nexus_string = "#nexus\n"
        new_nexus_string << "\n"
        new_nexus_string << "begin data;\n"
        new_nexus_string << "dimensions  ntax=#{nexus_ids_in_this_alignment.size} nchar=#{nexus_seqs_in_this_alignment[0].size};\n"
        new_nexus_string << "format datatype=DNA gap=- missing=?;\n"
        new_nexus_string << "matrix\n"
        nexus_seqs_in_this_alignment.size.times do |x|
            new_nexus_string << "#{nexus_ids_in_this_alignment[x].ljust(12)}#{nexus_seqs_in_this_alignment[x]}\n"
        end
        new_nexus_string << ";\n"
        new_nexus_string << "end;\n"

        # Write the new nexus file.
        new_nexus_file_name = "#{dir}.seq"
        new_nexus_file = File.open("#{alignment_directory_out}/#{new_nexus_file_name}","w")
        new_nexus_file.write(new_nexus_string)

        # Add to the selected nuclear exons info lines.
        filtered_nuclear_exons_info_lines.each do |l|
            selected_nuclear_exons_info_lines << l if l.include?(dir)
        end

    end
    is_first_dir = false

end

# Prepare the string for a concatenated phylip file.
phylip_string = "#{concatenated_seqs.size} #{concatenated_seqs[0].size}\n"
concatenated_seqs.size.times do |x|
    phylip_string << "#{nexus_ids_in_first_alignment[x].ljust(10)} #{concatenated_seqs[x]}\n"
end

# Write the concatenated phylip file.
phylip_file_name = "concatenated_#{concatenated_loci_count}.phy"
phylip_file = File.open("#{alignment_directory_out}/#{phylip_file_name}","w")
phylip_file.write(phylip_string)

# Prepare strings for concatenated nexus files of the first and second codon position.
cp1_nexus_string = "#nexus\n"
cp2_nexus_string = "#nexus\n"
cp1_nexus_string = "\n"
cp2_nexus_string = "\n"
cp1_nexus_string << "begin data;\n"
cp2_nexus_string << "begin data;\n"
cp1_nexus_string << "dimensions  ntax=#{cp1_concatenated_seqs.size} nchar=#{cp1_concatenated_seqs[0].size};\n"
cp2_nexus_string << "dimensions  ntax=#{cp2_concatenated_seqs.size} nchar=#{cp2_concatenated_seqs[0].size};\n"
cp1_nexus_string << "format datatype=DNA gap=- missing=?;\n"
cp2_nexus_string << "format datatype=DNA gap=- missing=?;\n"
cp1_nexus_string << "matrix\n"
cp2_nexus_string << "matrix\n"
cp1_concatenated_seqs.size.times do |x|
    cp1_nexus_string << "#{nexus_ids_in_first_alignment[x].ljust(12)}#{cp1_concatenated_seqs[x]}\n"
end
cp2_concatenated_seqs.size.times do |x|
    cp2_nexus_string << "#{nexus_ids_in_first_alignment[x].ljust(12)}#{cp2_concatenated_seqs[x]}\n"
end
cp1_concatenated_seqs << ";\n"
cp2_concatenated_seqs << ";\n"
cp1_concatenated_seqs << "end;\n"
cp2_concatenated_seqs << "end;\n"

# Write the concatenated file.
cp1_nexus_file_name = "concatenated_#{concatenated_loci_count}_cp1.nex"
cp2_nexus_file_name = "concatenated_#{concatenated_loci_count}_cp2.nex"
cp1_nexus_file = File.open("#{alignment_directory_out}/#{cp1_nexus_file_name}","w")
cp2_nexus_file = File.open("#{alignment_directory_out}/#{cp2_nexus_file_name}","w")
cp1_nexus_file.write(cp1_nexus_string)
cp2_nexus_file.write(cp2_nexus_string)

# Prepare the string for the selected nuclear exons file.
selected_nuclear_exons_info_string = "exon_id\tgene_id\ttranscript_id\ttranslation\tlength\tthreshold_bitscore\tcorrect_bitscores\tincorrect_bitscores\tmin(correct_bitscores)-max(incorrect_bitscores)\n"
selected_nuclear_exons_info_lines.each {|l| selected_nuclear_exons_info_string << l}

# Write the selected nuclear exons file to the info directory.
selected_nuclear_exons_info_file_name = "../../info/selected_nuclear_exons.txt"
selected_nuclear_exons_info_file = File.open(selected_nuclear_exons_info_file_name,"w")
selected_nuclear_exons_info_file.write(selected_nuclear_exons_info_string)
selected_nuclear_exons_info_file.close
