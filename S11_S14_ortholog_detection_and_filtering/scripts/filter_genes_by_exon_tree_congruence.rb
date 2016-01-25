## Michael Matschiner 2015-02-28.

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
number_of_cpus = ARGV[2].to_i
resources_directory = "./resources"

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
    Dir.mkdir(alignment_directory_out)
end

# Get the names of the directories with exon alignments for each gene.
alignment_directory_entries = Dir.entries(alignment_directory_in)
gene_alignment_directories = []
alignment_directory_entries.each {|e| gene_alignment_directories << e if e[0..6] == "ENSDARG"}

# Memorize the current directory
home = Dir.pwd

# Statistics.
number_of_genes = 0
total_alignment_length = 0

# Initiate an array for the ids and lengths of all filtered exons.
filtered_exon_ids = []
filtered_exon_lengths = []

# For each of the concaterpillar directories, do the following.
gene_alignment_directories.each do |dir|

    # Make sure this directory has not been analysed before.
    unless Dir.entries("#{alignment_directory_in}/#{dir}").include?("results.ccp")

        # Copy resources to that directory.
        system("cp #{resources_directory}/*.py* #{alignment_directory_in}/#{dir}")
        system("cp #{resources_directory}/raxmlHPC #{alignment_directory_in}/#{dir}")

        # Move into the analysis directory.
        Dir.chdir("#{alignment_directory_in}/#{dir}")

        # Start concaterpillar.
        system("python concaterpillar.py -t -c #{number_of_cpus} -m GTR")

        # Move back up.
        Dir.chdir(home)

        # Cleanup resources from that directory.
        system("rm #{alignment_directory_in}/#{dir}/*py*")
        system("rm #{alignment_directory_in}/#{dir}/raxmlHPC")

    end

    # Read the concaterpillar results file.
    results_file = File.open("#{alignment_directory_in}/#{dir}/results.ccp")
    results_lines = results_file.readlines

    # Parse the concaterpillar results file.
    congruent_sets = []
    in_sets = false
    results_lines.each do |l|
        if l.strip == "Finished concatenating sets.  The following files contain the final sets:"
            in_sets = true
        elsif in_sets == true
            congruent_sets << l.split(":")[1].strip[1..-2].gsub("'","").gsub(" ","").gsub(".seq","").split(",")
        end
    end

    # Sort the congruent sets by their size.
    if congruent_sets.size > 1
        sets_sorted = false
        while sets_sorted == false
            sets_sorted = true
            (congruent_sets.size-1).times do |x|
                if congruent_sets[x].size < congruent_sets[x+1].size
                    congruent_sets[x],congruent_sets[x+1] = congruent_sets[x+1],congruent_sets[x]
                    sets_sorted = false
                end
            end
        end
    end

    # Only choose the largest set if there is only one of its size, and if its size is at least three.
    selected_set = []
    if congruent_sets.size == 1 or congruent_sets[0].size > congruent_sets[1].size
        if congruent_sets[0].size >= 3
            congruent_sets[0].each do |i|
                selected_set << i
            end
        end
    end

    unless selected_set == []

        # Get the nexus file names for the selected set.
        nexus_file_names = []
        selected_set.each do |i|
            nexus_file_names << "#{i}.seq"
        end
        nexus_file_names.sort!
        
        # Get the ids and seqs from each alignment.
        nexus_ids_per_alignment = []
        nexus_seqs_per_alignment = []
        nexus_file_names.each do |f|
            nexus_file = File.open("#{alignment_directory_in}/#{dir}/#{f}")
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
            nexus_ids_per_alignment << nexus_ids_in_this_alignment
            nexus_seqs_per_alignment << nexus_seqs_in_this_alignment
            filtered_exon_ids << f.chomp(".seq")
            filtered_exon_lengths << nexus_seqs_in_this_alignment[0].size
        end

        # Make sure the IDs are the same in all alignments.
        1.upto(nexus_ids_per_alignment.size-1) do |x|
            if nexus_ids_per_alignment[0] != nexus_ids_per_alignment[x]
                raise "IDs differ between alignments!"
            end
        end

        # Make sure the same number of sequences is present in each set.
        1.upto(nexus_seqs_per_alignment.size-1) do |x|
            if nexus_seqs_per_alignment[0].size != nexus_seqs_per_alignment[x].size
                raise "Number of sequences differs between alignments!"
            end
        end

        # Produce the concatenated sequences.
        concatenated_seqs = []
        nexus_ids_per_alignment[0].size.times do |x|
            concatenated_seq_for_this_taxon = ""
            nexus_seqs_per_alignment.each do |seqs|
                concatenated_seq_for_this_taxon << seqs[x]
            end
            concatenated_seqs << concatenated_seq_for_this_taxon
        end
        
        # Get stats.
        number_of_genes += 1
        total_alignment_length += concatenated_seqs[0].size

        # Prepare the string for a new nexus file.
        Dir.mkdir("#{alignment_directory_out}/#{dir}") unless Dir.exists?("#{alignment_directory_out}/#{dir}")
        new_nexus_string = "#nexus\n"
        new_nexus_string << "\n"
        new_nexus_string << "begin data;\n"
        new_nexus_string << "dimensions  ntax=#{nexus_ids_per_alignment[0].size} nchar=#{concatenated_seqs[0].size};\n"
        new_nexus_string << "format datatype=DNA gap=- missing=?;\n"
        new_nexus_string << "matrix\n"
        nexus_ids_per_alignment[0].size.times do |x|
            new_nexus_string << "#{nexus_ids_per_alignment[0][x].ljust(12)}#{concatenated_seqs[x]}\n"
        end
        new_nexus_string << ";\n"
        new_nexus_string << "end;\n"
        
        # Write the new nexus file.
        new_nexus_file_name = "#{dir}.nex"
        new_nexus_file = File.open("#{alignment_directory_out}/#{dir}/#{new_nexus_file_name}","w")
        new_nexus_file.write(new_nexus_string)

    end

end

# Edit the info file ../../filtered_nuclear_exons.txt, removing exons that were filtered due to incongruence.
info_file_name = "../../info/filtered_nuclear_exons.txt"
info_file = File.open(info_file_name)
info_lines = info_file.readlines
info_file.close
new_info_lines = [info_lines[0]]
info_lines.each do |l|
    use_line = false
    len = 0
    filtered_exon_ids.size.times do |x|
        if l.split[0] == filtered_exon_ids[x]
            use_line = true
            len = filtered_exon_lengths[x]
            break
        end
    end
    if use_line
        line = l
        line_ary = line.split
        line_ary[4] = len.to_s
        line = ""
        line_ary.each do |i|
            line << "#{i.strip}\t"
        end
        line.chomp("\t")
        line << "\n"
        new_info_lines << line
    end
end
new_info_string = ""
new_info_lines.each {|l| new_info_string << l}
info_file = File.open(info_file_name,"w")
info_file.write(new_info_string)
info_file.close

# Report stats.
puts "Number of genes: #{number_of_genes}"
puts "Total alignment length: #{total_alignment_length}"
