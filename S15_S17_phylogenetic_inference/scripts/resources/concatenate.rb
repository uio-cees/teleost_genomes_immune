## this script takes a batch of (simple) nexus files and concatenates them into a single file of nexus, phylip, or fasta format.

# Define the valid formats.
valid_formats = ["nexus","phylip","fasta"]
format_raise_string = "A file format ("
valid_formats.each {|f| format_raise_string << "#{f}, "}
format_raise_string.chomp!(", ")
format_raise_string << ") must be specified with the option \"-f\"!"

# Define the help string.
help_string = ""
help_string << "\n"
help_string << "concatenate.rb\n"
help_string << "\n"
help_string << "Available options:\n"
help_string << "  Option   Value                    Comment\n"
help_string << "  -i       input file names       | Full or relative path plus file name (wildcards accepted)\n"
help_string << "  -o       output file name       | Full or relative path plus file name\n"
help_string << "  -f       format                 | Output format, can be "
valid_formats[0..-2].each {|f| help_string << "#{f}, "}
help_string << "or #{valid_formats[-1]}\n"
help_string << "  -c       -                      | Writes a cfg file for PartitionFinder in the output directory\n"
help_string << "\n"

# Read the arguments.
if ARGV == [] or ["-h","--help","-help"].include?(ARGV[0].downcase)
  puts help_string
  exit
end

# Read the option of PartitionFinder cfg files.
cfg_file = false
if ARGV.include?("-c")
  cfg_file = true
  ARGV[ARGV.index("-c")] = nil
  ARGV.compact!
end

# Read the specified format and remove the -f option from the ARGV array.
if ARGV.include?("-f")
  raise format_raise_string if ARGV[ARGV.index("-f")+1] == nil
  format = ARGV[ARGV.index("-f")+1].downcase
  raise format_raise_string unless valid_formats.include?(format)
else
  raise format_raise_string
end
ARGV[ARGV.index("-f")+1] = nil
ARGV[ARGV.index("-f")] = nil
ARGV.compact!

# Read the specified output file name.
if ARGV.include?("-o")
  raise "Please specify an output file name with the option \"-o\"!" if ARGV[ARGV.index("-o")+1] == nil
  output_file_name = ARGV[ARGV.index("-o")+1]
end
ARGV[ARGV.index("-o")+1] = nil
ARGV[ARGV.index("-o")] = nil
ARGV.compact!

# Read the specified input file names.
if ARGV.include?("-i")
  input_file_names = ARGV[ARGV.index("-i")+1..-1]
else
  raise "At least one input file name must be given with option \"-i\"!"
end

# Initiate arrays for output ids and output sequences.
out_ids = []
out_seqs = []

# Initiate arrays for the ids and lengths of the markers.
marker_ids = []
marker_lengths = []

input_file_names.each do |a|
  puts "Reading #{a}..."
  marker_id = a.split("/").last.chomp(".nex")
  matrix = false
  nexusIds = []
  nexusSequences = []
  nchar = 0
  ntax = 0
  taxCount = 0
  nexusLines = IO.readlines(a)
  raise "File is not in NEXUS format: #{a}" unless nexusLines[0].strip.downcase == "#nexus"
  nexusLines.size.times do |l|
    nexusLines[l].strip!
    if nexusLines[l].strip[0..9].downcase == "dimensions"
      nexus_line = nexusLines[l].strip.gsub(" = ","=")
      dimensionsAry = nexus_line.split(" ")
      raise "Dimensions could not be read properly: #{l}" unless dimensionsAry.size == 3
      ntax = dimensionsAry[1].split("=")[1].to_i
      nchar = dimensionsAry[2].split("=")[1].to_i
    end
    matrix = true if nexusLines[l-1].strip[0..5].downcase == "matrix"
    matrix = false if nexusLines[l].strip == ";"
    if matrix == true
      unless nexusLines[l].strip == ""
        raise "Taxon name in quotes detected in file #{a}" if nexusLines[l].include?("'")
        linesAry = nexusLines[l].split(" ")
        raise "Matrix line could not be read correctly: #{nexusLines[l]}: linesAry.size = #{linesAry.size}" unless linesAry.size == 2
        nexusIds << linesAry[0].strip
        nexusSequences << linesAry[1].strip
      end
    end
  end
  (nexusSequences.size-1).times do |n|
    raise "Two sequences have different lengths!" if nexusSequences[n].length != nexusSequences[n+1].length
  end
  if out_ids == []
    out_ids = nexusIds
    out_ids.size.times {out_seqs << ""}
  else
    raise "ID order differs!" if out_ids != nexusIds
  end
  nexusSequences.size.times do |x|
    out_seqs[x] << nexusSequences[x]
  end
  marker_ids << marker_id
  marker_lengths << nchar
end

# Determine the maximum id length.
max_id_length = out_ids[0].size
1.upto(out_ids.size-1) do |x|
  max_id_length = out_ids[x].size if out_ids[x].size > max_id_length
end

out_string = ""
if format == "phylip"
  out_string << "#{out_ids.size} #{out_seqs[0].length}\n"
  out_ids.size.times do |x|
    out_string << "#{out_ids[x].ljust(max_id_length)} #{out_seqs[x]}\n"
  end
elsif format == "nexus"
  out_string << "#nexus\n"
  out_string << "\n"
  out_string << "begin data;\n"
  out_string << "  dimensions ntax=#{out_ids.size} nchar=#{out_seqs[0].length};\n"
  out_string << "  format datatype=dna gap=- missing=?;\n"
  out_string << "  matrix\n"
  out_ids.size.times do |x|
    out_string << "  #{out_ids[x].ljust(max_id_length)} #{out_seqs[x]}\n"
  end
  out_string << "  ;\n"
  out_string << "end;\n"
elsif format == "fasta"
  out_ids.size.times do |x|
    out_string << ">#{out_ids[x]}\n"
    seq_to_write = out_seqs[x]
    while seq_to_write.size > 0 do
      out_string << "#{seq_to_write.slice!(0..59)}\n"
    end
  end
else
  raise "Unrecognized format #{format}!"
end
File.new(output_file_name,"w").puts out_string
puts "Wrote file #{output_file_name}."

out_seqs.size.times do |c|
  warn "No sequence information available for #{out_ids[c]}!" if out_seqs[c].gsub("-","").length == 0
end

# Write cfg file for PartitionFinder.
if cfg_file
  cfg_out_string = "## ALIGNMENT FILE ##\n"
  cfg_out_string << "alignment = #{output_file_name.split("/").last};\n"
  cfg_out_string << "\n"
  cfg_out_string << "## BRANCHLENGTHS: linked | unlinked ##\n"
  cfg_out_string << "branchlengths = unlinked;\n"
  cfg_out_string << "\n"
  cfg_out_string << "## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | beast | <list> ##\n"
  cfg_out_string << "##              for PartitionFinderProtein: all_protein | <list> ##\n"
  cfg_out_string << "models = beast;\n"
  cfg_out_string << "\n"
  cfg_out_string << "# MODEL SELECCTION: AIC | AICc | BIC #\n"
  cfg_out_string << "model_selection = BIC;\n"
  cfg_out_string << "\n"
  cfg_out_string << "## DATA BLOCKS: see manual for how to define ##\n"
  cfg_out_string << "[data_blocks]\n"
  current_pos = 0
  marker_ids.size.times do |m|
    cfg_out_string << "#{marker_ids[m]} = #{current_pos+1}-#{current_pos+marker_lengths[m]};\n"
    current_pos += marker_lengths[m]
  end
  cfg_out_string << "\n"
  cfg_out_string << "## SCHEMES, search: all | user | greedy ##\n"
  cfg_out_string << "[schemes]\n"
  cfg_out_string << "search = greedy;\n"
  cfg_out_string << "\n"
  cfg_out_string << "#user schemes go here if search=user. See manual for how to define.#\n"
  cfg_out_name = "#{output_file_name.chomp(output_file_name.split("/")[-1])}partition_finder.cfg"
  File.new(cfg_out_name,"w").puts cfg_out_string
end
