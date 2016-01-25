# This script uses tables like those produces by BLAST with the outfmt 6 option, and writes fasta files from the table's content.
# Run as
# ruby table_to_fasta.rb -f inputfilename -o outputfilename -id column_in_table_that_contains_the_ids -seq column_in_table_that_contains_the_sequences

if ARGV.include?("-f")
	table_file_name = ARGV[ARGV.index("-f")+1]
else
	raise "The tabular file name should be specified with the '-f' option!"
end

if ARGV.include?("-o")
	output_file_name = ARGV[ARGV.index("-o")+1]
else
	raise "The output file (in fasta format) name should be specified with the '-o' option!"
end

if ARGV.include?("-id")
	id_col = ARGV[ARGV.index("-id")+1].to_i - 1
else
	raise "The column of the id should be specified with the '-id' option!"
end

if ARGV.include?("-seq")
	seq_col = ARGV[ARGV.index("-seq")+1].to_i - 1
else
	raise "The column of the sequence should be specified with the '-seq' option!"
end

table_file = File.open(table_file_name)
table_lines = table_file.readlines

ids = []
seqs = []

table_lines.each do |l|
	line = l.strip
	unless line == ""
		line_ary = line.split
		ids << line_ary[id_col]
		seqs << line_ary[seq_col].gsub("-","")
	end
end

output_string = ""
ids.size.times do |x|
	output_string << ">#{ids[x]}\n"
	while seqs[x].size > 0 do 
		output_string << "#{seqs[x].slice!(0..59)}\n"
	end
	output_string << "\n"
end

output_file = File.new(output_file_name,"w")
output_file.write(output_string)
output_file.close

