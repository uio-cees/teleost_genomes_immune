queries_file_name=$1
subjects_file_name=$2
alignment_dir_out=${3}

# Make the alignment output directory if it doesn't exist yet.
mkdir -p $alignment_dir_out

# Run find_orthologs.py for each marker, with each of the 76 subjects.
python3 resources/find_orthologs.py -t --refine -s 1 $queries_file_name $subjects_file_name
cp *.fasta $alignment_dir_out

