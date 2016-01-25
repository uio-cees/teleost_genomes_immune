# Michael Matschiner, 2015-07-09

# Get the command line arguments.
tree_file_name=$1
table_file_name=$2
log_copy_number_phenogram_file_name=$3

# Get the directory names of the output files.
log_copy_number_phenogram_file_base_name=`basename ${log_copy_number_phenogram_file_name}`
log_copy_number_phenogram_dir=${log_copy_number_phenogram_file_name%${log_copy_number_phenogram_file_base_name}}

# Make the output directories if they don't exist yet.
mkdir -p ${log_copy_number_phenogram_dir}

# Produce the phenograms for copy number and log copy number.
rscript produce_phenogram.r ${tree_file_name} ${table_file_name} ${log_copy_number_phenogram_file_name}