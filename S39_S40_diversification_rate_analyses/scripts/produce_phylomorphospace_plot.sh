# Michael Matschiner, 2015-09-17

tree_file_name=$1
table_file_name=$2
plot_file_name=$3

# Make the output plot directory if it does not exist yet.
plot_file_base_name=`basename ${plot_file_name}`
plot_dir=${plot_file_name%${plot_file_base_name}}
mkdir -p ${plot_dir}

rscript produce_phylomorphospace_plot.r ${tree_file_name} ${table_file_name} ${plot_file_name}