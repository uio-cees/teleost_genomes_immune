# Michael Matschiner, 2015-07-15

# Get the command line arguments.
tree_file_name=$1
data_file_name=$2
plot_file_name=$3
results_file_name=$4

# Make the output plot directory if it does not exist yet.
plot_file_base_name=`basename ${plot_file_name}`
plot_dir=${plot_file_name%${plot_file_base_name}}
mkdir -p ${plot_dir}

# Make the results directory if it does not exist yet.
results_file_base_name=`basename ${results_file_name}`
results_dir=${results_file_name%${results_file_base_name}}
mkdir -p ${results_dir}

# Run trait evolution analyses.
rscript analyse_trait_evolution.r ${tree_file_name} ${data_file_name} ${plot_file_name} > ${results_file_name}