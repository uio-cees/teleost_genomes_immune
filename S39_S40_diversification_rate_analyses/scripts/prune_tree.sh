# Michael Matschiner, 2015-07-28

# Get the command line arguments.
input_table_file_name=$1
trait_table_file_name=$2
input_tree_file_name=$3
output_tree_file_name_w_relative_path=$4

# Make the output directories if it doesn't exist yet.
output_tree_file_name=`basename ${output_tree_file_name_w_relative_path}`
output_tree_directory=${output_tree_file_name_w_relative_path%${output_tree_file_name}}
mkdir -p ${output_tree_directory}

# Use prepare_for_pruning.rb to rename species ids in the tree and to prepare a table for tree pruning.
ruby prepare_for_pruning.rb ${input_table_file_name} ${trait_table_file_name} tmp.txt ${input_tree_file_name} tmp.tre

# Use prune_trees.r to prune the tree so that only one tip is left per clade.
rscript prune_tree.r tmp.tre ${output_tree_file_name_w_relative_path} tmp.txt

# Clean up.
rm tmp.tre
rm tmp.txt