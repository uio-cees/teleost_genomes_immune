# Michael Matschiner, 2015-07-21

# Get the command line arguments.
tree_file_name=$1
order_file_name=$2
slouch_input_dir=$3
slouch_output_dir=$4

# Run ancestral state reconstruction for all text files in the input directory.
for i in ${slouch_input_dir}/*.txt
do
	regime_code_w_ext=`basename ${i}`
	regime_code=${regime_code_w_ext%.txt}
	ruby calculate_ancestral_states.rb ${tree_file_name} ${slouch_input_dir}/${regime_code_w_ext} ${order_file_name} ${slouch_output_dir}/${regime_code}.results.txt ${slouch_output_dir}/${regime_code}.3d.svg ${slouch_output_dir}/${regime_code}.3d.txt
done