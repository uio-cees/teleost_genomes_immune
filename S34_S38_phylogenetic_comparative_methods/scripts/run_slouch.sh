# Michael Matschiner, 2015-07-17

# Get the command line argument.
slouch_input_dir=${1}
slouch_output_dir=${2}

# Make the output directory if it does not exist yet.
mkdir -p ${slouch_output_dir}

# For each text file in the slouch input directory, run slouch.
for i in ${slouch_input_dir}/*.txt
do
	regime_code_w_ext=`basename ${i}`
	regime_code=${regime_code_w_ext%.txt}
	echo ${regime_code}
	slouch_table_file_name=${i}
	slouch_result_file_name=${2}/${regime_code}.results.txt
	slouch_plot_file_name=${2}/${regime_code}.pdf
	# Run slouch.
	rscript run_slouch.r ${slouch_table_file_name} ${slouch_plot_file_name} &> ${slouch_result_file_name}
done
