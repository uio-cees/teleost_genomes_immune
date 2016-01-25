# Michael Matschiner, 2015-09-22

# Get the command line argument.
diversitree_dir=$1

# Run the diversitree script for each analysis directory.
for i in ${diversitree_dir}/t??_r???
do
	analysis_dir=${i}
	if [ ! -f ${analysis_dir}/results.txt ];
	then
		rscript run_diversitree.r ${analysis_dir}/unresolved.tre ${analysis_dir}/resolved.txt ${analysis_dir}/unresolved.txt > ${analysis_dir}/results.txt
	fi
done