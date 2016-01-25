# Michael Matschiner, 2015-07-07

# Define names and make directories.
astral_input_directory="../analysis/nuclear_single/astral/input"
astral_tree_file_name="combined.trees"
alignment_directory="../data/nuclear"
mkdir -p ${astral_input_directory}/beast2/mcc
mkdir -p ${astral_input_directory}/beast2/sample
mkdir -p ${astral_input_directory}/beast2/context
mkdir -p ${astral_input_directory}/beast2/extra

# Run treeannotator for the posterior tree distributions of the BEAST analysis done during ortholog detection.
# for i in ${alignment_directory}/*.nex
# do
# 	gene_id_w_ext=`basename ${i}`
# 	gene_id=${gene_id_w_ext%.nex}
# 	trees_file_directory="../../ortholog_identification/analysis/alignments/nuclear/09/${gene_id}"
# 	trees_file_name="${gene_id}.trees"
# 	mcc_tree_file_name=${trees_file_name%.trees}.tre
# 	sample_tree_file_name=${trees_file_name%.trees}.trees
# 	resources/treeannotator/bin/treeannotator -burnin 20 -heights mean ${trees_file_directory}/${trees_file_name} ${astral_input_directory}/beast2/mcc/${mcc_tree_file_name}
# 	resources/logcombiner.py --remove-comments -b 20 -n 100 ${trees_file_directory}/${trees_file_name} tmp.trees
# 	cat tmp.trees | grep "tree " | cut -d "=" -f 2 > ${astral_input_directory}/beast2/sample/${sample_tree_file_name}
# 	rm tmp.trees
# done

# Prepare a single file with all MCC trees.
for i in ${astral_input_directory}/beast2/mcc/ENSDARG*.tre
do
	cat ${i} | tail -r | tail -n +3 | tail -r > head.txt
	cat ${i} | grep "tree " | tail -r | tail -n +1 | tail -r >> body.txt
	cat ${i} | tail -n 1 > tail.txt
done
cat head.txt body.txt tail.txt > tmp1.trees
mv head.txt ${astral_input_directory}/beast2/context/head.txt
rm body.txt
mv tail.txt ${astral_input_directory}/beast2/context/tail.txt

# Remove comments with logcombiner.py.
resources/logcombiner.py --remove-comments -b 0 -n 0 tmp1.trees tmp2.trees
rm tmp1.trees
cat tmp2.trees | grep "tree " | cut -d " " -f 4 > ${astral_input_directory}/beast2/mcc/combined.tre
rm tmp2.trees

# Make a list of sample tree file names.
ls ${astral_input_directory}/beast2/sample/*.trees > ${astral_input_directory}/beast2/sample/combined_trees.txt

# Prepare a file with an extra tree (the MCC tree resulting from the concatenated BEAST analysis).
resources/logcombiner.py --remove-comments -b 0 -n 0 ../analysis/nuclear_concatenated/beast2/fossil_constrained/combined/76g_nucl_conc_fossils.combined.tre | grep "tree " | cut -d "=" -f 2 > ${astral_input_directory}/beast2/extra/76g_nucl_conc_fossils.combined.tre