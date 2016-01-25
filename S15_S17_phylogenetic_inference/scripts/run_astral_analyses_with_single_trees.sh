# Michael Matschiner, 2015-07-07

# Define names and make directories.
astral_input_directory=$1
astral_output_directory=$2
alignment_directory=$3
beast_input_tree_dir_stem=$4
raxml_input_tree_dir_stem=$5
mkdir -p ${astral_input_directory}/beast2/mcc
mkdir -p ${astral_input_directory}/beast2/sample
mkdir -p ${astral_input_directory}/beast2/context
mkdir -p ${astral_input_directory}/raxml/ml
mkdir -p ${astral_input_directory}/raxml/bootstraps

# Run treeannotator for the posterior tree distributions of the BEAST2 analysis done during ortholog detection.
for i in ${alignment_directory}/*.nex
do
	gene_id_w_ext=`basename ${i}`
	gene_id=${gene_id_w_ext%.nex}
	trees_file_directory="${beast_input_tree_dir_stem}/${gene_id}"
	trees_file_name="${gene_id}.trees"
	mcc_tree_file_name=${trees_file_name%.trees}.tre
	sample_tree_file_name=${trees_file_name%.trees}.trees
	if [ ! -f ${astral_input_directory}/beast2/sample/${sample_tree_file_name} ]
	then
		resources/treeannotator/bin/treeannotator -burnin 20 -heights mean ${trees_file_directory}/${trees_file_name} ${astral_input_directory}/beast2/mcc/${mcc_tree_file_name}
		resources/logcombiner.py --remove-comments -b 20 -n 100 ${trees_file_directory}/${trees_file_name} tmp.trees
		cat tmp.trees | grep "tree " | cut -d "=" -f 2 > tmp2.trees
		rscript drop_single_tip_from_newick_trees.r tmp2.trees 1 ${astral_input_directory}/beast2/sample/${sample_tree_file_name}
		rm tmp.trees
		rm tmp2.trees
	fi
done

# Prepare a single file with all BEAST2 MCC trees.
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
cat tmp2.trees | grep "tree " | cut -d " " -f 4 > tmp3.trees
rm tmp2.trees
rscript drop_single_tip_from_newick_trees.r tmp3.trees 1 ${astral_input_directory}/beast2/mcc/combined.tre
rm tmp3.trees

# Prepare a single file with all RAxML ML trees.
cat ${raxml_input_tree_dir_stem}/*.tre > ${astral_input_directory}/raxml/ml/combined.tre

# Copy all RAxML bootstrap trees to the Astral input directory for RAxML trees.
for i in ${raxml_input_tree_dir_stem}/*.bs.trees
do
	cp ${i} ${astral_input_directory}/raxml/bootstraps/
done

# Make lists of sample/bootstrap tree file names.
ls ${astral_input_directory}/beast2/sample/*.trees > ${astral_input_directory}/beast2/sample/combined_trees.txt
ls ${astral_input_directory}/raxml/bootstraps/*.trees > ${astral_input_directory}/raxml/bootstraps/combined_trees.txt

# Make the Astral output directory if it does not exist yet.
mkdir -p ${astral_output_directory}

# Run Astral with BEAST2 trees based on single gene alignments.
java -jar -Xmx8G resources/Astral/astral.4.7.8.jar -i ${astral_input_directory}/beast2/mcc/combined.tre -b ${astral_input_directory}/beast2/sample/combined_trees.txt -o body.txt
cat body.txt | tail -n 1 > body_tail.txt
rm body.txt
echo -n "tree astral = " > tree_line.txt
cat ${astral_input_directory}/beast2/context/head.txt tree_line.txt body_tail.txt ${astral_input_directory}/beast2/context/tail.txt > ${astral_output_directory}/beast2.single.tre
rm body_tail.txt
rm tree_line.txt

# Run Astral with RAxML trees based on single gene alignments.
java -jar -Xmx8G resources/Astral/astral.4.7.8.jar -i ${astral_input_directory}/raxml/ml/combined.tre -b ${astral_input_directory}/raxml/bootstraps/combined_trees.txt -o tmp.txt
cat tmp.txt | tail -n 1 > ${astral_output_directory}/raxml.single.tre
rm tmp.txt
