# Michael Matschiner, 2015-07-07

# Get the command line arguments.
input_tree_dir=$1
output_tree_dir=$2
raxml_tree_dir_stem=$3
beast2_tree_dir_stem=$4

# Make the Astral input and output directories if they don't exist yet.
mkdir -p ${input_tree_dir}
mkdir -p ${output_tree_dir}

# Prepare the input trees for the Astral analysis with BEAST trees.
if [ -f ${input_tree_dir}/beast2.trees ]
then
	rm ${input_tree_dir}/beast2.trees
fi
for i in ${beast2_tree_dir_stem}/bin_?/combined/bin_?.tre
do
	rscript drop_single_tip_from_tree.r ${i} Astmex tmp.tre
	cat tmp.tre >> ${input_tree_dir}/beast2.trees
	rm tmp.tre
done
for i in ${beast2_tree_dir_stem}/bin_?/combined/bin_?.100.trees
do
	trees_file_name=`basename $i`
	rscript drop_single_tip_from_trees.r ${i} Astmex ${input_tree_dir}/beast2.${trees_file_name}
done

# Run Astral with BEAST trees.
ls ${input_tree_dir}/beast2.bin_?.100.trees > tmp.txt
java -jar -Xmx8G resources/Astral/astral.4.7.8.jar -i ${input_tree_dir}/beast2.trees -b tmp.txt -o ${output_tree_dir}/beast2.trees
rm tmp.txt
cat ${output_tree_dir}/beast2.trees | tail -n 1 > ${output_tree_dir}/beast2.bins.tre
rm ${output_tree_dir}/beast2.trees

# Prepare the input trees for the Astral analysis with RAxML trees.
cat ${raxml_tree_dir_stem}/bin_?/concatenated.tre > ${input_tree_dir}/raxml.trees
for i in ${raxml_tree_dir_stem}/bin_?/concatenated.bs.trees
do
	bin=`echo $i | rev | cut -d '/' -f 2 | rev`
	cat $i | head -n 100 > ${input_tree_dir}/raxml.${bin}.100.trees
done

# Run Astral with RAxML trees.
ls ${input_tree_dir}/raxml.bin_?.100.trees > tmp.txt
java -jar -Xmx8G resources/Astral/astral.4.7.8.jar -i ${input_tree_dir}/raxml.trees -b tmp.txt -o ${output_tree_dir}/raxml.trees
rm tmp.txt
cat ${output_tree_dir}/raxml.trees | tail -n 1 > ${output_tree_dir}/raxml.bins.tre
rm ${output_tree_dir}/raxml.trees
