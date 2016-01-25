# Michael Matschiner, 2015-07-01.

# Get command line arguments
tree_file_name_w_relative_path=${1}
tree_file_name=`basename ${tree_file_name_w_relative_path}`
replicate_dir=${2}
number_of_reps=${3}

# Remove taxa without trait data from the tree.
rscript drop_ensembl_taxa.r ${tree_file_name_w_relative_path} tmp.tre

# Use add_lepidogalaxias.rb to add Lepidogalaxias to trees, with an age randomly chosen along the root branch leading to Neoteleostei.
for i in `seq 1 ${number_of_reps}`
do
	n=`printf "%3.3s\n" "$i" | tr " " "0"`
	ruby add_lepidogalaxias.rb tmp.tre ${replicate_dir}/r${n}/${tree_file_name%.tre}.added_r${n}.tre
done

# Clean up.
rm tmp.tre