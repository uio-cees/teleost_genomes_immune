# Michael Matschiner, 2015-04-07.

# Get the command line arguments.
replicate_dir=$1
combined_dir=$2

# Make the directory for the combined BEAST analyses.
mkdir -p ${combined_dir}

# Make a file with a list of log and trees files to be combined.
ls ${replicate_dir}/*/*.trees > ${combined_dir}/tree_files.txt
ls ${replicate_dir}/*/*.log > ${combined_dir}/log_files.txt

# Determine the name of the combined log and trees file.
first_replicate_file_name_w_path=`cat ${combined_dir}/tree_files.txt | head -n 1`
first_replicate_file_name=`basename ${first_replicate_file_name_w_path}`
combined_trees_file_name=${first_replicate_file_name%.trees}.combined.trees
mcc_tree_file_name=${first_replicate_file_name%.trees}.combined.tre
combined_log_file_name=${first_replicate_file_name%.trees}.combined.log

# Use logcombiner.py to combine the trees files.
python3 resources/logcombiner.py -b 10 ${combined_dir}/tree_files.txt ${combined_dir}/${combined_trees_file_name}
python3 resources/logcombiner.py -b 10 ${combined_dir}/log_files.txt ${combined_dir}/${combined_log_file_name}

# Run TreeAnnotator to produce a MCC tree from the tree sample.
java -Xms256m -Xmx2g -Djava.library.path="resources/lib" -cp "resources/lib/launcher.jar" beast.app.treeannotator.TreeAnnotatorLauncher -heights mean ${combined_dir}/${combined_trees_file_name} ${combined_dir}/${mcc_tree_file_name}