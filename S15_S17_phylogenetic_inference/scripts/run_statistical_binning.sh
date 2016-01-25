# Michael Matschiner, 2015-04-15

# Get command line arguments.
tree_dir_in=$1
tree_dir_out=$2
threshold=$3

# Clean up from previous runs if there were any.
rm $tree_dir_out/*

# Make commands.
bash ./resources/statistical_binning/makecommands.compatibility.sh ${tree_dir_in} ${threshold} ${tree_dir_out} ML.tre
# Run commands.
bash commands.compat.${threshold}

# Clean up.
rm temp*
rm commands.compat.${threshold}

# Build the bin definitions.
cd $tree_dir_out
ls| grep -v ge|sed -e "s/.${threshold}$//g" > genes
python ../../../../scripts/resources/statistical_binning/cluster_genetrees.py genes ${threshold}

# Clean up.
mkdir ../remove
mv * ../remove
mv ../remove/bin.*.txt .
rm -r ../remove