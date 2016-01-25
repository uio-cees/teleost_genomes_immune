# Michael Matschiner, 2015-07-08

# Get the command line arguments.
alignment_dir=$1
raxml_output_dir=$2

# Make the output directory if it does not exist yet.
mkdir -p ${raxml_output_dir}

# Run RAxML for each alignment.
for i in ${alignment_dir}/*
do
	marker_id_w_ext=`basename ${i}`
	marker_id=${marker_id_w_ext%.nex}
	nchar=`cat $i | grep nchar | cut -d "=" -f 3 | sed 's/;//g'`
	echo "DNA, cp1 = 1-${nchar}\2" > tmp.parts
	echo "DNA, cp2 = 2-${nchar}\2" >> tmp.parts
	python3 resources/run_raxml.py -b 100 -x Astmex -o Danrer -q tmp.parts --bootstrap-file ${raxml_output_dir}/${marker_id}.bs.trees ${i} ${raxml_output_dir}/${marker_id}.tre
done