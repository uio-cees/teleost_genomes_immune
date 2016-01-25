# Michael Matschiner, 2015-07-07

analysis_dir=$1
for i in ${analysis_dir}/bin_?
do
	bin=`basename $i`
	mkdir -p ${i}/combined
	ls ${i}/replicates/r??/*.log > ${i}/combined/log_file_names.txt
	ls ${i}/replicates/r??/*.trees > ${i}/combined/tree_file_names.txt
	resources/logcombiner.py -n 500 ${i}/combined/log_file_names.txt ${i}/combined/${bin}.log
	resources/logcombiner.py -n 500 ${i}/combined/tree_file_names.txt ${i}/combined/${bin}.trees
	resources/logcombiner.py -n 50 ${i}/combined/tree_file_names.txt ${i}/combined/${bin}.100.trees
	resources/treeannotator/bin/treeannotator -heights mean ${i}/combined/${bin}.trees ${i}/combined/${bin}.tre
done