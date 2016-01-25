# Michael Matschiner, 2015-06-23

alignment_dir=$1

for i in ${alignment_dir}/*
do
	python3 resources/run_raxml.py -b auto -x Astmex -o Danrer -q ${i}/concatenated.parts --bootstrap-file ${i}/concatenated.bs.trees ${i}/concatenated.phy ${i}/concatenated.tre
done