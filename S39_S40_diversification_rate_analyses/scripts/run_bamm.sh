# Michael Matschiner, 2015-09-09.

home=`pwd`
analysis_dir=${1}
control_file_name="control.txt"

for replicate_dir in ${analysis_dir}/r*
do
	cd ${replicate_dir}
	./bamm -c ${control_file_name}
	cd ${home}
done