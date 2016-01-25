# Michael Matschiner, 2015-07-07

analysis_dir=$1
home=`pwd`
for i in ${analysis_dir}/bin_?/replicates/r??
do
	cd ${i}
	sbatch start_accel.sh
	cd ${home}
done