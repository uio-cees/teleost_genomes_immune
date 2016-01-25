#!/bin/bash

echo USAGE: directory support outdirectory filename
test $# == 4 || exit 1

dir=$1
support=$2
out=$3
f=$4 

mkdir -p $out

pd=`pwd`
TMPDIR=$pd
all=`mktemp temp.XXXXX`
ls $dir/*/$f>$all.order
cat $all.order|xargs cat > $all.full

python $BIN_HOME/remove_edges_from_tree.py $all.full $support $all -strip-both

echo "+Group = \"GRAD\"
+Project = \"COMPUTATIONAL_BIOLOGY\"
+ProjectDescription = \"Binning\"

Universe = vanilla

Requirements = Arch == \"X86_64\" 

executable = $BIN_HOME/runcompat.sh

Log = logs/compatibility-$dir-$support.log

getEnv=True 
">condor.compat.$support

for x in `ls $dir`; do
    echo "
 Arguments = $dir $f $x $all $support
 Error=/dev/null
 Output=$out/$x.$support
 Queue">>condor.compat.$support
done
