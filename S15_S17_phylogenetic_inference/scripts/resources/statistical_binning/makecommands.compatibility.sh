#!/bin/bash

echo USAGE: directory support outdirectory filename
test $# == 4 || exit 1

BIN_HOME=${0%\/makecommands.compatibility.sh}
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

executable=$BIN_HOME/runcompat.sh
echo "#!/bin/bash" >commands.compat.$support
for x in `ls $dir`; do
    echo "
$executable  $dir $f $x $all $support 1>$out/$x.$support">>commands.compat.$support
done
