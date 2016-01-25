#!/bin/bash

test $# == 5 || echo USAGE: directory filename gene versus support

BIN_HOME=${0%\/runcompat.sh}
dir=$1
f=$2
x=$3
all=$4
supp=$5

if [ "$supp" == "" ]; then 
  temp=$dir/$x/$f
else
  temp=`mktemp /tmp/temp.XXXXXX`
  python $BIN_HOME/remove_edges_from_tree.py $dir/$x/$f $supp $temp -strip-both
fi

dir_edit=`echo "$dir" | tr '/' ':' | sed -e 's/:/\\\\\//g'`
$BIN_HOME/compareTrees.compatibility $temp $all |paste -d " " $all.order - | sed -e "s/${dir_edit}\//${x} /g" -e "s/\/[^ ]* / /g"
