#!/bin/bash

echo "USAGE: bindir align_home[full_path]"
test $# -gt 1 || exit 1

BINDIR=$1
GENEDIR="$2"
OUTDIR="supergenes"

EXT=fasta

mkdir -p $OUTDIR

for y in `wc -l $BINDIR/bin*txt|grep -v total|awk '{if ($1==1)print $2}'`; do 
  x=`echo $y|sed -e "s/.*bin/bin/g"`
  g=`cat $y`
  mkdir $OUTDIR/$x
  ln -fs $GENEDIR/$g/raxmlboot.gtrgamma $OUTDIR/$x  # this is for avian simulated
  echo "Done" > $OUTDIR/$x/.done.raxml.gtrgamma.1
  echo "Done" > $OUTDIR/$x/.done.raxml.gtrgamma.200.2
  echo $g > $OUTDIR/$x/$OUTFILE.part
done

for y in `wc -l $BINDIR/bin*txt|grep -v total|awk '{if ($1>1)print $2}'`; do 
  cat $y|xargs -I@ echo $GENEDIR/@/@.$EXT >.t;  
  x=`echo $y|sed -e "s/.*bin/bin/g"`
  mkdir $OUTDIR/$x
  $BIN_HOME/perl/concatenate_alignments.pl -i `pwd`/.t -o `pwd`/$OUTDIR/$x/$OUTFILE.fasta -p `pwd`/$OUTDIR/$x/$OUTFILE.part;
  tail -n1 `pwd`/$OUTDIR/$x/$OUTFILE.part 
  # convert_to_phylip.sh `pwd`/913supergenes/$x/sate.noout.fasta 913supergenes/$x/sate.noout.phylip; 
  echo $x done; 
done
