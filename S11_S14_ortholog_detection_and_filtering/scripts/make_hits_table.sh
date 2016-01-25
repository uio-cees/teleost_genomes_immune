echo exon_id | tr '\n' '\t'
cat ../analysis/alignments/nuclear/01/ENSDARE00000004235_nucl.fasta | grep ">" | sed 's/>//g' | sed 's/\[.*//g' | tr '\n' '\t'
echo
for i in `tail -n +2 ../../info/selected_nuclear_exons.txt | cut -f 1`
do 
  echo -n -e "${i}\t"; cat ../analysis/alignments/nuclear/01/${i}_nucl.fasta | grep ">" | sed 's/.*\[\&sseqid=//g' | sed 's/,bitscore.*//g' | tr '\n' '\t'
  echo
done