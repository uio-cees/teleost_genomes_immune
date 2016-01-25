#!/bin/bash -x
  
##SETUP INPUT VARIABLES
i=$1

##SETUP ENVIRONMENT AND COPY FILES NEEDED
module load blast+
module load mafft

##################################### UTG SEARCH #####################################

##ALIGN ALPHA3 REGIONS TO UNITIGS
tblastn \
	-query ZEA.fasta \
	-db ${i}.utg.fasta \
	-out ${i}_ZAE_hits \
	-evalue 1e-15 \
 	-outfmt '6 qseqid sseqid sstart send sseq evalue bitscore length' \
 	-max_target_seqs 30 \
	-num_threads 24 \

##GET ONLY UTGs THATH HIT WITH AT LEAST 75 AMINO ACIDS
cat ${i}_ZAE_hits |\
	awk '$8>70 {print $0}' |\
	cut -f2,5 \
	> ${i}_Good_UTGs_seq

##CONVERT BLASTOUTPUT TO FASTA FILE
cp ~/SCRIPTS/table_to_fasta.rb .
ruby table_to_fasta.rb \
	-f ${i}_Good_UTGs_seq \
	-o ${i}_Good_UTG_hits.fasta \
	-id 1 \
	-seq 2 \

##MAKE ALIGNMENT OF THE EXTRACTED UTG HIT SEQUENCES
mafft \
	--auto \
	${i}_Good_UTG_hits.fasta \
	> ${i}_Good_UTG_hits.aln

##MAKE CONSENSUS OF ALIGNED ALPHA3 SEQUENCES
cp ~/SCRIPTS/produce_majority_rule_consensus.rb .
ruby produce_majority_rule_consensus.rb \
	-f ${i}_Good_UTG_hits.aln \
	-o ${i}_Consensus.fasta