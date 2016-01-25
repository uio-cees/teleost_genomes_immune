#!/bin/bash -x

##SETUP INPUT VARIABLES
i=$1

##SETUP ENVIRONMENT AND COPY FILES NEEDED
module load blast+
module load mafft

##################################### UTG SEARCH #####################################

##ALIGN ALPHA3 REGIONS TO UNITIGS
tblastn \
	-query Target_seq_A3.fas \
	-db ${i}.utg.fasta \
	-out consensus-blastout_utg \
	-evalue 1e-5 \
 	-outfmt '6 qseqid sseqid sstart send sseq evalue bitscore length' \
 	-max_target_seqs 100000 \
	-num_threads 24 \

##GET ONLY UTGs THATH HIT WITH AT LEAST 75 AMINO ACIDS
cat consensus-blastout_utg |\
	awk '$8>70 {print $0}' |\
	cut -f2,5 \
	> Good_UTGs_seq

##CONVERT BLASTOUTPUT TO FASTA FILE
cp ~/SCRIPTS/table_to_fasta.rb .
ruby table_to_fasta.rb \
	-f Good_UTGs_seq \
	-o Good_UTG_hits.fasta \
	-id 1 \
	-seq 2 \

##MAKE ALIGNMENT OF THE EXTRACTED UTG HIT SEQUENCES
mafft \
	--auto \
	Good_UTG_hits.fasta \
	> Good_UTG_hits.aln

##MAKE CONSENSUS OF ALIGNED ALPHA3 SEQUENCES
cp ~/SCRIPTS/produce_majority_rule_consensus.rb .
ruby produce_majority_rule_consensus.rb \
	-f Good_UTG_hits.aln \
	-o Consensus_${i}.fasta