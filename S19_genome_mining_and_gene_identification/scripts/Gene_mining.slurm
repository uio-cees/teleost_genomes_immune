#!/bin/sh

#SBATCH --job-name=GM_2015
#SBATCH --account=XXX
#SBATCH --time=24:0:0
#SBATCH --mem-per-cpu=30000M --partition=hugemem
#SBATCH --cpus-per-task=1
source /cluster/bin/jobsetup

##SETUP INPUT VARIABLES
FISH=$1

##SETUP ENVIRONMENT AND COPY FILES NEEDED
module load blast+/2.2.26

##MAKE FOLDERS
cd DATABASES
mkdir BLASTS_2015
cd BLASTS_2015
mkdir UTG
mkdir SIN
echo "1) Files copied and formatted" \
	> ~/Progressreports_2015/Progress_$FISH 

###################################### UNITIGS SEARCH ######################################
##ALIGN GENELIST TO UNITIGS
cd UTG
tblastn \
	-query ~/Complete_genelist.fas \
	-db $FISH.utg.fasta \
	-out utg-blastout \
	-evalue 1e-10 \
 	-outfmt 6 \
	-max_target_seqs 10 \
	-num_threads 24
echo "2) Unitig BLAST finished" \
	>> ~/Progressreports_2015/Progress_$FISH

##PARSE BEST UNITIG HIT FOR EACH GENE TO NEW FILE
cut utg-blastout \
	-f1,2 | \
	sort -u \
	> All_uniqe_utg_hits_sorted
cat All_uniqe_utg_hits_sorted | \
	awk -F '_' '{print$1}' | \
	sort -u \
	> utg_hits
echo "3) Found initial hits (e-value < 1e-10) in unitigs for the following genes:" \
	>> ~/Progressreports_2015/Progress_$FISH
cat utg_hits \
	>> ~/Progressreports_2015/Progress_$FISH
mkdir UTGHITS
for f in `cat utg_hits`; do \
	grep $f All_uniqe_utg_hits_sorted \
	>> UTGHITS/$f; done
echo "4) Gene-containing utgs parsed to files" \
	>> ~/Progressreports_2015/Progress_$FISH

##EXTRACT READS IN FASTA FORMAT FOR EACH GENE
for f in UTGHITS/*; do \
	cut -f2 $f | \
	sort -u \
	> $f"_utg_reads_only" \
	; done
for f in UTGHITS/*utg_reads_only; do \
	fastagrep \
	-t \
	-f $f ../../$FISH.utg.fasta \
	> ${f%_only}"_fasta" \
	; done
cd UTGHITS
mkdir FASTA
mv *_fasta FASTA
echo "5) Utg-sequences extracted" \
	>> ~/Progressreports_2015/Progress_$FISH

##GET EACH UNITIG TO SEPARATE FILES BEFORE ORF PREDICTION
for f in FASTA/*; do \
	split -l 2 $f \
	${f%_utg_reads_fasta}"_splitted_utg_reads_" \
	; done

##PREDICT ORF FOR ALL GENES IN EACH UNITIG
for f in FASTA/*_splitted_utg_reads_*; do \
	genscan \
	HumanIso.smat $f > \
	$f"_" \
	; done
cd FASTA
mkdir ORFs
mv *_ ORFs

##EXTRACT UNITIG-ORFs SEQUENCE
for f in ORFs/*; do \
	extract_fasta $f \
	>> ${f%_splitted_utg_reads_???}"_merged_ORFs" \
	; done
cd ORFs
mkdir ORFs_FASTA
mv *_ORFs ORFs_FASTA
echo "6) Utg-ORFs predicted and parsed" \
	>> ~/Progressreports_2015/Progress_$FISH

##ALIGN SEQUENCE TO UNIPROT DATABASE
for f in ORFs_FASTA/*; do \
	blastp \
	-query $f \
	-db uniprot_complete_nospace.fasta \
	-out ${f%_merged_ORFs}"_reciprocal_utg_hits" \
	-evalue 1 \
	-outfmt 6 \
	-max_target_seqs 3 \
	-num_threads 24 \
	; done
cd ORFs_FASTA
mkdir RECIPROCAL_HITS
mv *_hits RECIPROCAL_HITS
echo "7) Predicted utg-ORFs aligned to UniProt database" \
	>> ~/Progressreports_2015/Progress_$FISH

##DETERMINE CONFIRMED HITS AND GET ANNOTATION
for f in RECIPROCAL_HITS/*; do \
	awk '$11 < 1e-10 {print $1,"\t",$11,"\t",$2}' \
	> ${f%_reciprocal_utg_hits}"_annotation" $f && \
	echo ${f%_reciprocal_utg_hits} | \
	awk -F '/' '{print$2}' \
	>> ${f%_reciprocal_utg_hits}"_annotation"\
	; done
cd RECIPROCAL_HITS
mkdir ANNOTATED_RECIPROCAL_HITS
mv *_annotation ANNOTATED_RECIPROCAL_HITS

##GET ANNOTATION TO ONE LINE
for f in ANNOTATED_RECIPROCAL_HITS/*; do \
	VAR=`tail -n1 $f` ;
	awk -v var="$VAR" '{print var"\t"$0}' $f \
	> $f"_oneline" \
	; done
cd ANNOTATED_RECIPROCAL_HITS
mkdir ONELINERS
mv *_oneline ONELINERS

##REPORT CONFIRMED HITS AND CONTINUE TO SEARCH THE SINGLETONS
for f in ONELINERS/*; do \
	cat $f | \
	awk -F "|" '{$2="";print $0}' \
	> $f"_beautified" \
	; done
for f in ONELINERS/*_beautified; do \
	head $f | \
	awk '/=/{print$1}' | \
	sort -u \
	>> List_of_genes_found_in_utgs \
	; done
echo "8) Genes present in UTGs, based on reciprocal hits (e-value < 1e-10):" \
	>> ~/Progressreports_2015/Progress_$FISH
for f in ONELINERS/*_beautified; do \
	cat $f | \
	awk '/=/{print$0}' \
	>> ~/Progressreports_2015/Progress_$FISH \
	; done

##################################### SINGLETON SEARCH #####################################

##MAKE A LIST OF GENES FOUND SO FAR
mv List_of_genes_found_in_utgs ${dataDir}/$FISH/CA/9-terminator/DATABASES/BLASTS_2015/SIN/
cd ${dataDir}/$FISH/CA/9-terminator/DATABASES/BLASTS_2015/SIN/

##SELECT ONLY THE GENES THAT WERE NOT FOUND IN SCAFFOLDS OR UNITIGS
grep \
	-f List_of_genes_found_in_utgs \
	-x \
	-v \
	~/Gene_names_only \
	> Still_no_good_hits
echo "9) Looking for the following genes in the singelton reads:" \
	>> ~/Progressreports_2015/Progress_$FISH
cat Still_no_good_hits \
	>> ~/Progressreports_2015/Progress_$FISH

##EKSTRACT FASTA SEQUENCE FOR THE SINGLETON BLAST SEARCH
fastagrep \
	-t \
	-f Still_no_good_hits \
	~/Complete_genelist.fas \
	> Still_no_good_hits_fas

##ALIGN REDUCED GENELIST TO SINGLETON READS
tblastn \
	-query Still_no_good_hits_fas \
	-db ../../$FISH.singleton.fasta \
	-out sin-blastout \
	-evalue 1e-1 \
	-outfmt 6 \
	-max_target_seqs 10 \
	-num_threads 24
echo "10) Singleton BLAST finished" \
	>> ~/Progressreports_2015/Progress_$FISH

##PARSE BEST HIT FOR EACH GENE TO NEW FILE
cut sin-blastout \
	-f1,2 | \
	sort -u \
	> All_uniqe_sin_hits_sorted
cat All_uniqe_sin_hits_sorted | \
	awk -F '_' '{print$1}' | \
	sort -u \
	> sin_hits
echo "11) Found initial hits in singletons (e-value < 1e-1) for the following genes:" \
	>> ~/Progressreports_2015/Progress_$FISH
cat sin_hits \
	>> ~/Progressreports_2015/Progress_$FISH
mkdir SINHITS
for f in `cat sin_hits`; do \
	grep $f All_uniqe_sin_hits_sorted \
	>> SINHITS/$f; done
echo "12) Gene-containing singletons parsed to files" \
	>> ~/Progressreports_2015/Progress_$FISH

##EXTRACT READS IN FASTA FORMAT FOR EACH GENE
for f in SINHITS/*; do \
	cut -f2 $f | \
	sort -u \
	> $f"_sin_reads_only" \
	; done
for f in SINHITS/*sin_reads_only; do \
	fastagrep \
	-t \
	-f $f ../../$FISH.singleton.fasta \
	> ${f%_only}"_fasta" \
	; done
cd SINHITS
mkdir FASTA
mv *_fasta FASTA
echo "13) Singleton-sequences extracted" \
	>> ~/Progressreports_2015/Progress_$FISH

##GET EACH SINGLETON TO SEPARATE FILES BEFORE ORF PREDICTION
for f in FASTA/*; do \
	split -l 2 $f \
	${f%_sin_reads_fasta}"_splitted_singletons_" \
	; done

##PREDICT ORF FOR ALL GENES IN EACH SINGLETON
for f in FASTA/*_splitted_singletons_*; do \
	genscan \
	/projects/454data/bin/GenScan/HumanIso.smat $f > \
	$f"_" \
	; done
cd FASTA
mkdir ORFs
mv *_ ORFs

##EXTRACT SEQUENCE FOR SINGLETON ORFS
for f in ORFs/*; do \
	extract_fasta $f \
	>> ${f%_splitted_singletons_???}"_merged_ORFs" \
	; done
cd ORFs
mkdir ORFs_FASTA
mv *_ORFs ORFs_FASTA
echo "14) Singleton-ORFs predicted and parsed" \
	>> ~/Progressreports_2015/Progress_$FISH

##ALIGN SINGLETONS TO UNIPROT DATABASE
for f in ORFs_FASTA/*; do \
	blastp \
	-query $f \
	-db ~/uniprot_complete_nospace.fasta \
	-out ${f%_merged_ORFs}"_reciprocal_sin_hits" \
	-evalue 1 \
	-outfmt 6 \
	-max_target_seqs 20 \
	-num_threads 24 \
	; done
cd ORFs_FASTA
mkdir RECIPROCAL_HITS
mv *_hits RECIPROCAL_HITS
echo "15) Predicted Singleton-ORFs aligned to UniProt database" \
	>> ~/Progressreports_2015/Progress_$FISH

##DETERMINE CONFIRMED HITS AND GET ANNOTATION
for f in RECIPROCAL_HITS/*; do \
	awk '$11 < 1e-1 {print $1,"\t",$11,"\t",$2}' \
	> ${f%_reciprocal_sin_hits}"_annotation" $f && \
	echo ${f%_reciprocal_sin_hits} | \
	awk -F '/' '{print$2}' \
	>> ${f%_reciprocal_sin_hits}"_annotation"\
	; done
cd RECIPROCAL_HITS
mkdir ANNOTATED_RECIPROCAL_HITS
mv *_annotation ANNOTATED_RECIPROCAL_HITS

##GET ANNOTATION TO ONE LINE
for f in ANNOTATED_RECIPROCAL_HITS/*; do \
	VAR=`tail -n1 $f` ;
	awk -v var="$VAR" '{print var"\t"$0}' $f \
	> $f"_oneline" \
	; done
cd ANNOTATED_RECIPROCAL_HITS
mkdir ONELINERS
mv *_oneline ONELINERS

##REPORT CONFIRMED HITS 
for f in ONELINERS/*; do \
	cat $f | \
	awk -F "|" '{$2="";print $0}' \
	> $f"_beautified" \
	; done
for f in ONELINERS/*_beautified; do \
	head $f | \
	awk '/=/{print$1}' \
	>> List_of_genes_found_in_singletons \
	; done
echo "16) Genes present in singletons, based on reciprocal hits (e-value < 1e-1):" \
	>> ~/Progressreports_2015/Progress_$FISH
for f in ONELINERS/*_beautified; do \
	head $f | \
	awk '/=/{print$0}' \
	>> ~/Progressreports_2015/Progress_$FISH \
	; done

##REPORT WHICH GENES ARE ABSENT 
mv List_of_genes_found_in_singletons ${dataDir}/$FISH/CA/9-terminator/DATABASES/BLASTS_2015/SIN/
cd ${dataDir}/$FISH/CA/9-terminator/DATABASES/BLASTS_2015/SIN/
cat List_of_genes_found_in_singletons \
	List_of_genes_found_in_utgs \
	> All_genes_present_in_$FISH
grep \
	-f All_genes_present_in_$FISH \
	-x \
	-v \
	~/Gene_names_only \
	> Genes_absent
echo "17) These genes are NOT present in the genome:" \
	>> ~/Progressreports_2015/Progress_$FISH
cat Genes_absent \
	>> ~/Progressreports_2015/Progress_$FISH

