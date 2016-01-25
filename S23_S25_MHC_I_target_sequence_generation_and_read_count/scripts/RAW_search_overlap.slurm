#!/bin/sh

#SBATCH --job-name=RS_MHCI
#SBATCH --account=nn9244k
#SBATCH --time=48:0:0
#SBATCH --mem-per-cpu=4000M 
#SBATCH --cpus-per-task=8
source /cluster/bin/jobsetup

  
##SETUP INPUT VARIABLES
FISH=$1


###################################   PREFORM BLAST  #####################################

tblastn \
	-query ~/MCMC_METODE/QUERIES/All_queries.fas \
	-db /work/users/martinma/Raw_db/${FISH}/${FISH}_raw_raw.fasta \
	-out raw-raw-blastout \
	-evalue 1 \
	-outfmt 6 \
	-max_target_seqs 100000 \
	-num_threads 24 

##################################   INDIVIDUAL SPECIES  #################################

for QUERY in $(cat Query_names) 
do
##EXTRACT UNIQUE RAW HITS FOR EACH QUERY 
cat raw-raw-blastout | grep ${QUERY} >> ${QUERY}_all
echo "TOTAL NUMBER OF UNIQUE RAW HITS WITH ${QUERY}:" \
		>> RAW_search_overlap_${FISH}.full_output
cat ${QUERY}_all |\
	cut -f2 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.full_output
	 
##GET UNIQUE RAW HITS FOR ${QUERY} WITH AT LEAST 75bp overlap
cat ${QUERY}_all | awk '{if ($9 > $10 && (($9-$10+1) > 75)) print $0; else if (($10-$9+1) > 75) print $0;}' >> ${QUERY}_overlap
echo "${QUERY} WITH AT LEAST 75bp OVERLAP:" \
		>> RAW_search_overlap_${FISH}.full_output
cat ${QUERY}_overlap |\
	cut -f2 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.full_output

cat ${QUERY}_overlap |\
	awk '$11 < 1e-5 {print $2}' \
	>> ${QUERY}_e-5
echo "${QUERY} HITS WITH 75bp OVERLAP AND E-VALUE BELOW 1e-5" \
		>> RAW_search_overlap_${FISH}.full_output
cat ${QUERY}_e-5 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.full_output

cat ${QUERY}_overlap |\
	awk '$11 < 1e-10 {print $2}' \
	>> ${QUERY}_e-10
echo "${QUERY} HITS WITH 75bp OVERLAP AND E-VALUE BELOW 1e-10" \
		>> RAW_search_overlap_${FISH}.full_output
cat ${QUERY}_e-10 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.full_output

cat ${QUERY}_overlap |\
	awk '$11 < 1e-15 {print $2}' \
	>> ${QUERY}_e-15
echo "${QUERY} HITS WITH 75bp OVERLAP AND E-VALUE BELOW 1e-15" \
		>> RAW_search_overlap_${FISH}.full_output
cat ${QUERY}_e-15 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.full_output

cat ${QUERY}_overlap |\
	awk '$11 < 1e-20 {print $2}' \
	>> ${QUERY}_e-20
echo "${QUERY} HITS WITH 75bp OVERLAP AND E-VALUE BELOW 1e-20" \
		>> RAW_search_overlap_${FISH}.full_output
cat ${QUERY}_e-20 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.full_output

done

###################################   COMBINE SPECIES   #################################

for TARGET in $(cat ~/MCMC_METODE/QUERIES/Targets) 
do

##SUM OF UNIQUE READS FOR ${TARGET} FROM OWN SPECIES + ALL ENSEMBLE SPECIES (-COD)
cat ${TARGET}_Danrer_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_Gasacu_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_Orylat_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_Orenil_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_Takrub_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_Tetnig_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_Poefor_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_Xipmac_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_Astmex_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all
cat ${TARGET}_${FISH}_all >> Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all

cat Own_and_ENSEMBLE_nocod_${TARGET}_uncut_all | awk '{if ($9 > $10 && (($9-$10+1) > 75)) print $0; else if (($10-$9+1) > 75) print $0;}' >> Own_and_ENSEMBLE_nocod_${TARGET}_all

cat Own_and_ENSEMBLE_nocod_${TARGET}_all |\
	awk '$11 < 1e-5 {print $2}' \
	>> Own_and_ENSEMBLE_nocod_${TARGET}-5
echo "${TARGET} HITS FOR OWN + ENSEMBLE (-COD) WITH E-VALUE BELOW 1e-5:" \
		>> RAW_search_overlap_${FISH}.out1
cat Own_and_ENSEMBLE_nocod_${TARGET}-5 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.out1

cat Own_and_ENSEMBLE_nocod_${TARGET}_all |\
	awk '$11 < 1e-10 {print $2}' \
	>> Own_and_ENSEMBLE_nocod_${TARGET}-10
echo "${TARGET} HITS FOR OWN + ENSEMBLE (-COD) WITH E-VALUE BELOW 1e-10:" \
		>> RAW_search_overlap_${FISH}.out1
cat Own_and_ENSEMBLE_nocod_${TARGET}-10 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.out1

cat Own_and_ENSEMBLE_nocod_${TARGET}_all |\
	awk '$11 < 1e-15 {print $2}' \
	>> Own_and_ENSEMBLE_nocod_${TARGET}-15
echo "${TARGET} HITS FOR OWN + ENSEMBLE (-COD) WITH E-VALUE BELOW 1e-15:" \
		>> RAW_search_overlap_${FISH}.out1
cat Own_and_ENSEMBLE_nocod_${TARGET}-15 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.out1
	 
cat Own_and_ENSEMBLE_nocod_${TARGET}_all |\
	awk '$11 < 1e-20 {print $2}' \
	>> Own_and_ENSEMBLE_nocod_${TARGET}-20
echo "${TARGET} HITS FOR OWN + ENSEMBLE (-COD) WITH E-VALUE BELOW 1e-20:" \
		>> RAW_search_overlap_${FISH}.out1
cat Own_and_ENSEMBLE_nocod_${TARGET}-20 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.out1


##SUM OF UNIQUE READS FOR ${TARGET} FROM USED SPECIES ONLY
cat ${TARGET}_fish_*_all >> All_used_species_${TARGET}_uncut_all
cat All_used_species_${TARGET}_uncut_all | awk '{if ($9 > $10 && (($9-$10+1) > 75)) print $0; else if (($10-$9+1) > 75) print $0;}' >> All_used_species_${TARGET}_all

cat All_used_species_${TARGET}_all |\
	awk '$11 < 1e-5 {print $2}' \
	>> All_used_species_${TARGET}-5
echo "${TARGET} HITS FOR ALL USED SPECIES WITH E-VALUE BELOW 1e-5:" \
		>> RAW_search_overlap_${FISH}.out2
cat All_used_species_${TARGET}-5 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.out2

cat All_used_species_${TARGET}_all |\
	awk '$11 < 1e-10 {print $2}' \
	>> All_used_species_${TARGET}-10
echo "${TARGET} HITS FOR ALL USED SPECIES WITH E-VALUE BELOW 1e-10:" \
		>> RAW_search_overlap_${FISH}.out2
cat All_used_species_${TARGET}-10 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.out2

cat All_used_species_${TARGET}_all |\
	awk '$11 < 1e-15 {print $2}' \
	>> All_used_species_${TARGET}-15
echo "${TARGET} HITS FOR ALL USED SPECIES WITH E-VALUE BELOW 1e-15:" \
		>> RAW_search_overlap_${FISH}.out2
cat All_used_species_${TARGET}-15 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.out2
	 
cat All_used_species_${TARGET}_all |\
	awk '$11 < 1e-20 {print $2}' \
	>> All_used_species_${TARGET}-20
echo "${TARGET} HITS FOR ALL USED SPECIES WITH E-VALUE BELOW 1e-20:" \
		>> RAW_search_overlap_${FISH}.out2
cat All_used_species_${TARGET}-20 |\
     sort -u |\
	 wc -l \
	 >> RAW_search_overlap_${FISH}.out2

done	
	 
###############################    CLEAN UP MESS     ####################################	 

##REMOVE EXTRA LINES FOR EASIER DOWNSTRAM HANDELING
sed -n 2~2p RAW_search_overlap_${FISH}.out1 > RAW_count_${FISH}.out1
sed -n 2~2p RAW_search_overlap_${FISH}.out2 > RAW_count_${FISH}.out2

##MOVE FINAL OUTPUT
#mv RAW_search_overlap_${FISH}.out1 ~/MCMC_METODE/OUTPUT/
#mv RAW_search_overlap_${FISH}.out2 ~/MCMC_METODE/OUTPUT/
mv RAW_count_${FISH}.out1 ~/MCMC_METODE/OUTPUT/
mv RAW_count_${FISH}.out2 ~/MCMC_METODE/OUTPUT/

##GATHER OUTPUT
mkdir RAW_overlap_counts
mv *_all RAW_overlap_counts
mv *_overlap RAW_overlap_counts
mv *-5 RAW_overlap_counts
mv *-10 RAW_overlap_counts
mv *-15 RAW_overlap_counts
mv *-20 RAW_overlap_counts


