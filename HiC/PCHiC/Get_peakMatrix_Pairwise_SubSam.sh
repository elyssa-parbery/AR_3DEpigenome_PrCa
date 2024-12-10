
CONFIGDIR=/scratch/tr07/PCHiC_Level_2_LNCaP_DHT_SubSam/config
OUTDIR=/scratch/tr07/ec2963/PhD_Proj1_AR_Sig/PCHiC_Level_3/BothDHT_SubSam/PCHi-C_3kb_Res

module load R/4.2.1 intel-compiler/2021.6.0

#update names file to include samples from both time courses

#Rscript makePeakMatrix.R \
#your_full_path_to/names_file.txt \
#hacat_myla_DS_peakMatrix
for TIMEPOINT in 05h 2h 4h 16h 05mins 10mins 20mins 30mins
do
	mkdir ${OUTDIR}/Pairwise_${TIMEPOINT}
	Rscript ${CONFIGDIR}/makePeakMatrix_EC.R --twopass --notrans --rda --var cd --scorecol score ${CONFIGDIR}/names-file_3kb_${TIMEPOINT} --cutoff 7 ${OUTDIR}//Pairwise_${TIMEPOINT}/peakMatrix_LNCaP_${TIMEPOINT}_DHT_SubSam_min7_3Kb.txt
done

