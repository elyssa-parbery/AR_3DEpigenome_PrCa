tmux new -s MakePeakMatrix_SubSam
#log in 04

qsub -I -P tr07 -q normal -l ncpus=28 -l mem=128GB -l jobfs=100GB -l storage=scratch/tr07+gdata/tr07 -l walltime=05:00:00 -l wd

cd /scratch/tr07/ec2963/PhD_Proj1_AR_Sig/PCHiC_Level_3/chicago

NAMESDIR=/scratch/tr07/PCHiC_Level_2/chicago
OUTDIR=/scratch/tr07/ec2963/PhD_Proj1_AR_Sig/PCHiC_Level_3/BothDHT_SubSam

module load R/4.2.1 intel-compiler/2021.6.0

#update names file to include samples from both time courses

#Rscript makePeakMatrix.R \
#your_full_path_to/names_file.txt \
#hacat_myla_DS_peakMatrix

Rscript makePeakMatrix_EC.R --twopass --notrans --rda --var cd --scorecol score ${NAMESDIR}/names-file_bothDHT_SubSam --cutoff 7 ${OUTDIR}/peakMatrix_LNCaP_Both_DHT_SubSam_min7_3Kb.txt

Rscript makePeakMatrix_EC.R --twopass --notrans --rda --var cd --scorecol score ${NAMESDIR}/names-file_bothDHT_SubSam --cutoff 5 ${OUTDIR}/peakMatrix_LNCaP_Both_DHT_SubSam_min5_3Kb.txt

Rscript makePeakMatrix_EC.R --twopass --notrans --rda --var cd --scorecol score ${NAMESDIR}/names-file_bothDHT_SubSam --cutoff 3 ${OUTDIR}/peakMatrix_LNCaP_Both_DHT_SubSam_min3_3Kb.txt


tmux capture-pane -pS - > /scratch/tr07/ec2963/PhD_Proj1_AR_Sig/PCHiC_Level_3/BothDHT_SubSam/tmux-buffer_makePeakMatrix_bothDHT_SubSam.log
