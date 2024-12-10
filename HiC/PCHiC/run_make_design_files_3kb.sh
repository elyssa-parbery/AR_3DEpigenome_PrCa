##Create design files from baitmap and rmap for use with CHiCAGO

cd /scratch/tr07/PCHiC_Level_2_LNCaP_DHT_SubSam/chicago_3kb

WORKDIR=/scratch/tr07/PCHiC_Level_2_LNCaP_DHT_SubSam/chicago
DESIGN=/scratch/tr07/PCHiC_Level_2_LNCaP_DHT_SubSam/chicago/DESIGN-DIR	

module load python2/2.7.16 

python2 ${WORKDIR}/makeDesignFiles.py --rmapfile=${DESIGN}/hg38_GATC_GANTC_fragments_final.rmap --baitmapfile=${DESIGN}/hg38_baitmap_final_sort.uniq.5col.baitmap --outfilePrefix=hg38_GATC_GANTC_fragments_5kb --minFragLen=25 --maxFragLen=1200 --maxLBrownEst=75000 --binsize=3000 --removeb2b=True --removeAdjacent=True
