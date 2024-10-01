module load bedtools/2.28.0 
module load java/jdk-8.40 
CHROMHMM='/g/data/tr07/ec2963/Compute_Assets/Software/ChromHMM'
cd /scratch/tr07/ec2963/PhD_Proj1_AR_Sig/ChromHMM
PWD=$(pwd)
# for TIME in 5mins 10mins 20mins 30mins
# do
#       mkdir ${PWD}/Binarized_data_RNA_${TIME}
#       java -d64 -Xmx16G -jar $CHROMHMM/ChromHMM.jar BinarizeBed -center $CHROMHMM/CHROMSIZES/hg38.txt ${PWD}/data ${PWD}/
BinarizeBed_RNA_${TIME}.config ${PWD}/Binarized_data_RNA_${TIME}
# done
for TIME in 5mins
do
        mkdir ${PWD}/Binarized_data_RNA_${TIME}
        java -d64 -Xmx16G -jar $CHROMHMM/ChromHMM.jar BinarizeBed -center $CHROMHMM/CHROMSIZES/hg38.txt ${PWD}/data ${PWD}/
BinarizeBed_RNA_${TIME}.config ${PWD}/Binarized_data_RNA_${TIME}
done