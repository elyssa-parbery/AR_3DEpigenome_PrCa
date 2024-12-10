## Running the nfcore cutandrun pipeline from a tmux session on Gadi
## Author: Elyssa Campbell
## Date written: 5/09/2023
## User:
## Date run:


# log in to Gadi via ssh or the ARE web app
# cd to working dir on login node

# record the login node number youre in as you'll need to be on the same node to access yout tmux session later
## login node: 

tmux new -s yoursessionname

module purge
module load nextflow
module load singularity

# change path to your working Dir:
WORKDIR=/scratch/tr07/ec2963/CutRun_Level_2_hg38_LNCaP_EarlyDHT
# change path to your samplesheet:
SAMPLESHEET=/scratch/tr07/ec2963/CutRun_Level_2_hg38_LNCaP_EarlyDHT/samplesheet_Gadi_EarlyDHT_SubSam_CTCF.csv
# use my iGenomes references dir or change path to your iGenomes reference dir:
iGENOMES=/scratch/tr07/ec2963/assets/genomes/iGenomes/references
# use my cut and run blacklist or change path to your own copy:
BLACKLIST=/scratch/tr07/ec2963/assets/CutandRun_hg38_suspectlist_ANordinetal2023.bed
# set the BATCH variable to identify your unique batch eg, if you've run batches per target use the name of the target.
BATCH="CTCF"

####
# only pull the first time you run the cutandrun pipeline to test nf core is working for you:
nextflow pull cutandrun/main.nf
####

nextflow run cutandrun/main.nf \
     -profile singularity,nci_gadi \
     --input ${SAMPLESHEET} --outdir ${WORKDIR}/nfcore_results_${BATCH} \
     --genome GRCh38 --igenomes_base ${iGENOMES} --blacklist ${BLACKLIST} \
     --normalisation_mode "CPM" --replicate_threshold 2 --only_alignment -resume

#keep a record of your run name, you'll need this if you have to resume your pipe:
## runName : 

# press control + b then press d to detach from the tmux session 
# you can return to this job on the same login node using the following command:
tmux attach -t yoursessionname


