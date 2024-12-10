#!/bin/bash -e

##qsub run chicago (2020 updated R script)

echo input files: ${SAMPLE_FOLDERS}_1_100M.SubSam ${SAMPLE_FOLDERS}_2_100M.SubSam

echo TimePoint: ${SAMPLE_NAME} ${TIMEPOINT}

echo Output: ${OUTDIR}

echo feature file: ${FEATFILE}
cat ${FEATFILE}

echo design dir: ${DESIGNDIR}

cd ${OUTDIR}

module load R/4.2.1

Rscript ${OUTDIR}/runChicago_2020_EC.R --design-dir ${DESIGNDIR} --en-feat-list ${FEATFILE} --rda --export-format interBed,washU_text,washU_track ${SAMPLE_FOLDERS}_1_100M.SubSam/${SAMPLE_NAME}_1_100M.SubSam.chinput,${SAMPLE_FOLDERS}_2_100M.SubSam/${SAMPLE_NAME}_2_100M.SubSam.chinput PCHi-C_${TIMEPOINT}_SubSam_chicago_3kb_pooled
