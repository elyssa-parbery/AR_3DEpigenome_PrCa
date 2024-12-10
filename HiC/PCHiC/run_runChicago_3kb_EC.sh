#!/bin/bash -e

##qsub run chicago (2020 updated R script)

echo folder: ${SAMPLE_FOLDER}
echo Sample: ${SAMPLE_NAME}
echo file: ${FILE_NAME}
echo Output: ${OUTDIR}
echo feature file: ${FEATFILE}
cat ${FEATFILE}

echo design dir: ${DESIGNDIR}

cd ${OUTDIR}

module load R/4.2.1

Rscript ${OUTDIR}/runChicago_2020_EC.R --design-dir ${DESIGNDIR} --en-feat-list ${FEATFILE} --rda --export-format interBed,washU_text,washU_track ${SAMPLE_FOLDER}/${FILE_NAME}.chinput ${SAMPLE_NAME}_SubSam_chicago_3kb
