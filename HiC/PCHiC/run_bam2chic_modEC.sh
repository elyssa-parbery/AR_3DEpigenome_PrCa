#!/bin/bash -e

##qsub run bam2chic
echo folder: ${SAMPLE_INPUT}
echo Sample: ${SAMPLE_NAME}
echo Output: ${OUTDIR}
echo baitmap: ${BAITMAP}
echo design dir: ${DESIGNDIR}

module load bedtools/2.28.0

cd ${OUTDIR}

${OUTDIR}/bam2chicago_modEC.sh ${SAMPLE_INPUT} ${BAITMAP} ${DESIGNDIR}/hg38_GATC_GANTC_fragments_3kb_final.rmap ${SAMPLE_NAME} [nodelete]




