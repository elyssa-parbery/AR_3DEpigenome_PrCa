##calculate number of "on target" reads per sample

cd ${OUTDIR}

echo INPUT : ${SAMPLE_INPUT}

echo SAMPLE NAME : ${SAMPLE_NAME}

echo OUT DIR : ${OUTDIR}

###Step 1
module load samtools/1.12

CMD="samtools view ${SAMPLE_INPUT} | wc -l"
echo $CMD && eval $CMD

CMD="samtools view ${SAMPLE_INPUT} | wc -l > ${OUTDIR}/${SAMPLE_NAME}.hicup.bam.length.txt"
eval $CMD


##Step 3
module load bedtools/2.28.0

CMD="bedtools intersect -u -bed -a ${SAMPLE_INPUT} -b ${BAITMAP} > ${OUTDIR}/${SAMPLE_NAME}.on-targeted.bed"
echo $CMD && eval $CMD

##Step 4
CMD="cut -f4 ${OUTDIR}/${SAMPLE_NAME}.on-targeted.bed | cut -f1 -d"/" | sort | uniq | wc -l"
echo $CMD && eval $CMD

CMD="cut -f4 ${OUTDIR}/${SAMPLE_NAME}.on-targeted.bed | cut -f1 -d"/" | sort | uniq | wc -l > ${OUTDIR}/${SAMPLE_NAME}.on-target.txt"
eval $CMD

echo DONE ${SAMPLE_NAME}
