. /home/eservant/bioprog/promline/bench.config


mkdir $BENCH
cd $BENCH

# get > 30X 
samtools depth $EXBAM -Q 5 | awk -v OFS='\t' '$3>30 {print $1,$2-1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -i stdin > 30X_${EXOME}.bed

# filter 30X and bad GT
bcftools filter -i 'FORMAT/DP>30 & FORMAT/GT!="mis" & FORMAT/GT!="0/0"' $EXVCF -o 30X_${EXOME}.vcf.gz -Oz

#liftover 30X.bed 19 to 38
/home/eservant/bioprog/liftover/liftOver 30X_${EXOME}.bed /home/eservant/bioprog/liftover/chains/hg19ToHg38.over.chain.gz 38_30X_${EXOME}.bed unMapped30X_38to19.bed

#liftovervcf
/home/eservant/bioprog/gatk-4.4.0.0/gatk LiftoverVcf \
    -I 30X_${EXOME}.vcf.gz \
    -O 38_30X_${EXOME}.vcf.gz \
    -CHAIN /home/eservant/bioprog/liftover/chains/hg19ToHg38.over.chain.gz \
    --REJECT rejected_variants30X_38to19.vcf \
    -R ${HGREF} \
    --CREATE_INDEX true \
    -WMC true


# filter promethion VCF with regions in merged AND in gencode CDS
#gencode
bcftools filter $VCFPROM -R /home/eservant/bioprog/Annotations/hg38/gencode-v41.basic.annotation.CDS.lvl1-2-3.sortU.pad75.merged.bed  -e'FILTER="LowQual"' -o temp_gencodeFilt.vcf.gz -Oz
tabix temp_gencodeFilt.vcf.gz
#30X bed
bcftools filter temp_gencodeFilt.vcf.gz -R 38_30X_${EXOME}.bed -o ExFiltered_${SAMPLE}.vcf.gz -Oz
rm temp_gencodeFilt.vcf.gz

TRUTH=`readlink -f 38_30X_${EXOME}.vcf.gz`
VCF=`readlink -f ExFiltered_${SAMPLE}.vcf.gz`
BED=`readlink -f 38_30X_${EXOME}.bed`
THREADS=16

mkdir $OUTPUT_DIR

# Run hap.py
docker run -it \
-v "${TRUTH}":"${TRUTH}" \
-v "${VCF}":"${VCF}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${HGREF}":"${HGREF}" \
-v "${BED}":"${BED}" \
-v "/media/eservant/HDD_12To_31/ref/hg38_exclusion/":"/media/eservant/HDD_12To_31/ref/hg38_exclusion/" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
"${TRUTH}" \
"${VCF}" \
-r "${HGREF}" \
-o "${OUTPUT_DIR}/happy_${SAMPLE}" \
-f "${BED}" \
--pass-only \
--engine=xcmp \
--threads="${THREADS}"


