#!/usr/bin/sudo bash

#  pipeline.sh
#  
#
#  Created by Euphrasie Servant on 17/06/2022.
#
basesh=`dirname $0`
. ~/miniconda3/etc/profile.d/conda.sh
pipeconfig=`realpath ${basesh}/pipeline.config`
. $pipeconfig


conda activate promline
echo ""
echo "==================== PIPELINE ====================="
echo ""
echo ""
echo `date`
echo ""
echo "sample: $FLAGS_sample "
echo ""
echo ""
mkdir -p $BASE

TORCH=$(cat <<EOF
import torch

print(torch.cuda.is_available())
EOF
)

GPU=`python -c "$TORCH"`
 
if [ $GPU = "True" ]
	then 
	echo "GPU(s) available, and will be use with Guppy, Dorado, and PMDV"
else
	echo "No GPU available"
fi

echo $0

##### 0 FAST5 to POD5 #####

if [[ "$FLAGS_pod5" = "" ]] 
then 
    echo ""
    echo "=========== no FAST5 to POD5 conversion============"
    echo ""   
elif [[ ! -d "$FLAGS_pod5" || ! "$(ls -A ${FLAGS_pod5})" ]] 
then
echo ""
echo "=========== FAST5 to POD5 ============"
echo ""
    echo -e "$FLAGS_pod5 is empty or doesn't exist: converting fast5 to pod5."
    pod5-convert-from-fast5 `find $FAST5 -name "*.fast5"` $FLAGS_pod5
else
    echo ""
    echo ""
    echo "No pod5 conversion."
    echo ""
fi



##### 1 BASECALLING #####


echo ""
echo "=========== Basecalling for ${FLAGS_sample} ============"
echo ""
echo ""


## GUPPY
if [[ "$FLAGS_basecalling" = "guppy" || "$FLAGS_basecalling" = "all" ]]
then 
    if [[ $GPU = "True" ]]
    then
        echo "GUPPY 6"
        echo ""

        mkdir $GUPPY

        time guppy_basecaller \
        -i $FAST5 \
        -s $GUPPY \
        --records_per_fastq 0 \
        -r \
        -c ${MODEL_GUPPY} \
        --device 'auto' \
        --compress_fastq \
        --num_callers ${FLAGS_threads} \
        --chunk_size 1000 \
        --gpu_runners_per_device 4 \
        --chunks_per_runner 512 \
        --disable_pings

        cat $GUPPY/pass/* > $FASTQ

    elif [[ $GPU = "True" ]]
    then
        echo "GUPPY 6"
        echo ""

        mkdir $GUPPY

        time guppy_basecaller \
        -i $FAST5 \
        -s $GUPPY \
        --records_per_fastq 0 \
        -r \
        -c ${MODEL_GUPPY} \
        --compress_fastq \
        --num_callers ${FLAGS_threads} \
        --disable_pings

        cat $GUPPY/pass/* > $FASTQ
    fi
    
fi

## DORADO

if [[ "$FLAGS_basecalling" = "dorado" || "$FLAGS_basecalling" = "all" ]]
then 
    echo ""
    echo "=========== DORADO ==========="
    echo ""

    mkdir ${BASE}/dorado/

    time ${dorado}/bin/dorado basecaller -b 512 ${dorado}/${MODEL_DORADO} $FLAGS_pod5/ | samtools view -Sh -@ 6 - > $DORADOBAM
fi

if [[ "$FLAGS_basecalling" = "none" ]]
then 
    echo ""
    echo "=========== no basecalling ==========="
    echo ""
fi


##### 2 MAPPING MMI ######

## MMI 

if [[ "$FLAGS_basecalling" = "guppy" || "$FLAGS_basecalling" = "all" ]]
then 
    echo ""
    echo "=========== MMI GUPPY============"
    echo ""

    mkdir $MMI

    time minimap2 -t $FLAGS_threads \
    -ax map-ont \
    $REFMMI \
    --MD \
    -Y \
    $FASTQ | samtools sort -o $BAM 

    samtools index $BAM
fi


## 2 MMI DORADO

if [[ "$FLAGS_basecalling" = "dorado" || "$FLAGS_basecalling" = "all" ]]
then
    echo ""
    echo "=========== MMI DORADO  ============"
    echo ""

    samtools bam2fq $DORADOBAM | $minimap2 -t 10 -ax map-ont --MD $REFMMI - | samtools sort -@ 4 -o $DORADOMMI

    samtools index $DORADOMMI
fi

## MMI2

if [[ "$FLAGS_basecalling" = "none" ]]
then 
    echo ""
    echo "=========== MMI  ============"
    echo ""

    mkdir $MMI
    
    if [ -d $FASTQSFOLD ]
    then
        cat $FASTQSFOLD/pass/* > $FASTQ
    else 
        FASTQ = "$FASTQSFOLD"
    fi

    time minimap2 -t $FLAGS_threads \
    -ax map-ont \
    $REFMMI \
    --MD \
    -Y \
    $FASTQ | samtools sort -o $BAM 

    samtools index $BAM
fi


##### QC ######


if [[ "$FLAGS_basecalling" = "guppy" || "$FLAGS_basecalling" = "all" ]]
then 
    echo ""
    echo "=========== QC GUPPY-MMI ==========="
    echo ""

    pycoQC -f $GUPPY/sequencing_summary.txt -a $BAM -o $BASE/${FLAGS_sample}_QC.html

elif [[ "$FLAGS_basecalling" = "dorado" || "$FLAGS_basecalling" = "all" ]]
then
    echo ""
    echo "=========== QC DORADO-MMI ==========="
    echo ""
    echo "coming soon..."
else
    echo ""
    echo "=========== QC MMI ==========="
    echo ""
    if [ -d $FASTQSFOLD ]
    then
        pycoQC -f $FASTQSFOLD/sequencing_summary.txt -a $BAM -o $BASE/${FLAGS_sample}_QC.html
    else 
        pycoQC -f $BASE/sequencing_summary.txt -a $BAM -o $BASE/${FLAGS_sample}_QC.html
    fi
fi




##### VC #####


echo ""
echo "=========== VC ==========="
echo ""


#### SV ####

### sniffles ###

echo "SV with sniffles"
echo ""

mkdir -p $BASE/vc/sniffles
time sniffles -i $BAM \
	--vcf $VCF_SNF \
	--tandem-repeats ${basesh}/human_GRCh38_no_alt_analysis_set.trf.bed \
	--reference $REF \
    --long-del-coverage 5 \
    --long-dup-coverage 0.5 \
	-t $FLAGS_threads 

time sniffles -i $BAM \
	--vcf $NOQC_VCF_SNF \
	--tandem-repeats ${basesh}/human_GRCh38_no_alt_analysis_set.trf.bed \
	--reference $REF \
    --no-qc \
	-t $FLAGS_threads 


# write BND mates

echo "writing BND mates..."

echo "grepping..."
grep '\[N' $VCF_SNF > grep.txt 

while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 2, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "[" $1 ":" $2 "[N", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep '\]N' $VCF_SNF > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 2, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "N[" $1 ":" $2 "[", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep 'N\]' $VCF_SNF > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 3, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "N]" $1 ":" $2 "]", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep 'N\[' $VCF_SNF > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 3, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "]" $1 ":" $2 "]N", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt 

if [ -s BND.txt ]; then echo "BNDs in BND.txt"; fi

cat $BASE/vc/sniffles/${FLAGS_sample}_sniffles.vcf BND.txt > unsorted.txt
grep -v '#' unsorted.txt | sort -k1,1V -k2,2n > sorted.txt
grep '#' unsorted.txt > header.txt
#echo '##BND mates lines added' >> header.txt
cat header.txt sorted.txt > $BASE/vc/sniffles/${FLAGS_sample}_sniffles_BND.vcf 

if [ -s $BASE/vc/sniffles/${FLAGS_sample}_sniffles_BND.vcf  ]; then echo "BND mating succesful."; else echo "WARNING: BND mating not successful";fi

rm header.txt
rm unsorted.txt
rm BND.txt
rm grep.txt 
rm sorted.txt



#### small variants calling ####
#### norm ####
vtf()
{
vcf=${1%.*}
vt decompose_blocksub ${1} | \
vt decompose -s - | \
vt normalize -r $REF -o ${vcf}_raw.norm.vcf.gz -
tabix ${vcf}_raw.norm.vcf.gz
}


### pmdv ###

if [[ "$FLAGS_snp_caller" = "pmdv" || "$FLAGS_snp_caller" = "all" ]]
then 
    if [[ $GPU = "True" ]]
    then 
        echo ""
        echo ""
        echo "SMALL VARIANTS with PMDV - GPU"
        echo "" 

        ## Create local directory structure
        mkdir -p "${OUTPUT_DIR_PMDV}"
        BAMDIR=`dirname $BAM`
        echo $BAMDIR

        docker run --ipc=host \
	        --gpus all \
	        -v "${BAMDIR}":"${BAMDIR}" \
	        -v "${OUTPUT_DIR_PMDV}":"${OUTPUT_DIR_PMDV}" \
	        -v "${REF}":"${REF}" \
	        kishwars/pepper_deepvariant:r0.8-gpu \
	        run_pepper_margin_deepvariant call_variant \
	        -o "${OUTPUT_DIR_PMDV}" \
	        -b "${BAM}" \
	        -f "${REF}" \
	        -p "${PMDV_PREFIX}" \
	        -t "${FLAGS_threads}" \
	        -g \
	        ${MODEL_PMDV}

        echo "remove refcalls with bcftools..."
        echo ""

        bcftools filter -e'FILTER="refCall"' ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}.vcf.gz -o ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz -Oz
        tabix ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz
        vtf ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz

    elif [[ $GPU = "False" ]]
    then 
        echo ""
        echo ""
        echo "SMALL VARIANTS with PMDV - no GPU"
        echo "${MODEL_PMDV}"

        ## Create local directory structure
        mkdir -p "${OUTPUT_DIR_PMDV}"
        BAMDIR=`dirname $BAM`
        echo $BAMDIR
        echo 
        docker run --ipc=host \
	        -v "${BAMDIR}":"${BAMDIR}" \
	        -v "${OUTPUT_DIR_PMDV}":"${OUTPUT_DIR_PMDV}" \
	        -v "${REF}":"${REF}" \
	        kishwars/pepper_deepvariant:r0.8 \
	        run_pepper_margin_deepvariant call_variant \
	        -o "${OUTPUT_DIR_PMDV}" \
	        -b "${BAM}" \
	        -f "${REF}" \
	        -p "${PMDV_PREFIX}" \
	        -t "{FLAGS_threads}" \
	        -g \
	        "${MODEL_PMDV}"

        echo "remove refcalls with bcftools..."
        echo ""

        bcftools filter -e'FILTER="refCall"' ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}.vcf.gz -o ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz -Oz
        tabix ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz
        vtf ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz
    fi
fi

### clair3 ###

if [[ "$FLAGS_snp_caller" = "clair3" || "$FLAGS_snp_caller" = "all" ]]
then 

    ### CLAIR3 ###
    echo ""
    echo ""
    echo "SMALL VARIANTS with CLAIR3"
    echo ""

    run_clair3.sh \
	    --bam_fn=${BAM} \
	    --ref_fn=${REF} \
	    --threads=${FLAGS_threads} \
	    --platform="ont" \
	    --model_path=${CONDA_PREFIX}/bin/models/${MODEL_CLAIR} \
	    --output=${CLAIR_OUT} \
	    --remove_intermediate_dir

    # python /home/euphrasie/miniconda3/envs/promline/bin/clair3.py SwitchZygosityBasedOnSVCalls \
    #   --bam_fn ${BAM} \
    #   --clair3_vcf_input ${CLAIR_OUT}/merge_output.vcf.gz \
    #   --sv_vcf_input $VCF_SNF \
    #   --vcf_output ${CLAIR_OUT}_merge_output_switch.vcf \
    #   --threads ${FLAGS_threads}

    vtf ${CLAIR_OUT}/merge_output.vcf
fi



echo ""
echo ""
echo "job done"
echo ""
echo ""