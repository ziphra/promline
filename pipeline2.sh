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
	echo "GPU(s) available, and will be use with Dorado, and PMDV"
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
    pod5-convert-from-fast5 -p $FLAGS_threads `find $FAST5 -name "*.fast5"` $FLAGS_pod5
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

fsam_t=$(($FLAGS_threads/4))
isam_t=$(printf "%.0f" "$fsam_t")

## DORADO
## v03 simplex, mod bases 
 if [[ $FLAGS_modified -eq ${FLAGS_TRUE} ]]  && [[ $FLAGS_basecalling -eq ${FLAGS_TRUE} ]]
     then
     echo ""
     echo "=========== DORADO MODIFIED BASES ==========="
     echo ""
     mkdir ${BASE}/dorado/                     
     echo ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO}
     time dorado basecaller ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} $FLAGS_pod5/ --modified-bases 5mCG_5hmCG -r  > ${DORADOBAM}
fi

## v03 simplex
if [[ $FLAGS_basecalling -eq ${true} ]] && [[ $FLAGS_duplex -eq ${FLAGS_FALSE} ]] && [[ $FLAGS_modified -eq ${FLAGS_FALSE} ]]
    then 
    echo ""
    echo "=========== DORADO ==========="
    echo ""

    mkdir ${BASE}/dorado/
    # simplex
    echo ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO}
    time dorado basecaller ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} $FLAGS_pod5/ -r > ${DORADOBAM}

fi 



## v03 duplex
if [[ $FLAGS_duplex -eq ${FLAGS_TRUE} ]] 
    then
    echo ""
    echo "=========== DORADO DUPLEX ==========="
    echo ""
    mkdir ${BASE}/dorado/
    echo $FLAGS_pod5/
    echo ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO}
    echo ${MODEL_DORADO}
    dorado duplex ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} ${FLAGS_pod5} -r -t $FLAGS_threads > ${DORADOBAM}

fi


if [[ $FLAGS_basecalling -eq ${FLAGS_FALSE} ]]
then 
    echo ""
    echo "=========== no basecalling ==========="
    echo ""
fi

dorado summary $DORADOBAM > ${BASE}/sequencing_summary.txt

##### 2 MAPPING MMI ######

## MMI2

if [[ "$FLAGS_basecalling" -eq ${FLAGS_FALSE} ]] && [[ $FLAGS_alignment -eq ${FLAGS_TRUE} ]]
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
    -Y -y \
    $FASTQ | samtools sort -@ $isam_t -o $BAM 

    samtools index $BAM
fi


if [[ $FLAGS_alignment -eq ${FLAGS_FALSE} ]]
then 
    echo ""
    echo "=========== NO ALIGNMENT  ============"
    echo ""
fi



if [[ $FLAGS_alignment -eq ${FLAGS_TRUE} ]]
then 
    mkdir $MMI
    echo ""
    echo "=========== MMI  ============"
    echo ""
    echo $REFMMI
    time samtools fastq -T MM,ML $DORADOBAM | minimap2 -t $FLAGS_threads \
    -ax map-ont \
    $REFMMI \
    --MD \
    -Y -y - | samtools sort -@ $isam_t -o $BAM 

    samtools index $BAM
fi

if [ -n "$FLAGS_bam" ]
then
    BAM=${FLAGS_bam}
fi

##### QC ######
#### SEQUENCING SUMM #####

echo ""
echo "=========== QC MMI ==========="
echo ""

if [ -n "$FASTQSFOLD" ]
then
    echo "fastqfold"
    echo $FASTQSFOLD
    pycoQC -f $FASTQSFOLD/sequencing_summary.txt -a $BAM -o $BASE/${FLAGS_sample}_QC.html
elif [ -n "${FLAGS_summary}" ]
then 
    pycoQC -f $FLAGS_summary -a $BAM -o $BASE/${FLAGS_sample}_QC.html
else
    pycoQC -f $BASE/sequencing_summary.txt -a $BAM -o $BASE/${FLAGS_sample}_QC.html
fi




##### VC #####
echo ""
echo "=========== SV calling ==========="
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
# conda deactivate 
# conda activate clair3

if [[ "$FLAGS_snp_caller" = "clair3" || "$FLAGS_snp_caller" = "all" ]]
then 

    ### CLAIR3 ###
    echo ""
    echo ""
    echo "SMALL VARIANTS with CLAIR3 AND PHASING"
    echo ""

    BAMDIR=`dirname $BAM`
    REFDIR=`dirname $REF`

    docker run -it \
        -v ${REFDIR}:${REFDIR} \
        -v ${CLAIR_OUT}:${CLAIR_OUT} \
        -v ${BAMDIR}:${BAMDIR} \
        -v ${CONDA_PREFIX}/bin/clair3_models/${MODEL_CLAIR}:${CONDA_PREFIX}/bin/clair3_models/${MODEL_CLAIR} \
        hkubal/clair3:v1.0.4 \
        /opt/bin/run_clair3.sh \
        --bam_fn=${BAM} \
        --ref_fn=${REF} \
        --threads=$FLAGS_threads \
        --platform="ont" \
        --model_path="${CONDA_PREFIX}/bin/clair3_models/$MODEL_CLAIR" \
        --output=${CLAIR_OUT} \
        --enable_phasing \
        --use_whatshap_for_final_output_haplotagging

    vtf ${CLAIR_OUT}/merge_output.vcf.gz
 fi

conda deactivate


#### CNV ####
### cnvpytor ###

conda activate cnvpytor

mkdir $BASE/vc/cnvpytor/

#read depth information from bam 
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -j $FLAGS_threads -rd $BAM


# root object
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -j $FLAGS_threads -his 1000 5000

#partition
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -j $FLAGS_threads -partition 1000 5000

#call
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -j $FLAGS_threads -call 1000 > $BASE/vc/cnvpytor/calls.1000.tsv
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -j $FLAGS_threads -call 5000 > $BASE/vc/cnvpytor/calls.5000.tsv


#import SNP from vcf 
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -j $FLAGS_threads -snp ${CLAIR_OUT}/merge_output.vcf.gz -sample SAMPLE 
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -j $FLAGS_threads -baf 1000 5000


# interactive mode 
echo "interactive mode..."
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -view 1000 <<ENDL
set Q0_range -1 0.5
set p_range 0 0.0001
set p_N 0 0.5
set size_range 50000 inf
set print_filename $BASE/vc/cnvpytor/cnvpytor_1000.xlsx
set annotate
print calls
ENDL

echo "interactive mode..."
cnvpytor -root $BASE/vc/cnvpytor/file.pytor -view 5000 <<ENDL
set Q0_range -1 0.5
set p_range 0 0.0001
set p_N 0 0.5
set size_range 50000 inf
set print_filename $BASE/vc/cnvpytor/cnvpytor_5000.xlsx
set annotate
print calls
ENDL

echo ""
echo ""
echo "job done"
echo ""
echo ""