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

fsam_t=$(($FLAGS_threads/4))
isam_t=$(printf "%.0f" "$fsam_t")

## DORADO

 if [[ $FLAGS_modified -eq ${FLAGS_TRUE} ]]  && [[ $FLAGS_basecalling -eq ${FLAGS_TRUE} ]]
     then
     echo ""
     echo "=========== DORADO MODIFIED BASES ==========="
     echo ""
     mkdir ${BASE}/dorado/                     
     echo ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO}
     echo ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO_MOD}
     time dorado basecaller -r 4 ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} $FLAGS_pod5/ --modified-bases-models ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO_MOD} | samtools view -bSh -@ isam_t - > $DORADOBAM
fi

# elif [[ $FLAGS_modified -eq ${FLAGS_FALSE} ]] && [[ $FLAGS_basecalling -eq ${true} ]]
#     then 
#     echo ""
#     echo "=========== DORADO ==========="
#     echo ""

#     mkdir ${BASE}/dorado/

#     time dorado basecaller -r 4 -b 256 ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} $FLAGS_pod5/ | samtools view -bSh -@ $isam_t - > $DORADOBAM

# fi 



if [[ $FLAGS_basecalling -eq ${true} ]] && [[ $FLAGS_duplex -eq ${FLAGS_FALSE} ]] && [[ $FLAGS_modified -eq ${FLAGS_FALSE} ]]
    then 
    echo ""
    echo "=========== DORADO ==========="
    echo ""

    mkdir ${BASE}/dorado/

    time dorado basecaller -r 4 -b 256 ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} $FLAGS_pod5/ | samtools view -bSh -@ $isam_t - > $DORADOBAM

fi 






if [[ $FLAGS_duplex -eq ${FLAGS_TRUE} ]] 
        then
        echo ""
        echo "=========== DORADO DUPLEX ==========="
        echo ""
        mkdir ${BASE}/dorado/ 
        mkdir ${BASE}/dorado/duplex                    
        echo ${CONDA_PREFIX}/bin/dorado_models/dna_r10.4.1_e8.2_${FLAGS_bps}bps_fast@v4.0.0
        echo $FLAGS_pod5/
        # simplex basecall with dorado
        time dorado basecaller ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} $FLAGS_pod5/ --emit-moves | samtools sort -@ $isam_t -o  ${BASE}/dorado/duplex/unmapped_reads_with_moves.bam
        samtools index ${BASE}/dorado/duplex/unmapped_reads_with_moves.bam
        
        # find duplex pairs 
        time duplex_tools pair --output_dir ${BASE}/dorado/duplex/pairs_from_bam ${BASE}/dorado/duplex/unmapped_reads_with_moves.bam

        #find addditional duplex in non split reads
        time duplex_tools split_pairs ${BASE}/dorado/duplex/unmapped_reads_with_moves.bam $FLAGS_pod5/ ${BASE}/dorado/duplex/pod5s_splitduplex/
        
        cat ${BASE}/dorado/duplex/pod5s_splitduplex/*_pair_ids.txt > ${BASE}/dorado/duplex/split_duplex_pair_ids.txt
        
        #stereo dupelx all
        time dorado duplex ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} $FLAGS_pod5/ --pairs ${BASE}/dorado/duplex/pairs_from_bam/pair_ids_filtered.txt | samtools sort -@ $isam_t -o  ${BASE}/dorado/duplex/duplex_orig.bam
        samtools index ${BASE}/dorado/duplex/duplex_orig.bam
        
        time dorado duplex ${CONDA_PREFIX}/bin/dorado_models/${MODEL_DORADO} ${BASE}/dorado/duplex/pod5s_splitduplex/ --pairs ${BASE}/dorado/duplex/split_duplex_pair_ids.txt | samtools sort -@ $isam_t -o  ${BASE}/dorado/duplex/duplex_splitduplex.bam
        samtools index ${BASE}/dorado/duplex/duplex_splitduplex.bam

        samtools cat ${BASE}/dorado/duplex/duplex_orig.bam ${BASE}/dorado/duplex/duplex_splitduplex.bam ${BASE}/dorado/duplex/unmapped_reads_with_moves.bam -@ $isam_t -o $DORADOBAM
fi


if [[ $FLAGS_basecalling -eq ${FLAGS_FALSE} ]]
then 
    echo ""
    echo "=========== no basecalling ==========="
    echo ""
fi

##### 2 MAPPING MMI ######

## 2 MMI DORADO

if [[ $FLAGS_basecalling -eq ${FLAGS_TRUE} ]] && [[ $FLAGS_alignment -eq ${FLAGS_TRUE} ]]
then
    echo ""
    echo "=========== MMI DORADO  ============"
    echo ""
    mkdir $MMI
    
    samtools bam2fq $DORADOBAM | minimap2 -Y -t 10 -ax map-ont --MD -y $REFMMI - | samtools sort -@ $isam_t -o $DORADOMMI

    samtools index $DORADOMMI
    BAM=${DORADOMMI}
fi

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
    samtools view $DORADOBAM | awk -v OFS='\t' '{print $1, $12}' | sed 's/qs:i://g' > qscore.txt
    echo -e 'read_id\tmean_qscore_template' | cat - qscore.txt > h_qscore.txt
    samtools view $DORADOBAM | awk -v OFS='\t' '{print $1, $10}' > read.txt
    echo -e 'read_id\tread' | cat - read.txt > h_read.txt
    python `realpath ${basesh}/sequencingsumm.py` $FAST5
    awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$17,$16}' presequencing_summary.txt > $BASE/sequencing_summary.txt 
    rm presequencing_summary.txt h_qscore.txt h_read.txt qscore.txt read.txt
    pycoQC -f $BASE/sequencing_summary.txt -a $BAM -o $BASE/${FLAGS_sample}_QC.html
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

    vtf ${CLAIR_OUT}/merge_output.vcf.gz
fi



echo ""
echo ""
echo "job done"
echo ""
echo ""