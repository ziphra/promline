#!/usr/bin/sudo bash

#  pipeline.sh
#  
#
#  Created by Euphrasie Servant on 17/06/2022.
#
. ~/miniconda3/etc/profile.d/conda.sh


. /media/euphrasie/Alienware_May202/Promethion/promline/annotation.flag



vepf() 
{

echo ""
echo "SMALL VARIANTS ANNOTATIONS WITH VEP"
echo ""

conda activate VEP

mkdir ${FLAGS_base}/VEP/

vcflong=`basename $1`
vcf=VEP_${vcflong}

time vep --cache \
	--offline \
	--fork 16 \
	--species homo_sapiens \
	--assembly GRCh38 \
	-i ${1} \
 	-o  ${FLAGS_base}/VEP/${vcf} \
	--fasta $FLAGS_ref \
	--use_given_ref \
	--format "vcf" \
	--compress_output bgzip \
	--vcf \
	--ccds \
	--shift_hgvs=1 \
	--hgvs --hgvsg \
	--symbol \
	--numbers \
	--canonical \
	--tsl \
	--appris \
	--pubmed \
	--variant_class \
	--mane \
    --regulatory \
	--pick \
	--gene_phenotype \
	--pick_order rank,mane,biotype,canonical \
	--plugin UTRannotator,$HOME/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
	--plugin SpliceAI,snv=/media/euphrasie/Alienware_May202/Annotations/spliceai/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/media/euphrasie/Alienware_May202/Annotations/spliceai/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz \
	--plugin REVEL,/media/euphrasie/Alienware_May202/Annotations/revel/new_tabbed_revel_grch38.tsv.gz \
	--plugin CADD,/media/euphrasie/Alienware_May202/Annotations/CADD_1.6_hg38/whole_genome_SNVs.tsv.gz,/media/euphrasie/Alienware_May202/Annotations/CADD_1.6_hg38/gnomad.genomes.r3.0.indel.tsv.gz \
    --custom /media/euphrasie/Alienware_May202/Annotations/clinvar/apr2022/GRCh38/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN \
    --custom /media/euphrasie/Alienware_May202/Annotations/OMIM/apr2022/omiminh_apr2022.bed.gz,OMIMINH,bed,overlap \
    --custom /media/euphrasie/Alienware_May202/Annotations/OMIM/apr2022/omimphen_apr2022.bed.gz,OMIMPHEN,bed,overlap \
    --custom /media/euphrasie/Alienware_May202/Annotations/HGMD/HGMD_2021.4/hgmd_pro_2021.4_hg38.vcf.gz,HGMD,vcf,exact,PHEN \
    --custom /media/euphrasie/Alienware_May202/Annotations/genesACC/genesACC_Aug20192.bed.gz,ACC,bed,overlap \
	--custom /media/euphrasie/Alienware_May202/Annotations/gnomAD_merge_exom2.1.1_genom3.1.2/gnomad38_211_312.vcf.gz,gnomAD,vcf,exact,0,homhem,het


tabix ${FLAGS_base}/VEP/${vcf}

xlsx=${vcf::-7}.txt
# change gnomad fields from string to float
bcftools +split-vep -c gnomAD_het,gnomAD_homhem:FLOAT ${FLAGS_base}/VEP/${vcf} -o ${FLAGS_base}/VEP/float_${vcf::-3}

sed -i 's/^##INFO=<ID=gnomAD_het,Number=.,Type=String/##INFO=<ID=gnomAD_het,Number=.,Type=Float/g' ${FLAGS_base}/VEP/float_${vcf::-3}
sed -i 's/^##INFO=<ID=gnomAD_homhem,Number=.,Type=String/##INFO=<ID=gnomAD_homhem,Number=.,Type=Float/g' ${FLAGS_base}/VEP/float_${vcf::-3}


bcftools filter -i'FILTER="PASS"' ${FLAGS_base}/VEP/float_${vcf::-3} | bcftools filter -e 'DP<5 || QUAL<10 || gnomAD_het>2000 || gnomAD_homhem>1' -Oz -o ${FLAGS_base}/VEP/FILT_${vcf::-3}

bcftools +split-vep -f '%CHROM %POS %CSQ\n' -A tab -d ${FLAGS_base}/VEP/float_${vcf::-3} > noHead_${xlsx}
bcftools +split-vep ${FLAGS_base}/VEP/float_${vcf::-3} -l | awk '{print $2}' | tr "\n" "\t" > header
cat header noHead_${xlsx} > $xlsx

#rm ${FLAGS_base}/VEP/${vcf} ${FLAGS_base}/VEP/${vcf}.tbi header 
conda deactivate

}


mkdir -p $FLAGS_base

vepf ${FLAGS_small_vcf}


### ANNOTSV ###

conda deactivate  

echo "STRUCTURAL VARIANTS ANNOTATIONS WITH ANNOTSV"
echo ""
mkdir -p ${FLAGS_base}/ANNOTSV

$ANNOTSV/bin/AnnotSV -SvinputFile ${FLAGS_sv_vcf} \
    -outputDir ${FLAGS_base}/ANNOTSV \
	-candidateSnvIndelFiles ${FLAGS_small_vcf} \
	-snvIndelFiles ${FLAGS_small_vcf}

$KNOTANNOTSV --genomeBuild hg38 \
    --annotSVfile ${FLAGS_base}/ANNOTSV/*.annotated.tsv \
    --outDir ${FLAGS_base}/ANNOTSV/ \
    --configFile $ANNOTSV_YAML

$KNOTANNOTSV2XL --genomeBuild hg38 \
    --annotSVfile ${FLAGS_base}/ANNOTSV/*.annotated.tsv \
    --outDir ${FLAGS_base}/ANNOTSV/ \
    --configFile $ANNOTSV_YAML
