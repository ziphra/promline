bash promline/pipeline.sh -w Promethion_SQKLSK110_17102022 -s 6622CY001026 -r /media/euphrasie/DATA/reference_genome/hg38_exclusion/GRCh38.p13.GRC_exclusions_T2Tv2.fa -f /media/euphrasie/Expansion/Promethion_SQKLSK110_17102022/6622CY001026/20221024_1201_2C_PAM61860_63f212a3/fast5 -p /media/euphrasie/Expansion/Promethion_SQKLSK110_17102022/pod5 -c all -b all -m r9 -t /home/euphrasie/Documents/lr_test3/sniffles/human_GRCh38_no_alt_analysis_set.trf.bed 2>&1 | tee Promethion_SQKLSK110_17102022/log.txt


bash promline/annotation.sh -w Promethion_SQKLSK110_17102022 -r /media/euphrasie/DATA/reference_genome/hg38_exclusion/GRCh38.p13.GRC_exclusions_T2Tv2.fa -v /media/euphrasie/Alienware_May202/Promethion_SQKLSK110_17102022/vc/sniffles/6622CY001026_sniffles_BND.vcf -t /media/euphrasie/Alienware_May202/Promethion_SQKLSK110_17102022/vc/clair3/merge_output_raw.norm.vcf.gz 2>&1 | tee Promethion_SQKLSK110_17102022/logann.txt


bash promline/annotation.sh -w TRASH -r /media/euphrasie/DATA/reference_genome/hg38_exclusion/GRCh38.p13.GRC_exclusions_T2Tv2.fa -v /media/euphrasie/Alienware_May202/TRASH/subprom.vcf.gz -t /media/euphrasie/Alienware_May202/TRASH/subprom.vcf.gz 


bash promline/pipeline.sh -w run/Promethion_SQKLSK110_1411_17102022 \
-s "6622CY001026bis" \
-r /media/euphrasie/DATA/reference_genome/hg38_exclusion/GRCh38.p13.GRC_exclusions_T2Tv2.fa \
-f /media/euphrasie/Expansion/Promethion_SQK-LSK110_14112022/6622CY001026/20221115_1050_1H_PAM59472_9c810272/fast5 \
-p /media/euphrasie/Expansion/Promethion_SQK-LSK110_14112022/pod5 \
-q run/Promethion_SQKLSK110_1411_17102022/6622CY001026.fastq.gz \
-c all \
-b none \
-m r9 \
-t /home/euphrasie/Documents/lr_test3/sniffles/human_GRCh38_no_alt_analysis_set.trf.bed 2>&1 | tee Promethion_SQKLSK110_1411_17102022/log.txt


bash Promethion/promline/pipeline.sh -w test \
-s "sampletest" \
-r /media/euphrasie/DATA/reference_genome/hg38_exclusion/GRCh38.p13.GRC_exclusions_T2Tv2.fa \
-f /media/euphrasie/Expansion/subset5 \
-p /media/euphrasie/Expansion/subsetpod5 \
-c pmdv \
-b guppy \
-m r9 \
-t 16 | tee log.txt


promline -w /media/eservant/HDD_12To_3/Promethion/Promethion_SQKLSK110_fusion \
-s "6622CY001026_dell" \
-r /media/eservant/HDD_12To_3/ref/hg38_exclusion/GRCh38.p13.GRC_exclusions_T2Tv2.fa \
-f /media/eservant/Expansion/Promethion_SQKLSK110_fusion \
-p /media/eservant/HDD_12To_3/Promethion/Promethion_SQKLSK110_fusion/pod5 \
-c all \
-b guppy \
-m r9 \
-t 40 | tee Promethion_SQKLSK110_1411_17102022/log.txt


promline -w /media/eservant/HDD_12To_3/Promethion/Promethion_SQKLSK110_fusion \
-s "6622CY001026_dell" \
-r /media/eservant/HDD_12To_3/ref/hg38_exclusion/GRCh38.p13.GRC_exclusions_T2Tv2.fa \
-f /media/eservant/Expansion/Promethion_SQKLSK110_fusion \
-p /media/eservant/HDD_12To_3/Promethion/Promethion_SQKLSK110_fusion/pod5 \
-c all \
-b all \
-m r9 \
-t 40 | tee Promethion_SQKLSK110_1411_17102022/log.txt


promline -w /media/eservant/HDD_12To_3/Promethion/Promethion_SQKLSK110_fusion -s "6622CY001026_dell" -r /media/eservant/HDD_12To_3/ref/hg38_exclusion/GRCh38.p13.GRC_exclusions_T2Tv2.fa -f /media/eservant/Expansion/Promethion_SQKLSK110_fusion -p /media/eservant/HDD_12To_3/Promethion/Promethion_SQKLSK110_fusion/pod5 -c all -b all -m r9 -t 40 | tee /media/eservant/HDD_12To_3/Promethion/Promethion_SQKLSK110_fusion/log.txt
