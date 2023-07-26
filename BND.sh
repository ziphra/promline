# write BND mates
VCF_SNF="/media/eservant/HDD_12To_31/Promethion/SQK-LSK114_062023/dorado3.0/vc/sniffles/17CY001482_sniffles2_t2t.vcf"
BASE="/media/eservant/HDD_12To_31/Promethion/SQK-LSK114_062023/dorado3.0"
FLAGS_sample="17CY001482_t2t"
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

cat $VCF_SNF BND.txt > unsorted.txt
grep -v '#' unsorted.txt | sort -k1,1V -k2,2n > sorted.txt
grep '#' unsorted.txt > header.txt
#echo '##BND mates lines added' >> header.txt
cat header.txt sorted.txt > $BASE/vc/sniffles/${FLAGS_sample}_sniffles_BND.vcf 

if [ -s $BASE/vc/sniffles/${FLAGS_sample}_sniffles_BND.vcf  ]; then echo "BND mating succesful."; else echo "WARNING: BND mating not successful";fi

# rm header.txt
# rm unsorted.txt
# rm BND.txt
# rm grep.txt 
# rm sorted.txt