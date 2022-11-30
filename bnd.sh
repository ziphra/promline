# write BND mates

echo "writing BND mates..."

echo "grepping..."
grep '\[N' Promethion/run/Promethion_SQKLSK110_17102022/vc/sniffles/Y_6622CY001026_sniffles.vcf > grep.txt 

while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 2, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "[" $1 ":" $2 "[N", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep '\]N' Promethion/run/Promethion_SQKLSK110_17102022/vc/sniffles/Y_6622CY001026_sniffles.vcf > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 2, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "N[" $1 ":" $2 "[", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep 'N\]' Promethion/run/Promethion_SQKLSK110_17102022/vc/sniffles/Y_6622CY001026_sniffles.vcf > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 3, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "N]" $1 ":" $2 "]", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep 'N\[' Promethion/run/Promethion_SQKLSK110_17102022/vc/sniffles/Y_6622CY001026_sniffles.vcf > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 3, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "]" $1 ":" $2 "]N", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt 

if [ -s BND.txt ]; then echo "BNDs in BND.txt"; fi

cat Promethion/run/Promethion_SQKLSK110_17102022/vc/sniffles/Y_6622CY001026_sniffles.vcf BND.txt > unsorted.txt
grep -v '#' unsorted.txt | sort -k1,1V -k2,2n > sorted.txt
grep '#' unsorted.txt > header.txt
echo '##BND mates lines added' >> header.txt
cat header.txt sorted.txt > Promethion/run/Promethion_SQKLSK110_17102022/vc/sniffles/Y_6622CY001026_sniffles_BND.vcf 

rm header.txt
rm unsorted.txt
rm BND.txt
rm grep.txt 
rm sorted.txt