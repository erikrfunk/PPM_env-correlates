#!/bin/bash

indCount=$1 # number of individuals total
fin=$2
fout=$3

zcat $fin | cut -f1,2,3 > $fout.txt
sed -i 's/\(HiC_scaffold_[0-9]*\)_/\1 /g' $fout.txt
sed -i 's/marker/chr\tpos/g' $fout.txt

for i in $(seq 1 $indCount); do
echo $i
der=$((($i+1)*3)) # Add 1 to skip the pos and allele cols
het=$(($der-1))
zcat $2 | cut -f $der,$het | awk 'FNR==1 {print $1} FNR>1 {print $1/2+$2}' > $fin"_IND"$i.txt
done 

paste $fout".txt" $fin"_IND"* > $fout"_joined.txt"
