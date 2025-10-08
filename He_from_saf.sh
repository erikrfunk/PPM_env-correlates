#!/bin/bash

samples=$1

while read -r id; do
echo "Starting: " ${id}
realSFS ${id}_*.saf.idx -maxiter 2000 -tole 1e-16 > ${id}_He_noTrans.ml
awk -v samp=${id} '{print samp"\t"$2/($2+$1)}' ${id}_He_noTrans.ml >> heterozygosity.txt
done < $samples


