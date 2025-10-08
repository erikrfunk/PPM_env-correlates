#!/bin/bash

start=$1
end=$2

for i in $(eval echo "{$start..$end}"); do
echo "Starting iteration" ${i} "..."
mkdir iteration${i}
cd iteration${i}

grep pxrecal ../coded_lcwps_inds_bams.txt | shuf | awk 'FNR<4 {print $0 >> "bam1"} FNR>=4 {print $0 >> "bam2"}'
grep -v pxrecal ../coded_lcwps_inds_bams.txt | shuf | awk 'FNR<43 {print $0 >> "bam1"} FNR >=43 {print $0 >> "bam2"}'

mkdir fst
cd fst

../../generate-fst-from-bams.sh 15000 10000 iteration${i}_fst15kbWin10kbSli_14cands.txt

cd ../..
done
