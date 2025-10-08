#!/bin/bash
reads=$1
bam=$2
filename=$(basename $bam .bam)

echo 'writing downsampled bam to:' ${filename}_subsampled${reads}.bam

fraction=$(samtools idxstats $bam | cut -f3 | awk -v ct=$reads 'BEGIN {total=0} {total += $1} END {print ct/total}')
samtools view -b -s ${fraction} $bam > ${filename}_subsampled${reads}.bam
