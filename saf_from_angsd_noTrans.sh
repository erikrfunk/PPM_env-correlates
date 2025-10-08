#!/bin/bash

samples=$1

while read -r id; do
echo "Starting: " ${id}
angsd -i ${id}_*.bam -minQ 20 -minmapQ 20 -C 50 -noTrans 1 -doSaf 1 -GL 1 -out ${id}_He_noTrans -P 16 -anc /home/centos/USS/erik/PPM/PPM_reference.fasta -ref /home/centos/USS/erik/PPM/PPM_reference.fasta
done<$samples


