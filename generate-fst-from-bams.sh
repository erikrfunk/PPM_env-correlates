#!/bin/bash

win=$1
step=$2
out=$3

angsd -b ../bam1 -ref /home/centos/USS/conservation-genetics/PPM/HiC_assembly/PPM_HiRise_rh_HiC_final.fasta \
   -anc /home/centos/USS/conservation-genetics/PPM/HiC_assembly/PPM_HiRise_rh_HiC_final.fasta \
   -P 6 -dosaf 1 -GL 1 \
   -sites /home/centos/USS/erik/PPM/joined_historical_and_contemporary/realigned/newSnpCalls/SnpList_minInd140_200DP1500_0.02MAF.txt \
   -rf ../../14_overlap_gene_coords_regions.txt \
   -minQ 20 -out bam1

angsd -b ../bam2 -ref /home/centos/USS/conservation-genetics/PPM/HiC_assembly/PPM_HiRise_rh_HiC_final.fasta \
    -anc /home/centos/USS/conservation-genetics/PPM/HiC_assembly/PPM_HiRise_rh_HiC_final.fasta \
    -P 6 -dosaf 1 -GL 1 \
    -sites /home/centos/USS/erik/PPM/joined_historical_and_contemporary/realigned/newSnpCalls/SnpList_minInd140_200DP1500_0.02MAF.txt \
    -rf ../../14_overlap_gene_coords_regions.txt \
    -minQ 20 -out bam2

#calculate the 2dsfs prior
realSFS bam1.saf.idx bam2.saf.idx -fold 1 >bam1.bam2.ml
#prepare the fst for easy window analysis etc
realSFS fst index bam1.saf.idx bam2.saf.idx -sfs bam1.bam2.ml -fstout fst_intermediate
#get the global estimate
realSFS fst stats fst_intermediate.fst.idx 
#below is not tested that much, but seems to work
realSFS fst stats2 fst_intermediate.fst.idx -win $win -step $step > $out
