#!/bin/bash
#$ -cwd -V
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile
$out=$1'.bedgraph'
genomeCoverageBed -bg -ibam $1/accepted_hits.bam -g hg19.genome > bedgraph/$out

