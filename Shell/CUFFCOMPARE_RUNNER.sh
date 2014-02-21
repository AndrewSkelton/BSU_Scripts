#!/bin/bash
#$ -cwd -V
#$ -pe smp 2
#$ -l h_vmem=25G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cuffcompare -s /opt/databases/genomes/hg19/bt2_th_hg19_ucsc/genome.fa -r /opt/databases/genomes/hg19/bt2_th_hg19_ucsc/genes.gtf ../Merge/merged.gtf

