#!/bin/bash
#$ -cwd -V
#$ -pe smp 5
#$ -l h_vmem=20G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

cufflinks -p 4 -q --max-bundle-frags 30000000 -u -b /opt/databases/genomes/hg19/bt2_th_hg19_ucsc/genome.fa \
		  --GTF-guide /opt/databases/genomes/hg19/bt2_th_hg19_ucsc/genes.gtf \
		  -o ../Cufflinks/$1 $1/accepted_hits.bam