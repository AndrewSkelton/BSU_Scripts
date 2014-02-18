#!/bin/bash
#$ -cwd -V
#$ -pe smp 6
#$ -l h_vmem=10G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

tophat -G /opt/databases/genomes/hg19/andrew/genes.gtf --tmp-dir $TMPDIR \
    --transcriptome-index=$TMPDIR/transcriptome_index_hg19/known \
    /opt/databases/bowtie2/hg19

mv $TMPDIR/transcriptome_index_hg19 ~/WORKING_DATA

