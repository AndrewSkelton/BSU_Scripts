#!/bin/bash
#$ -cwd -V
#$ -pe smp 5
#$ -l h_vmem=10G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

tophat -r $3 --phred64-quals --tmp-dir $TMPDIR -p 4 --transcriptome-index ~/WORKING_DATA/RNA_Seq/Analysis/transcriptome_index_hg19/known \
		--b2-very-sensitive -o ~/WORKING_DATA/RNA_Seq/Analysis/Tophat/$4 -G /opt/databases/genomes/hg19/bt2_th_hg19_ucsc/genes.gtf \
		/opt/databases/genomes/hg19/bt2_th_hg19_ucsc/genome	$1 $2