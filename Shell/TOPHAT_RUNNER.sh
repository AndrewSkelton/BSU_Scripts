#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=50G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

tophat -r $3 --mate-std-dev $4 --phred64-quals --tmp-dir $TMPDIR -p 10 --transcriptome-index /opt/databases/genomes/Ensembl/transcriptome_index/Ensembl \
		--b2-very-sensitive -o ./Tophat/$5 -G /opt/databases/genomes/Ensembl/Homo_sapiens.GRCh37.Ensembl.gtf \
		/opt/databases/genomes/Ensembl/Homo_sapiens.GRCh37.Ensembl	$1 $2

mv ./Tophat/$5/accepted_hits.bam ./Tophat/accepted_hits/$5'.bam'
samtools index ./Tophat/accepted_hits/$5'.bam'
