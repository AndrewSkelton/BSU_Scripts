#!/bin/bash
#$ -cwd -V
#$ -pe smp 10
#$ -l h_vmem=50G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile

dir="${1%.bam}"

cufflinks -p 10 -q --max-bundle-frags 30000000 --GTF-guide /opt/databases/genomes/Ensembl/Homo_sapiens.GRCh37.Ensembl.gtf \
		  -o ./Cufflinks/$dir ./Tophat/accepted_hits/$1

mv ./Cufflinks/$dir/transcripts.gtf ./Cufflinks/transcripts/$dir'.gtf'
