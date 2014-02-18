#!/bin/bash
#$ -cwd -V
#$ -pe smp 6
#$ -l h_vmem=10G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile
filename=$(basename "$1")
filename="${filename%.*}"
output=$filename'.sorted'
samtools sort $1 ./sorted/$output