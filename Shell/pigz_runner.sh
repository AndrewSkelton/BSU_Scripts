#!/bin/bash
#$ -cwd -V
#$ -pe smp 5
#$ -l h_vmem=8G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile
pigz -d -p 4 $1