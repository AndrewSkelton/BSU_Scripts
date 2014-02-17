#!/bin/bash
#$ -cwd -V
#$ -pe smp 5
#$ -l h_vmem=10G
source ~/.bash_profile
fastx_trimmer -t 20 -z -i $1 | pigz -p 4 $TMPDIR/$1.gz
mv $TMPDIR/$1.gz ~/WORKING_DATA/