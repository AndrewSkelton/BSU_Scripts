#!/bin/bash
#$ -cwd -V
#$ -pe smp 6
#$ -l h_vmem=10G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile
INPUTA='Flowcell_A_'$1'_Forward.fastq.gz'
INPUTB='Flowcell_A_'$1'_Reverse.fastq.gz'
OUTPUT='Flowcell_A_'$1'.sam'

bowtie2 -p 5 -u 3000000 --phred64-quals --very-sensitive -x /opt/databases/bowtie2/hg19 -1 $INPUTA -2 $INPUTB -S $TMPDIR/$OUTPUT
mv $TMPDIR/$OUTPUT ~/WORKING_DATA/

