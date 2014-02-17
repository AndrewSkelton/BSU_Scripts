#!/bin/bash
source ~/.bash_profile
##'USE FASTX_TRIMMER FOR EACH SET OF READS
##'HARD 20 TRIM
##'COMPRESS OUTPUT TO $TMPDIR
##'SET EACH COMMAND AS A CHILD PROCESS
fastx_trimmer -t 20 -z -i Flowcell_A_1_Forward.fastq -o $TMPDIR/Flowcell_A_1_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_2_Forward.fastq -o $TMPDIR/Flowcell_A_2_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_3_Forward.fastq -o $TMPDIR/Flowcell_A_3_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_4_Forward.fastq -o $TMPDIR/Flowcell_A_4_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_5_Forward.fastq -o $TMPDIR/Flowcell_A_5_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_6_Forward.fastq -o $TMPDIR/Flowcell_A_6_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_7_Forward.fastq -o $TMPDIR/Flowcell_A_7_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_8_Forward.fastq -o $TMPDIR/Flowcell_A_8_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_1_Reverse.fastq -o $TMPDIR/Flowcell_A_1_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_2_Reverse.fastq -o $TMPDIR/Flowcell_A_2_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_3_Reverse.fastq -o $TMPDIR/Flowcell_A_3_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_4_Reverse.fastq -o $TMPDIR/Flowcell_A_4_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_5_Reverse.fastq -o $TMPDIR/Flowcell_A_5_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_6_Reverse.fastq -o $TMPDIR/Flowcell_A_6_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_7_Reverse.fastq -o $TMPDIR/Flowcell_A_7_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_A_8_Reverse.fastq -o $TMPDIR/Flowcell_A_8_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_1_Forward.fastq -o $TMPDIR/Flowcell_B_1_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_2_Forward.fastq -o $TMPDIR/Flowcell_B_2_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_3_Forward.fastq -o $TMPDIR/Flowcell_B_3_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_4_Forward.fastq -o $TMPDIR/Flowcell_B_4_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_5_Forward.fastq -o $TMPDIR/Flowcell_B_5_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_6_Forward.fastq -o $TMPDIR/Flowcell_B_6_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_7_Forward.fastq -o $TMPDIR/Flowcell_B_7_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_8_Forward.fastq -o $TMPDIR/Flowcell_B_8_Forward.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_1_Reverse.fastq -o $TMPDIR/Flowcell_B_1_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_2_Reverse.fastq -o $TMPDIR/Flowcell_B_2_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_3_Reverse.fastq -o $TMPDIR/Flowcell_B_3_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_4_Reverse.fastq -o $TMPDIR/Flowcell_B_4_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_5_Reverse.fastq -o $TMPDIR/Flowcell_B_5_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_6_Reverse.fastq -o $TMPDIR/Flowcell_B_6_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_7_Reverse.fastq -o $TMPDIR/Flowcell_B_7_Reverse.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell_B_8_Reverse.fastq -o $TMPDIR/Flowcell_B_8_Reverse.fastq.gz &
##'WAIT FOR CHILD PROCESSES TO RETURN
wait
##'MOVE ALL OUTPUTS TO WORKING_DATA DIRECTORY
mv $TMPDIR/*.fastq.gz ~/WORKING_DATA