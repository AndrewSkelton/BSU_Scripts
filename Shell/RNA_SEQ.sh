##'FASTQC FOR QUALITY CONTROL REPORTS
##'USE PIGZ FOR PARALLEL DECOMPRESSION

##'I USED A HARD TRIM AT 20, HOWEVER SMART TRIMMERS CAN BE USED SUCH AS FASTQ_QUALITY_TRIMMER
##' -o DESIGNATED THE OUTPUT AND -i DESIGNATES THE INPUT. -t IS THE HARD TRIM LENGTH.
fastx_trimmer -t 20 -z -i Flowcell1_1_Forward_Sequence.fastq -o Flowcell1_1_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_2_Forward_Sequence.fastq -o Flowcell1_2_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_3_Forward_Sequence.fastq -o Flowcell1_3_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_4_Forward_Sequence.fastq -o Flowcell1_4_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_5_Forward_Sequence.fastq -o Flowcell1_5_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_6_Forward_Sequence.fastq -o Flowcell1_6_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_7_Forward_Sequence.fastq -o Flowcell1_7_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_8_Forward_Sequence.fastq -o Flowcell1_8_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_1_Reverse_Sequence.fastq -o Flowcell1_1_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_2_Reverse_Sequence.fastq -o Flowcell1_2_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_3_Reverse_Sequence.fastq -o Flowcell1_3_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_4_Reverse_Sequence.fastq -o Flowcell1_4_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_5_Reverse_Sequence.fastq -o Flowcell1_5_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_6_Reverse_Sequence.fastq -o Flowcell1_6_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_7_Reverse_Sequence.fastq -o Flowcell1_7_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell1_8_Reverse_Sequence.fastq -o Flowcell1_8_Reverse_Sequence.fastq.gz &

fastx_trimmer -t 20 -z -i Flowcell2_1_Forward_Sequence.fastq -o Flowcell2_1_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_2_Forward_Sequence.fastq -o Flowcell2_2_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_3_Forward_Sequence.fastq -o Flowcell2_3_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_4_Forward_Sequence.fastq -o Flowcell2_4_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_5_Forward_Sequence.fastq -o Flowcell2_5_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_6_Forward_Sequence.fastq -o Flowcell2_6_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_7_Forward_Sequence.fastq -o Flowcell2_7_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_8_Forward_Sequence.fastq -o Flowcell2_8_Forward_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_1_Reverse_Sequence.fastq -o Flowcell2_1_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_2_Reverse_Sequence.fastq -o Flowcell2_2_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_3_Reverse_Sequence.fastq -o Flowcell2_3_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_4_Reverse_Sequence.fastq -o Flowcell2_4_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_5_Reverse_Sequence.fastq -o Flowcell2_5_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_6_Reverse_Sequence.fastq -o Flowcell2_6_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_7_Reverse_Sequence.fastq -o Flowcell2_7_Reverse_Sequence.fastq.gz &
fastx_trimmer -t 20 -z -i Flowcell2_8_Reverse_Sequence.fastq -o Flowcell2_8_Reverse_Sequence.fastq.gz &
##'END OF QC

##'FOR ANALYSIS, I USED THE TUXEDO TOOLS. 
##'STARTING WITH TOPHAT. 
##'TOPHAT IS A SPLICE JUNCTION MAPPER FOR RNA SEQ READS. IT ALIGNS RNA SEQ READS TO MAMMEL-SIZED GENOMES
##'USING THE SHORT READ ALIGNER BOWTIE. 

###############################################################################################################################
##'THE BELOW SCRIPT IS BASH (LINUX COMMAND LINE SCRIPT) AND EACH COMMAND IS ENCAPSULATED WITHIN A "SCREEN" CALL               #
##'"SCREEN" IS A PROGRAM THAT CREATES A COMMAND LINE SESSION THAT IS DETACHABLE. THE ENCLOSED COMMANDS ARE FOR TOPHAT.        #
##'																															  #
##'FLAGS: -r THE MATE INNER DISTANCE (DEF 50), SHOULD BE CALCULATED -> (READ LENGTH - (2 * SEQUENCE LENGTH)) |PAIRED END ONLY #
##'       --phred64-quals IS WHAT KIND OF QUALITY SCORES ARE USED IN THE DATA - THIS CAN BE FOUND IN THE FASTQC REPORT        #
##'       --transcriptome-index REFERENCE TRANSCRIPTOME INDEX FROM THE HUMAN GENOME ASSISTS WITH MAPPING                      #
##'       --b2-very-sensitive MAKES BOWTIE2 DO A SENSATIVE ALIGNMENT, OFTEN SLOWER                                            #
##'       -o OUTPUT DIRECTORY (DOES NOT HAVE TO EXIST PRIOR TO COMMAND EXECUTION) 											  #
##'       -G PROVIDES TOPHAT WITH A GTF FILE CONTAINING TRANSCRIPT ANNOTATION												  #
##'																															  #
##'ONCE ALL FLAGS ARE SET, PRIOR TO THE INPUT FILES, YOU MUST INPUT A PATH TO THE REFERENCE GENOME YOU'RE PROVIDING TOPHAT.   #
##'FINALLY, SET THE PATHS FOR YOUR FORWARD AND REVERSE READS (PAIRED END ONLY) SEPERATED BY A SPACE.                          #
##'"; exec bash" AT THE END OF THE COMMAND SIMPLY FORCES THE SCREEN SESSION TO REMAIN ACTIVE, EVEN IF THE COMMAND FINISHES    #
##'																															  #
##'IF YOU ARE RUNNING THIS ON THE CLUSTER, THEN YOU WILL NOT NEED THE "SCREEN" COMMAND SURROUNDING THE TOPHAT COMMAND         #
##'																															  #
##'PLEASE READ THE COMMAND LIST BEFORE PROCEDING - http://tophat.cbcb.umd.edu/manual.shtml                                    #
##'                                              - http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml                     #
###############################################################################################################################

#!/bin/bash
screen -dmS Flowcell_A_1 bash -c 'tophat -r 44 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_A_1 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_A_1_1.fastq.gz Flowcell_A_1_2.fastq.gz; exec bash' 
screen -dmS Flowcell_A_2 bash -c 'tophat -r 44 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_A_2 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_A_2_1.fastq.gz Flowcell_A_2_2.fastq.gz; exec bash' 
screen -dmS Flowcell_A_3 bash -c 'tophat -r 44 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_A_3 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_A_3_1.fastq.gz Flowcell_A_3_2.fastq.gz; exec bash' 
screen -dmS Flowcell_A_4 bash -c 'tophat -r 44 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_A_4 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_A_4_1.fastq.gz Flowcell_A_4_2.fastq.gz; exec bash' 
screen -dmS Flowcell_A_5 bash -c 'tophat -r 44 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_A_5 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_A_5_1.fastq.gz Flowcell_A_5_2.fastq.gz; exec bash' 
screen -dmS Flowcell_A_6 bash -c 'tophat -r 44 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_A_6 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_A_6_1.fastq.gz Flowcell_A_6_2.fastq.gz; exec bash' 
screen -dmS Flowcell_A_7 bash -c 'tophat -r 44 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_A_7 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_A_7_1.fastq.gz Flowcell_A_7_2.fastq.gz; exec bash' 
screen -dmS Flowcell_A_8 bash -c 'tophat -r 44 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_A_8 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_A_8_1.fastq.gz Flowcell_A_8_2.fastq.gz; exec bash' 
screen -dmS Flowcell_B_1 bash -c 'tophat -r 94 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_B_1 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_B_1_1.fastq.gz Flowcell_B_1_2.fastq.gz; exec bash' 
screen -dmS Flowcell_B_2 bash -c 'tophat -r 94 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_B_2 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_B_2_1.fastq.gz Flowcell_B_2_2.fastq.gz; exec bash' 
screen -dmS Flowcell_B_3 bash -c 'tophat -r 94 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_B_3 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_B_3_1.fastq.gz Flowcell_B_3_2.fastq.gz; exec bash' 
screen -dmS Flowcell_B_4 bash -c 'tophat -r 94 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_B_4 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_B_4_1.fastq.gz Flowcell_B_4_2.fastq.gz; exec bash' 
screen -dmS Flowcell_B_5 bash -c 'tophat -r 94 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_B_5 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_B_5_1.fastq.gz Flowcell_B_5_2.fastq.gz; exec bash' 
screen -dmS Flowcell_B_6 bash -c 'tophat -r 94 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_B_6 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_B_6_1.fastq.gz Flowcell_B_6_2.fastq.gz; exec bash' 
screen -dmS Flowcell_B_7 bash -c 'tophat -r 94 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_B_7 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_B_7_1.fastq.gz Flowcell_B_7_2.fastq.gz; exec bash' 
screen -dmS Flowcell_B_8 bash -c 'tophat -r 94 --phred64-quals --transcriptome-index /data/genomes/grch37/genes --b2-very-sensitive -o /data/customers/Young/Project/Tophat/Flowcell_B_8 -G /data/genomes/grch37/genes.gtf /data/genomes/grch37/grch37 Flowcell_B_8_1.fastq.gz Flowcell_B_8_2.fastq.gz; exec bash' 
echo 'Running...'

##'IF YOU DID USE THE SCREEN APPROACH THEN, KILLING THE SCREEN SESSIONS ONCE YOU'VE FINISHED CAN BE VERY TEDIOUS. I HAVE A SCRIPT THAT WILL
##'KILL THEM ALL, SO DO ASK IF YOU NEED A HAND.

##'LOOK AT THE EXAMPLE WORKFLOWS TO IDENTIFY EXACTLY WHAT YOU WANT TO ACHIEVE BEFORE GOING ANY FURTHER. http://cufflinks.cbcb.umd.edu/tutorial.html

###############################################################################################################################
##'THE BELOW SCRIPT IS BASH (LINUX COMMAND LINE SCRIPT) AND EACH COMMAND IS ENCAPSULATED WITHIN A "SCREEN" CALL               #
##'"SCREEN" IS A PROGRAM THAT CREATES A COMMAND LINE SESSION THAT IS DETACHABLE. THE ENCLOSED COMMANDS ARE FOR CUFFLINKS.     #
##'																															  #
##'FLAGS:  --multi-read-correct TELLS CUFFLINKS TO DO AN INITIAL ESTIMATION TO MORE ACCURATELY WEIGH READS MAPPING TO MULTIPLE#
##'								LOCATIONS IN THE GENOME	- PROVIDE REFERENCE GENOME FASTA FILE								  #
##'        --GTF-guide PROVIDES REFERENCE ANNOTATION																		  #
##'        -o OUTPUT DIRECTORY																								  #
##'        THE LAST PARAMETER SHOULD BE THE accepted_hits.bam FILE FROM TOPHAT                                                # 
##'																															  #
##'PLEASE READ THE COMMAND LIST BEFORE PROCEDING - http://cufflinks.cbcb.umd.edu/manual.html	                              #
###############################################################################################################################

#!/bin/bash
screen -dmS Flowcell_A_1 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_A_1/Cufflinks2 Flowcell_A_1/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_A_2 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_A_2/Cufflinks2 Flowcell_A_2/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_A_3 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_A_3/Cufflinks2 Flowcell_A_3/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_A_4 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_A_4/Cufflinks2 Flowcell_A_4/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_A_5 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_A_5/Cufflinks2 Flowcell_A_5/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_A_6 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_A_6/Cufflinks2 Flowcell_A_6/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_A_7 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_A_7/Cufflinks2 Flowcell_A_7/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_A_8 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_A_8/Cufflinks2 Flowcell_A_8/accepted_hits.bam; exec bash'  
screen -dmS Flowcell_B_1 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_B_1/Cufflinks2 Flowcell_B_1/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_B_2 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_B_2/Cufflinks2 Flowcell_B_2/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_B_3 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_B_3/Cufflinks2 Flowcell_B_3/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_B_4 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_B_4/Cufflinks2 Flowcell_B_4/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_B_5 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_B_5/Cufflinks2 Flowcell_B_5/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_B_6 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_B_6/Cufflinks2 Flowcell_B_6/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_B_7 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_B_7/Cufflinks2 Flowcell_B_7/accepted_hits.bam; exec bash' 
screen -dmS Flowcell_B_8 bash -c 'cufflinks --multi-read-correct /data/genomes/grch37/grch37.fa --GTF-guide /data/genomes/grch37/genes.gtf -o Flowcell_B_8/Cufflinks2 Flowcell_B_8/accepted_hits.bam; exec bash'  

##'CREATE A FILE CALLED assemblies.txt, IT SHOULD LOOK SIMILAR TO THIS 
Flowcell_A_1/Cufflinks2/transcripts.gtf
Flowcell_A_2/Cufflinks2/transcripts.gtf
Flowcell_A_3/Cufflinks2/transcripts.gtf
Flowcell_A_4/Cufflinks2/transcripts.gtf
Flowcell_A_5/Cufflinks2/transcripts.gtf
Flowcell_A_6/Cufflinks2/transcripts.gtf
Flowcell_A_8/Cufflinks2/transcripts.gtf
Flowcell_A_7/Cufflinks2/transcripts.gtf
Flowcell_B_1/Cufflinks2/transcripts.gtf
Flowcell_B_2/Cufflinks2/transcripts.gtf
Flowcell_B_3/Cufflinks2/transcripts.gtf
Flowcell_B_4/Cufflinks2/transcripts.gtf
Flowcell_B_5/Cufflinks2/transcripts.gtf
Flowcell_B_6/Cufflinks2/transcripts.gtf
Flowcell_B_8/Cufflinks2/transcripts.gtf
Flowcell_B_7/Cufflinks2/transcripts.gtf

##'CUFFMERGE MERGES SEVERAL CUFFLINKS ASSEMBLIES THROUGH THE transcripts.gtf FILE 
screen -dmS Merrrrrrrge1 bash -c 'cuffmerge -s /data/genomes/grch37/grch37.fa -g /data/genomes/grch37/genes.gtf -o ./Merge1 assemblies.txt; exec bash' 

##'FROM THIS STAGE THERE ARE TWO VERSIONS OF CUFFDIFF THAT I RAN - USING THE CUFFMERGE RESULT AND USING THE HUMAN GENOME GTF REFERENCE FILE

##'CUFFMERGE RESULTS
##' -u MULTI-READ CORRECTION, AS SEEN IN CUFFLINKS
##' -q QUIET - CREATES A READABLE OUTPUT
##' -b FRAG BIAS CORRECTION WHICH USES THE ORIGINAL FASTA FILE THAT THE READS WERE MAPPED TO, TO IMPROVE ACCURACY 
##' -o OUTPUT DIRECTORY
##' FINAL ARGUMENT BEFORE INPUT FILES SHOULD BE THE ANNOTATION FILE
##' FIRST SAMPLES OF THE SAME TYPE SHOULD BE SEPERATED BY ',' AND THE DIFFERENT SAMPLE TYPES SHOULD BE SEPERATED BY A SPACE.  

#!/bin/bash
screen -dmS Cuffdiff_Merge bash -c 'cuffdiff -q -u -b /data/genomes/grch37/grch37.fa -o ./cuffdiff_out_merge ./Merge1/merged.gtf \
Flowcell_A_2/accepted_hits.bam,Flowcell_A_5/accepted_hits.bam,Flowcell_A_7/accepted_hits.bam,Flowcell_B_2/accepted_hits.bam,Flowcell_B_5/accepted_hits.bam,Flowcell_B_5/accepted_hits.bam \
Flowcell_A_1/accepted_hits.bam,Flowcell_A_3/accepted_hits.bam,Flowcell_A_4/accepted_hits.bam,Flowcell_A_6/accepted_hits.bam,Flowcell_A_8/accepted_hits.bam,\
Flowcell_B_1/accepted_hits.bam,Flowcell_B_3/accepted_hits.bam,Flowcell_B_4/accepted_hits.bam,Flowcell_B_7/accepted_hits.bam,Flowcell_B_8/accepted_hits.bam &> cuffdiff_merge.log; exec bash'

#!/bin/bash
screen -dmS Cuffdiff_Genes bash -c 'cuffdiff -q -u -b /data/genomes/grch37/grch37.fa -o ./cuffdiff_out_genes /data/genomes/grch37/genes.gtf \
Flowcell_A_2/accepted_hits.bam,Flowcell_A_5/accepted_hits.bam,Flowcell_A_7/accepted_hits.bam,Flowcell_B_2/accepted_hits.bam,Flowcell_B_5/accepted_hits.bam,Flowcell_B_5/accepted_hits.bam \
Flowcell_A_1/accepted_hits.bam,Flowcell_A_3/accepted_hits.bam,Flowcell_A_4/accepted_hits.bam,Flowcell_A_6/accepted_hits.bam,Flowcell_A_8/accepted_hits.bam,\
Flowcell_B_1/accepted_hits.bam,Flowcell_B_3/accepted_hits.bam,Flowcell_B_4/accepted_hits.bam,Flowcell_B_7/accepted_hits.bam,Flowcell_B_8/accepted_hits.bam &> cuffdiff_genes.log; exec bash'