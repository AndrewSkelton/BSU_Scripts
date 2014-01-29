#!/bin/bash
screen -dmS Cuffdiff_Genes bash -c 'cuffdiff -q -u -b /data/genomes/grch37/grch37.fa -o ./cuffdiff_out_genes /data/genomes/grch37/genes.gtf \
Flowcell_A_2/accepted_hits.bam,Flowcell_A_5/accepted_hits.bam,Flowcell_A_7/accepted_hits.bam,Flowcell_B_2/accepted_hits.bam,Flowcell_B_5/accepted_hits.bam,Flowcell_B_5/accepted_hits.bam \
Flowcell_A_1/accepted_hits.bam,Flowcell_A_3/accepted_hits.bam,Flowcell_A_4/accepted_hits.bam,Flowcell_A_6/accepted_hits.bam,Flowcell_A_8/accepted_hits.bam,\
Flowcell_B_1/accepted_hits.bam,Flowcell_B_3/accepted_hits.bam,Flowcell_B_4/accepted_hits.bam,Flowcell_B_7/accepted_hits.bam,Flowcell_B_8/accepted_hits.bam &> cuffdiff_genes.log; exec bash'
