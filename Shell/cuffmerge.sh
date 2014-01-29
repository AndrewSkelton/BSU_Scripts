#!/bin/bash
#screen -dmS Merrrrrrrge bash -c 'cuffmerge -s /data/genomes/grch37/grch37.fa assemblies.txt; exec bash' 

screen -dmS Merrrrrrrge1 bash -c 'cuffmerge -s /data/genomes/grch37/grch37.fa -g /data/genomes/grch37/genes.gtf -o ./Merge1 assemblies.txt; exec bash' 
