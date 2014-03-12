#!/bin/bash
declare -a arr=("A" "B")

mkdir accepted_hits

for i in "${arr[@]}"
do
        for j in {1..8}
        do
                i1='Flowcell_'$i'_'$j
                #cp $i1/accepted_hits.bam ./accepted_hits/$i1'.bam'
                cp $i1/accepted_hits.bam.bai ./accepted_hits/$i1'.bam.bai'
        done
done