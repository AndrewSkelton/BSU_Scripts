#!/bin/bash
declare -a arr=("A" "B")

for i in "${arr[@]}"
do
        for j in {1..8}
        do
                i1='Flowcell_'$i'_'$j
                qsub ./CUFFLINKS_RUNNER.sh $i1
				echo $i1
        done
done