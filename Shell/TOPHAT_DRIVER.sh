#!/bin/bash
declare -a arr=("A" "B")
declare -a mate_inner_dist=(29 35 28 28 40 39 51 43 36 42 37 41 49 27 53 54 )
declare -a std_dev=(33 35 38 42 39 36 34 33 36 42 39 61 54 41 50 48 )

pointer=0

#echo "${mate_inner_dist[$pointer]}"
#echo "${std_dev[$pointer]}"

for i in "${arr[@]}"
do
        for j in {1..8}
        do
                i1='Flowcell_'$i'_'$j'_Forward.fastq.gz'
                i2='Flowcell_'$i'_'$j'_Reverse.fastq.gz'
                dir='Flowcell_'$i'_'$j 
				
                qsub ./TOPHAT_RUNNER.sh $i1 $i2 ${mate_inner_dist[$pointer]} ${std_dev[$pointer]} $dir

                pointer=$pointer+1
        done
done

