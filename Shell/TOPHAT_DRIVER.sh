#!/bin/bash
declare -a arr=("A" "B")
declare -a mate_inner_dist=(29.40 35.14 27.94 27.93 39.71 38.59 50.61 43.47 36.19 41.65 36.66 40.56 49.01 27.27 52.62 53.53)
declare -a std_dev=(32.71 35.31 38.49 41.83 38.62 35.75 33.57 33.32 36.00 41.90 38.77 60.84 53.84 40.54 49.83 47.89)

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

