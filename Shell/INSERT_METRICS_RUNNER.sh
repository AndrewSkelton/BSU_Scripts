#!/bin/bash
#$ -cwd -V
#$ -pe smp 6
#$ -l h_vmem=5G
#$ -e ~/log/
#$ -o ~/log
source ~/.bash_profile
filename=$(basename "$1")
filename="${filename%.*}"
output=$filename'.metrics'
hist=$filename'.hist'
java -Xms5g -jar /opt/software/bsu/bin/CollectInsertSizeMetrics.jar INPUT=./$1 OUTPUT=metrics/$output HISTOGRAM_FILE=metrics/$hist
#mean=`cut metrics/$output -f 5 | head`
#std_dev=`cut metrics/$output -f 6 | head`
#str_out=$filename" "$mean" "$std_dev
#echo $str_out >> insert_metrics.log
