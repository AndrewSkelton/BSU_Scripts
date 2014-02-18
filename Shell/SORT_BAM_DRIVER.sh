#!/bin/bash

for current in *.bam
do
    qsub ./SORT_BAM_RUNNER.sh $current
done


