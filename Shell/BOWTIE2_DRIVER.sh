#!/bin/bash
for i in {1..8}
do
    qsub ./BOWTIE2_RUNNER.sh $i 
done