#!/bin/bash

for current in *.bam
do
    qsub ./INSERT_METRICS_RUNNER.sh $current
done

