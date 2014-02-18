#!/bin/bash

for current in *.sam
do
    qsub ./SAM_TO_BAM_RUNNER.sh $current
done


