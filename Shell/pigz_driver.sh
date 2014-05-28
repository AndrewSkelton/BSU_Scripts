#!/bin/bash
for i in *.gz
do
	qsub ./pigz_runner.sh $i
done