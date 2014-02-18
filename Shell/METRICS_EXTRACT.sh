#!/bin/bash
k="Name       Mean_Insert_Size  Standard_Deviation"
echo $k > metrics.log
for i in *.metrics
do 
	mean=`cut $i -f 5 | awk '{ if (NF > 0) { last = $NF } } END { print last }' "$@"`
	std_dev=`cut $i -f 6 | awk '{ if (NF > 0) { last = $NF } } END { print last }' "$@"`
	echo $i $mean $std_dev >> metrics.log
done
