#!/bin/bash
echo "Start : " `date` > speed.log
cache=$((wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz) 2>&1 | tail -2 | head -1 | awk '{print $3 $4 }')
echo "Finish : " `date` >> speed.log
echo "Download speed : $cache " >> speed.log
