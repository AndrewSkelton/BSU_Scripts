#!/bin/bash
for files in *.html
do
 mv "$files" "${files%.txt}.fastq"
done
