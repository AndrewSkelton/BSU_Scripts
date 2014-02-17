for current_fastq in *.fastq
do
    qsub ./CLUSTER_FASTX_TRIM_RUNNER.sh $current_fastq
done