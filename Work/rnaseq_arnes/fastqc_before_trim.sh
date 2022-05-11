#!/bin/bash

rnapath=/NextGenSeqData/project-data/sejun/rnaseq_arnes/
datapath=$rnapath/files/
out=`echo "qc"`
mkdir $datapath/$out
outpath=$datapath/$out

# quality check before trimming

for file in $datapath/*.fastq.gz
do
    fastqc -o $outpath $file
    
done