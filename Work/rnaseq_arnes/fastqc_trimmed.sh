#!/bin/bash

rnapath=/NextGenSeqData/project-data/sejun/rnaseq_arnes/
datapath=$rnapath/files/
trim=$datapath/trimmed/
out=`echo "qc_trimmed"`
mkdir $datapath/$out
outpath=$datapath/$out


for file in $trim/*.fa
do
    fastqc -o $outpath $file
    
done