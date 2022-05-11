#!/bin/bash

call_peak(){
    bamfile="$1"
    size=$2 # size for human = hs, mouse = mm
    
    samplename="${bamfile:33:11}"

    # Peak calling

    #echo -e "${darkgreen}Peak calling with MACS2: ${lightblue}"

    peak=`echo "Peak_calling"`
    mkdir $peak

    macs2 callpeak -t $bamfile --outdir $peak -f BAMPE -g $size -n "$samplename" \
    -q 0.01 --nomodel --shift 0

 
}