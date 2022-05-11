#!/bin/bash
darkred='\e[0;31m'
white='\e[1;37m'
lightblue='\e[1;34m'
darkgreen='\e[0;32m'
pink='\e[1;35m'

indexpath=/k/genomes/hg38/index/STAR_full
rnapath=/NextGenSeqData/project-data/sejun/rnaseq/star_mapping
datapath=/NextGenSeqData/project-data/sejun/rnaseq/trimmed


STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD001_trimmed.fa \
 --outFileNamePrefix $rnapath/LD001 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD002_trimmed.fa \
 --outFileNamePrefix $rnapath/LD002 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD005_trimmed.fa \
 --outFileNamePrefix $rnapath/LD005 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD006_trimmed.fa \
 --outFileNamePrefix $rnapath/LD006 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD007_trimmed.fa \
 --outFileNamePrefix $rnapath/LD007 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD008_trimmed.fa \
 --outFileNamePrefix $rnapath/LD008 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD009_trimmed.fa \
 --outFileNamePrefix $rnapath/LD009 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD010_trimmed.fa \
 --outFileNamePrefix $rnapath/LD010 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD011_trimmed.fa \
 --outFileNamePrefix $rnapath/LD011 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

