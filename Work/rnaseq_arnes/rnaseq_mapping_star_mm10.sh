#!/bin/bash

# set the file and program paths
indexpath=/k/genomes/mm10/index/STAR_full
rnapath=/NextGenSeqData/project-data/sejun/rnaseq_arnes/
datapath=$rnapath/files/
map=`echo "star_mapping"`
# mkdir $datapath/$map
star=$datapath/$map
trim=$datapath/trimmed/


# do transcriptome mapping and convert sam format to bam format in a loop

for number in 001 002 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019
do
    r1=${trim}/KL${number}_R1.fa
    r2=${trim}/KL${number}_R2.fa
    printf "mapping R1 %s\n" "$r1"
    printf "mapping R2 %s\n" "$r2"

    STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $r1 $r2 \
    --outFileNamePrefix $star/KL${number} --outFilterIntronMotifs RemoveNoncanonicalUnannotated

    samtools view --threads 16 -Sb $star/KL${number}Aligned.out.sam > $star/KL${number}.bam
    
done


# do KL003 transcriptome mapping and convert sam format to bam format because it is only single-end file
STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $trim/KL003_R2.fa \
--outFileNamePrefix $star/KL003 --outFilterIntronMotifs RemoveNoncanonicalUnannotated

samtools view --threads 16 -Sb $star/KL003Aligned.out.sam > $star/KL003.bam

# move bam file to separate with paired-end file for featurecounts
mv $star/KL003.bam $datapath/.