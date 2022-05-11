#!/bin/bash
darkred='\e[0;31m'
white='\e[1;37m'
lightblue='\e[1;34m'
darkgreen='\e[0;32m'
pink='\e[1;35m'
##################################################################################################
#INSTRUCTIONS
echo -e "${darkred}INSTRUCTIONS OF THIS SCRIPT :"
echo -e "Reference genome is default: GRCh38"
echo -e "Decompressing of fastq.gz file is not required."
echo -e "The raw data should be paired-end data."
echo -e "The adapter sequences used for library preparation is required: TruSeq3-PE-2.fa."
echo -e "The STAR index of reference genome is required."
echo -e "${darkgreen}To interrupt this script, type Ctrl+C"
######################################################################################################
####################################################################################################
#Files with the reads
fastqR1=()
fastqR2=()
fastqcutR1=()
fastqcutR2=()
echo -e "${darkred}Please inform the the absolute path and the name of raw fastq file. Since ATAC-seq is in paired-end, you need to check there are two fastq files called R1 and R2.${lightblue}"
    echo -e "${darkred}Indicate the information of R1 file ${pink}"
    read fastqR1
    echo -e "${darkred}Is your file name ${lightblue} $fastqR1 ${darkred}? If not, press Ctrl+C${white}"
    fastqR1=$fastqR1
    fastqcutR1=`awk -F'/' '{print $(NF)}' <<< $fastqR1`
    echo -e "${darkred}Indicate the information of R2 file ${pink}"
    read fastqR2
    echo -e "${darkred}Is your file name ${lightblue} $fastqR2 ${darkred}? If not, press Ctrl+C${white}"
    fastqR2=$fastqR2
    fastqcutR2=`awk -F'/' '{print $(NF)}' <<< $fastqR2`
########################################################################################################
###################################################################################################
# Define sample name used for the file name
echo -e "${darkred}Enter the name of sample${pink}"
read sample
####################################################################################################
###################################################################################################
# File with the list of all adapters used in the libraries preparation
adapter=/home/audrey/reference_genomes/adapters/TruSeq3-PE-2.fa
####################################################################################################
gunzip *
# Remove adapter sequences
echo -e "${darkgreen}Remove adapter sequences :${lightblue}"
trim1=`echo "$fastqR1.trim.fq"`
trim2=`echo "$fastqR2.trim.fq"`
untrim1=`echo "$fastqR1.untrim.fq"`
untrim2=`echo "$fastqR2.untrim.fq"`
TrimmomaticPE -threads 4 -phred33 $fastqR1 $fastqR2 $trim1 $untrim1 $trim2 $untrim2 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:8
#mapping reads to reference genome
STAR --genomeDir /home/audrey/reference_genomes/hg38 --readFilesIn $trim1 $trim2 --runThreadN 4 --outFileNamePrefix $sample. --outFilterIntronMotifs RemoveNoncanonicalUnannotated
#convert sam to bam
sambamba view -f bam -S -o $sample.bam $sample.Aligned.out.sam
#sort by chromosome
sambamba sort -m 3GB -p -t 4 -o sortc_$sample.bam $sample.bam
#remove duplicated reads
picard.jar MarkDuplicates I=sortc_$sample.bam O=rmdup_$sample.bam M=$sample.rmdup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 
#sort by name
sambamba sort -n -m 3GB -p -t 4 -o sortn_rmdup_$sample.bam rmdup_$sample.bam
#Counting of reads
htseq-count -f bam -r name -s no --nonunique all sortn_rmdup_$sample.bam /home/audrey/reference_genomes/features/Homo_sapiens.GRCh38.87.gtf > $sample.count.txt

