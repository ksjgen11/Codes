#!/bin/bash

# set the file and program paths
rnapath=/NextGenSeqData/project-data/sejun/rnaseq_arnes/
datapath=$rnapath/files/
trimmomaticpath=~/data/Trimmomatic-0.39
trim=$datapath/trimmed/
untrim=$datapath/untrimmed/
adapter=/k/genomes/adapters/fa/adapters.fa

# do trimming after quality check
# I don't like hardcording like this....but I couldn't find the way.
java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL001_CGATGT_L005_R1_001.fastq.gz $datapath/KL001_CGATGT_L005_R2_001.fastq.gz \
$trim/KL001_R1.fa $untrim/KL001_R1_untrimmed.fastq.gz \
$trim/KL001_R2.fa $untrim/KL001_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL002_TGACCA_L005_R1_001.fastq.gz $datapath/KL002_TGACCA_L005_R2_001.fastq.gz \
$trim/KL002_R1.fa $untrim/KL002_R1_untrimmed.fastq.gz \
$trim/KL002_R2.fa $untrim/KL002_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# only KL003 is single-end file
java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 16 -phred33 \
$datapath/KL003_ACAGTG_L005_R2_001.fastq.gz \
$trim/KL003_R2.fa \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL004_GCCAAT_L006_R1_001.fastq.gz $datapath/KL004_GCCAAT_L006_R2_001.fastq.gz \
$trim/KL004_R1.fa $untrim/KL004_R1_untrimmed.fastq.gz \
$trim/KL004_R2.fa $untrim/KL004_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL005_CAGATC_L006_R1_001.fastq.gz $datapath/KL005_CAGATC_L006_R2_001.fastq.gz \
$trim/KL005_R1.fa $untrim/KL005_R1_untrimmed.fastq.gz \
$trim/KL005_R2.fa $untrim/KL005_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL006_CTTGTA_L006_R1_001.fastq.gz $datapath/KL006_CTTGTA_L006_R2_001.fastq.gz \
$trim/KL006_R1.fa $untrim/KL006_R1_untrimmed.fastq.gz \
$trim/KL006_R2.fa $untrim/KL006_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL007_AGTCAA_L007_R1_001.fastq.gz $datapath/KL007_AGTCAA_L007_R2_001.fastq.gz \
$trim/KL007_R1.fa $untrim/KL007_R1_untrimmed.fastq.gz \
$trim/KL007_R2.fa $untrim/KL007_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL008_AGTTCC_L007_R1_001.fastq.gz $datapath/KL008_AGTTCC_L007_R2_001.fastq.gz \
$trim/KL008_R1.fa $untrim/KL008_R1_untrimmed.fastq.gz \
$trim/KL008_R2.fa $untrim/KL008_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL009_ATGTCA_L007_R1_001.fastq.gz $datapath/KL009_ATGTCA_L007_R2_001.fastq.gz \
$trim/KL009_R1.fa $untrim/KL009_R1_untrimmed.fastq.gz \
$trim/KL009_R2.fa $untrim/KL009_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL010_CCGTCC_L008_R1_001.fastq.gz $datapath/KL010_CCGTCC_L008_R2_001.fastq.gz \
$trim/KL010_R1.fa $untrim/KL010_R1_untrimmed.fastq.gz \
$trim/KL010_R2.fa $untrim/KL010_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL011_GTCCGC_L008_R1_001.fastq.gz $datapath/KL011_GTCCGC_L008_R2_001.fastq.gz \
$trim/KL011_R1.fa $untrim/KL011_R1_untrimmed.fastq.gz \
$trim/KL011_R2.fa $untrim/KL011_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL012_GTGAAA_L008_R1_001.fastq.gz $datapath/KL012_GTGAAA_L008_R2_001.fastq.gz \
$trim/KL012_R1.fa $untrim/KL012_R1_untrimmed.fastq.gz \
$trim/KL012_R2.fa $untrim/KL012_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL013_CGATGT_L007_R1_001.fastq.gz $datapath/KL013_CGATGT_L007_R2_001.fastq.gz \
$trim/KL013_R1.fa $untrim/KL013_R1_untrimmed.fastq.gz \
$trim/KL013_R2.fa $untrim/KL013_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL014_TGACCA_L007_R1_001.fastq.gz $datapath/KL014_TGACCA_L007_R2_001.fastq.gz \
$trim/KL014_R1.fa $untrim/KL014_R1_untrimmed.fastq.gz \
$trim/KL014_R2.fa $untrim/KL014_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL015_ACAGTG_L007_R1_001.fastq.gz $datapath/KL015_ACAGTG_L007_R2_001.fastq.gz \
$trim/KL015_R1.fa $untrim/KL015_R1_untrimmed.fastq.gz \
$trim/KL015_R2.fa $untrim/KL015_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL016_GCCAAT_L008_R1_001.fastq.gz $datapath/KL016_GCCAAT_L008_R2_001.fastq.gz \
$trim/KL016_R1.fa $untrim/KL016_R1_untrimmed.fastq.gz \
$trim/KL016_R2.fa $untrim/KL016_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL017_CAGATC_L008_R1_001.fastq.gz $datapath/KL017_CAGATC_L008_R2_001.fastq.gz \
$trim/KL017_R1.fa $untrim/KL017_R1_untrimmed.fastq.gz \
$trim/KL017_R2.fa $untrim/KL017_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL018_CTTGTA_L008_R1_001.fastq.gz $datapath/KL018_CTTGTA_L008_R2_001.fastq.gz \
$trim/KL018_R1.fa $untrim/KL018_R1_untrimmed.fastq.gz \
$trim/KL018_R2.fa $untrim/KL018_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
$datapath/KL019_AGTCAA_L006_R1_001.fastq.gz $datapath/KL019_AGTCAA_L006_R2_001.fastq.gz \
$trim/KL019_R1.fa $untrim/KL019_R1_untrimmed.fastq.gz \
$trim/KL019_R2.fa $untrim/KL019_R2_untrimmed.fastq.gz \
CROP:95 HEADCROP:15 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

