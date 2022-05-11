trimmomaticpath=~/data/Trimmomatic-0.39
datapath=/NextGenSeqData/project-data/sejun/rnaseq/
trim=/NextGenSeqData/project-data/sejun/rnaseq/trimmed/
untrim=/NextGenSeqData/project-data/sejun/rnaseq/untrimmed/
adapter=/k/genomes/adapters/fa/adapters.fa

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD001_ACAGTG_L003_R1_001.fastq.gz $trim/LD001_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD002_GCCAAT_L003_R1_001.fastq.gz $trim/LD002_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD003_CAGATC_L004_R1_001.fastq.gz $trim/LD003_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD004_ACTTGA_L004_R1_001.fastq.gz $trim/LD004_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD005_GATCAG_L004_R1_001.fastq.gz $trim/LD005_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD006_TAGCTT_L004_R1_001.fastq.gz $trim/LD006_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD007_GGCTAC_L004_R1_001.fastq.gz $trim/LD007_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD008_CTTGTA_L004_R1_001.fastq.gz $trim/LD008_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD009_ACTTGA_L007_R1_001.fastq.gz $trim/LD009_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD010_GATCAG_L007_R1_001.fastq.gz $trim/LD010_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trimmomaticpath/trimmomatic-0.39.jar SE -threads 4 -phred33 \
$datapath/LD011_TAGCTT_L007_R1_001.fastq.gz $trim/LD011_trimmed.fa \
 CROP:100 HEADCROP:13 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# mapping on the genome


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


# convert sam file to bam
samtools view -Sb $rnapath/LD001Aligned.out.sam > $rnapath/LD001.bam
samtools view -Sb $rnapath/LD002Aligned.out.sam > $rnapath/LD002.bam
samtools view -Sb $rnapath/LD005Aligned.out.sam > $rnapath/LD005.bam
samtools view -Sb $rnapath/LD006Aligned.out.sam > $rnapath/LD006.bam
samtools view -Sb $rnapath/LD007Aligned.out.sam > $rnapath/LD007.bam
samtools view -Sb $rnapath/LD008Aligned.out.sam > $rnapath/LD008.bam
samtools view -Sb $rnapath/LD009Aligned.out.sam > $rnapath/LD009.bam
samtools view -Sb $rnapath/LD010Aligned.out.sam > $rnapath/LD010.bam
samtools view -Sb $rnapath/LD011Aligned.out.sam > $rnapath/LD011.bam
