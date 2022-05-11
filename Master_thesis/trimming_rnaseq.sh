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
