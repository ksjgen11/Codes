
programpath=~/downloads/Trimmomatic-0.39
filepath=/mnt/d/work/ChIPseq/
trim=$filepath/trimmed/
untrim=$filepath/untrimmed/

java -jar $programpath/trimmomatic-0.39.jar PE -threads 8 -phred33 $filepath/CTCF-ab_1.fastq.gz $filepath/CTCF-ab_2.fastq.gz $trim/CTCF-ab_1_trimmed.fa $untrim/CTCF-ab_1_untrimmed.fa \
 $trim/CTCF-ab_2_trimmed.fa $untrim/CTCF-ab_2_untrimmed.fa HEADCROP:10 ILLUMINACLIP:$filepath/references/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
java -jar $programpath/trimmomatic-0.39.jar PE -threads 8 -phred33 $filepath/Ig-G_1.fastq.gz $filepath/Ig-G_2.fastq.gz $trim/Ig-G_1_trimmed.fa $untrim/Ig-G_1_untrimmed.fa \
 $trim/Ig-G_2_trimmed.fa $untrim/Ig-G_2_untrimmed.fa HEADCROP:10 ILLUMINACLIP:$filepath/references/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
java -jar $programpath/trimmomatic-0.39.jar PE -threads 8 -phred33 $filepath/Input_1.fastq.gz $filepath/Input_2.fastq.gz $trim/Input_1_trimmed.fa $untrim/Input_1_untrimmed.fa \
 $trim/Input_2_trimmed.fa $untrim/Input_2_untrimmed.fa HEADCROP:10 ILLUMINACLIP:$filepath/references/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
