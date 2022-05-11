indexpath=/k/genomes/hg38/index/bowtie2_full/
rnapath=/NextGenSeqData/project-data/sejun/rnaseq
datapath=/NextGenSeqData/project-data/sejun/rnaseq/trimmed


bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD001_trimmed.fa \
 -S $rnapath/mapping/LD001.sam --no-unal --phred33 --local -t
bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD002_trimmed.fa \
 -S $rnapath/mapping/LD002.sam --no-unal --phred33 --local -t
bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD005_trimmed.fa \
 -S $rnapath/mapping/LD005.sam --no-unal --phred33 --local -t
bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD006_trimmed.fa \
 -S $rnapath/mapping/LD006.sam --no-unal --phred33 --local -t
bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD007_trimmed.fa \
 -S $rnapath/mapping/LD007.sam --no-unal --phred33 --local -t
bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD008_trimmed.fa \
 -S $rnapath/mapping/LD008.sam --no-unal --phred33 --local -t
bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD009_trimmed.fa \
 -S $rnapath/mapping/LD009.sam --no-unal --phred33 --local -t
bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD010_trimmed.fa \
 -S $rnapath/mapping/LD010.sam --no-unal --phred33 --local -t
bowtie2 -p 8 -x $indexpath/hg38 \
 -U $rnapath/trimmed/LD011_trimmed.fa \
 -S $rnapath/mapping/LD011.sam --no-unal --phred33 --local -t
