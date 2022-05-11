indexpath=/k/genomes/hg38/index/STAR_full
rnapath=/NextGenSeqData/project-data/sejun/rnaseq/star_mapping
datapath=/NextGenSeqData/project-data/sejun/rnaseq/trimmed


STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD001_trimmed.fa \
 --outFileNamePrefix $rnapath/LD001

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD002_trimmed.fa \
 --outFileNamePrefix $rnapath/LD002

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD003_trimmed.fa \
 --outFileNamePrefix $rnapath/LD003

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD004_trimmed.fa \
 --outFileNamePrefix $rnapath/LD004

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD005_trimmed.fa \
 --outFileNamePrefix $rnapath/LD005

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD006_trimmed.fa \
 --outFileNamePrefix $rnapath/LD006

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD007_trimmed.fa \
 --outFileNamePrefix $rnapath/LD007

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD008_trimmed.fa \
 --outFileNamePrefix $rnapath/LD008

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD009_trimmed.fa \
 --outFileNamePrefix $rnapath/LD009

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD010_trimmed.fa \
 --outFileNamePrefix $rnapath/LD010

STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn $datapath/LD011_trimmed.fa \
 --outFileNamePrefix $rnapath/LD011
