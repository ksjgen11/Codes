rnapath=/NextGenSeqData/project-data/sejun/rnaseq/
trans_index=$rnapath/index/
genome_index=/k/genomes/hg38/index/bowtie2_full/hg38

bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD001_ACAGTG_L003_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD001.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD002_GCCAAT_L003_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD002.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD003_CAGATC_L004_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD003.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD004_ACTTGA_L004_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD004.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD005_GATCAG_L004_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD005.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD006_TAGCTT_L004_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD006.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD008_CTTGTA_L004_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD008.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD009_ACTTGA_L007_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD009.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD010_GATCAG_L007_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD010.sam --no-unal --local -t
bowtie2 -p 8 -x $trans_index/hg38.refGene \
 -U $rnapath/trimmed/LD011_TAGCTT_L007_R1_001_trimmed.fa \
 -S $rnapath/mapping/LD011.sam --no-unal --local -t
