indexpath=/NextGenSeqData/project-data/sejun/chipseq/2/index/ebv
chippath=/NextGenSeqData/project-data/sejun/chipseq/2
datapath=/NextGenSeqData/project-data/sejun/chipseq/2/2/trimmed

bowtie2 -p 8 -x $indexpath/ebv \
 -1 $datapath/CTCF-ab_1_trimmed.fa -2 $datapath/CTCF-ab_2_trimmed.fa \
 -S $chippath/CTCF_ebv_only.sam --no-unal --local --phred33 -t
