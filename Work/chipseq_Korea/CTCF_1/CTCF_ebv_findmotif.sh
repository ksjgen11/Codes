#set the directory paths
indexpath=/NextGenSeqData/project-data/sejun/chipseq/2/index/hg19_ebv/
chippath=/NextGenSeqData/project-data/sejun/chipseq/2
datapath=/NextGenSeqData/project-data/sejun/chipseq/2/2/trimmed
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/2
linkpath=/webdata/bric-data/sejun/chipseq/2

# find motifs
nice findMotifsGenome.pl $chippath/ctcf_ebv_input.txt $indexpath/hg19_ebv.fa $chippath/ctcf_ebv_motifs/ \
 -size 200 -mask -p 4

