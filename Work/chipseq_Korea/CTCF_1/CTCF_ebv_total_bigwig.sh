#set the directory paths
indexpath=/NextGenSeqData/project-data/sejun/chipseq/2/index/hg19_ebv/
chippath=/NextGenSeqData/project-data/sejun/chipseq/2
datapath=/NextGenSeqData/project-data/sejun/chipseq/2/2/trimmed
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/2
linkpath=/webdata/bric-data/sejun/chipseq/2

# make UCSC files
nice makeUCSCfile $chippath/CTCF_ebv/ -tbp 1 -style chipseq

nice makeUCSCfile $chippath/Input_ebv/ -tbp 1 -style chipseq

# make bigwig files
nice makeBigWig.pl $chippath/CTCF_ebv/ -normal -force -url $httppath -webdir $linkpath

nice makeBigWig.pl $chippath/Input_ebv/ -normal -force -url $httppath -webdir $linkpath

