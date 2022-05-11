#set the directory paths
chippath=/NextGenSeqData/project-data/sejun/chipseq/
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/
linkpath=/webdata/bric-data/sejun/chipseq/

# make UCSC files
nice makeUCSCfile $chippath/1/VSF_a/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/VSF_b/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/VSF_c/ -tbp 1 -style chipseq

# make BigWig files
nice makeBigWig.pl $chippath/1/VSF_a/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/VSF_b/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/VSF_c/ mm10 -normal -force -url $httppath -webdir $linkpath
