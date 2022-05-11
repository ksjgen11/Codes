#set the directory paths
chippath=/NextGenSeqData/project-data/sejun/chipseq/
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/
linkpath=/webdata/bric-data/sejun/chipseq/

# remove outofbounds
nice removeOutOfBoundsReads.pl $chippath/1/IN_1/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/IN_2/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/IN_3/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/IN_4/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/IN_5/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/IN_6/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/VSF_7/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/VSF_8/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/VSF_9/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/VSF_10/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/VSF_11/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/VSF_12/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/TF_13/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/TF_14/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/TF_15/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/TF_16/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/TF_17/ mm10
nice removeOutOfBoundsReads.pl $chippath/1/TF_18/ mm10

# make UCSC files
nice makeUCSCfile $chippath/1/IN_1/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/IN_2/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/IN_3/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/IN_4/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/IN_5/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/IN_6/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/VSF_7/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/VSF_8/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/VSF_9/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/VSF_10/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/VSF_11/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/VSF_12/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/TF_13/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/TF_14/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/TF_15/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/TF_16/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/TF_17/ -tbp 1 -style chipseq
nice makeUCSCfile $chippath/1/TF_18/ -tbp 1 -style chipseq


# make BigWig files
nice makeBigWig.pl $chippath/1/IN_1/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/IN_2/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/IN_3/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/IN_4/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/IN_5/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/IN_6/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/VSF_7/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/VSF_8/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/VSF_9/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/VSF_10/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/VSF_11/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/VSF_12/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/TF_13/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/TF_14/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/TF_15/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/TF_16/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/TF_17/ mm10 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $chippath/1/TF_18/ mm10 -normal -force -url $httppath -webdir $linkpath
