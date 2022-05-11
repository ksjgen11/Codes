#set the directory paths
rnapath=/NextGenSeqData/project-data/sejun/rnaseq/
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/rnaseq/
linkpath=/webdata/bric-data/sejun/rnaseq/


# remove outofbounds
nice removeOutOfBoundsReads.pl $rnapath/LD001_mapping/ hg38
nice removeOutOfBoundsReads.pl $rnapath/LD002_mapping/ hg38
nice removeOutOfBoundsReads.pl $rnapath/LD005_mapping/ hg38
nice removeOutOfBoundsReads.pl $rnapath/LD006_mapping/ hg38
nice removeOutOfBoundsReads.pl $rnapath/LD007_mapping/ hg38
nice removeOutOfBoundsReads.pl $rnapath/LD008_mapping/ hg38
nice removeOutOfBoundsReads.pl $rnapath/LD009_mapping/ hg38
nice removeOutOfBoundsReads.pl $rnapath/LD010_mapping/ hg38
nice removeOutOfBoundsReads.pl $rnapath/LD011_mapping/ hg38

# make UCSC files
nice makeUCSCfile $rnapath/LD001_mapping/ -tbp 1 -fragLength given -style chipseq
nice makeUCSCfile $rnapath/LD002_mapping/ -tbp 1 -fragLength given -style chipseq
nice makeUCSCfile $rnapath/LD005_mapping/ -tbp 1 -fragLength given -style chipseq
nice makeUCSCfile $rnapath/LD006_mapping/ -tbp 1 -fragLength given -style chipseq
nice makeUCSCfile $rnapath/LD007_mapping/ -tbp 1 -fragLength given -style chipseq
nice makeUCSCfile $rnapath/LD008_mapping/ -tbp 1 -fragLength given -style chipseq
nice makeUCSCfile $rnapath/LD009_mapping/ -tbp 1 -fragLength given -style chipseq
nice makeUCSCfile $rnapath/LD010_mapping/ -tbp 1 -fragLength given -style chipseq
nice makeUCSCfile $rnapath/LD011_mapping/ -tbp 1 -fragLength given -style chipseq


# make BigWig files
nice makeBigWig.pl $rnapath/LD001_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD002_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD005_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD006_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD007_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD008_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD009_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD010_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD011_mapping/ hg38 -normal -force -url $httppath -webdir $linkpath
