#set the directory paths
indexpath=/NextGenSeqData/project-data/sejun/chipseq/2/index/hg19_ebv/
rnapath=/NextGenSeqData/project-data/sejun/rnaseq
datapath=/NextGenSeqData/project-data/sejun/chipseq/2/2/trimmed
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/rnaseq
linkpath=/webdata/bric-data/sejun/rnaseq

removeOutOfBoundsReads.pl $rnapath/LD001 hg38
removeOutOfBoundsReads.pl $rnapath/LD002 hg38
removeOutOfBoundsReads.pl $rnapath/LD005 hg38
removeOutOfBoundsReads.pl $rnapath/LD006 hg38
removeOutOfBoundsReads.pl $rnapath/LD007 hg38
removeOutOfBoundsReads.pl $rnapath/LD008 hg38
removeOutOfBoundsReads.pl $rnapath/LD009 hg38
removeOutOfBoundsReads.pl $rnapath/LD010 hg38
removeOutOfBoundsReads.pl $rnapath/LD011 hg38


# make UCSC files
nice makeUCSCfile $rnapath/LD001 -tbp 1 -style rnaseq -o auto
nice makeUCSCfile $rnapath/LD002 -tbp 1 -style rnaseq -o auto
nice makeUCSCfile $rnapath/LD005 -tbp 1 -style rnaseq -o auto
nice makeUCSCfile $rnapath/LD006 -tbp 1 -style rnaseq -o auto
nice makeUCSCfile $rnapath/LD007 -tbp 1 -style rnaseq -o auto
nice makeUCSCfile $rnapath/LD008 -tbp 1 -style rnaseq -o auto
nice makeUCSCfile $rnapath/LD009 -tbp 1 -style rnaseq -o auto
nice makeUCSCfile $rnapath/LD010 -tbp 1 -style rnaseq -o auto
nice makeUCSCfile $rnapath/LD011 -tbp 1 -style rnaseq -o auto

nice makeBigWig.pl $rnapath/LD001 hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD002 hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD005 hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD006 hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD007 hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD008 hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD009 hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD010 hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD011 hg38 -strand -force -url $httppath -webdir $linkpath
