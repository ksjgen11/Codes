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

# merge Inputs
nice makeTagDirectory $chippath/1/IN_a -d $chippath/1/IN_1/ $chippath/1/IN_2/ -tbp 1 -unique -genome mm10
nice makeTagDirectory $chippath/1/IN_b -d $chippath/1/IN_3/ $chippath/1/IN_4/ -tbp 1 -unique -genome mm10
nice makeTagDirectory $chippath/1/IN_c -d $chippath/1/IN_5/ $chippath/1/IN_6/ -tbp 1 -unique -genome mm10

# findPeaks
nice findPeaks $chippath/1/VSF_7/ -style factor -o $chippath/1/VSF_7_input.txt -i $chippath/1/IN_a -tbp 1
nice findPeaks $chippath/1/VSF_8/ -style factor -o $chippath/1/VSF_8_input.txt -i $chippath/1/IN_a -tbp 1
nice findPeaks $chippath/1/VSF_9/ -style factor -o $chippath/1/VSF_9_input.txt -i $chippath/1/IN_b -tbp 1
nice findPeaks $chippath/1/VSF_10/ -style factor -o $chippath/1/VSF_10_input.txt -i $chippath/1/IN_b -tbp 1
nice findPeaks $chippath/1/VSF_11/ -style factor -o $chippath/1/VSF_11_input.txt -i $chippath/1/IN_c -tbp 1
nice findPeaks $chippath/1/VSF_12/ -style factor -o $chippath/1/VSF_12_input.txt -i $chippath/1/IN_c -tbp 1
nice findPeaks $chippath/1/TF_13/ -style factor -o $chippath/1/TF_13_input.txt -i $chippath/1/IN_a -tbp 1
nice findPeaks $chippath/1/TF_14/ -style factor -o $chippath/1/TF_14_input.txt -i $chippath/1/IN_a -tbp 1
nice findPeaks $chippath/1/TF_15/ -style factor -o $chippath/1/TF_15_input.txt -i $chippath/1/IN_b -tbp 1
nice findPeaks $chippath/1/TF_16/ -style factor -o $chippath/1/TF_16_input.txt -i $chippath/1/IN_b -tbp 1
nice findPeaks $chippath/1/TF_17/ -style factor -o $chippath/1/TF_17_input.txt -i $chippath/1/IN_c -tbp 1
nice findPeaks $chippath/1/TF_18/ -style factor -o $chippath/1/TF_18_input.txt -i $chippath/1/IN_c -tbp 1
