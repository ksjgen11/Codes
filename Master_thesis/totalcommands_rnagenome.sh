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
nice makeUCSCfile $rnapath/LD001_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq
nice makeUCSCfile $rnapath/LD002_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq
nice makeUCSCfile $rnapath/LD005_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq
nice makeUCSCfile $rnapath/LD006_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq
nice makeUCSCfile $rnapath/LD007_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq
nice makeUCSCfile $rnapath/LD008_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq
nice makeUCSCfile $rnapath/LD009_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq
nice makeUCSCfile $rnapath/LD010_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq
nice makeUCSCfile $rnapath/LD011_mapping/ -tbp 1 -fragLength given -strand separate -style rnaseq


# make BigWig files
nice makeBigWig.pl $rnapath/LD001_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD002_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD005_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD006_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD007_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD008_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD009_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD010_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath
nice makeBigWig.pl $rnapath/LD011_mapping/ hg38 -strand -force -url $httppath -webdir $linkpath

# findPeaks
nice findPeaks $rnapath/LD005_mapping -style region -o $rnapath/LD005_input_histone.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD006_mapping -style region -o $rnapath/LD006_input_histone.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD010_mapping -style region -o $rnapath/LD010_input_histone.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD007_mapping -style region -o $rnapath/LD007_input_histone.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD008_mapping -style region -o $rnapath/LD008_input_histone.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD011_mapping -style region -o $rnapath/LD011_input_histone.txt -i $rnapath/LD001_mapping -tbp 1

nice findPeaks $rnapath/LD005_mapping -style factor -o $rnapath/LD005_input_TF.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD006_mapping -style factor -o $rnapath/LD006_input_TF.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD010_mapping -style factor -o $rnapath/LD010_input_TF.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD007_mapping -style factor -o $rnapath/LD007_input_TF.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD008_mapping -style factor -o $rnapath/LD008_input_TF.txt -i $rnapath/LD001_mapping -tbp 1
nice findPeaks $rnapath/LD011_mapping -style factor -o $rnapath/LD011_input_TF.txt -i $rnapath/LD001_mapping -tbp 1

# find motifs
nice findMotifsGenome.pl $rnapath/LD005_input_TF.txt hg38 -size 200 -mask -p 4
nice findMotifsGenome.pl $rnapath/LD006_input_TF.txt hg38 -size 200 -mask -p 4
nice findMotifsGenome.pl $rnapath/LD010_input_TF.txt hg38 -size 200 -mask -p 4
nice findMotifsGenome.pl $rnapath/LD007_input_TF.txt hg38 -size 200 -mask -p 4
nice findMotifsGenome.pl $rnapath/LD008_input_TF.txt hg38 -size 200 -mask -p 4
nice findMotifsGenome.pl $rnapath/LD011_input_TF.txt hg38 -size 200 -mask -p 4
