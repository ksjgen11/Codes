#set the directory paths
indexpath=/NextGenSeqData/project-data/sejun/chipseq/2/index/hg19_ebv/
chippath=/NextGenSeqData/project-data/sejun/chipseq/2
datapath=/NextGenSeqData/project-data/sejun/chipseq/2/2/trimmed
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/2
linkpath=/webdata/bric-data/sejun/chipseq/2

# making tagdirectory
nice makeTagDirectory $chippath/CTCF_ebv_only $chippath/CTCF_ebv_only.sam -format sam \
 -tbp 1 -unique

nice makeTagDirectory $chippath/Input_ebv_only $chippath/Input_ebv_only.sam -format sam \
 -tbp 1 -unique


# make UCSC files
nice makeUCSCfile $chippath/CTCF_ebv_only/ -tbp 1 -style chipseq -o auto

nice makeUCSCfile $chippath/Input_ebv_only/ -tbp 1 -style chipseq

# make bigwig files
nice makeBigWig.pl $chippath/CTCF_ebv_only/ $indexpath/ebv -normal -force -url $httppath -webdir $linkpath

nice makeBigWig.pl $chippath/Input_ebv_only/ $indexpath/ebv -normal -force -url $httppath -webdir $linkpath

# find peaks
nice findPeaks $chippath/CTCF_ebv_only/ -style factor -o $chippath/ctcf_ebv_only_input.txt -i $chippath/Input_ebv -tbp 1\
-fragLength 100

# find motifs
nice findMotifsGenome.pl $chippath/ctcf_ebv_only_input.txt $indexpath/ebv $chippath/ctcf_ebv_motifs_2/ \
 -size 200 -mask -p 4

