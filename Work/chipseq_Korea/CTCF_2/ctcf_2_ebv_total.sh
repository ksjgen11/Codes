#set the directory paths
indexpath=/NextGenSeqData/project-data/sejun/chipseq/2/index/ebv
chippath=/NextGenSeqData/project-data/sejun/chipseq/ctcf_2
datapath=/NextGenSeqData/project-data/sejun/chipseq/ctcf_2/trimmed
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/ctcf_2
linkpath=/webdata/bric-data/sejun/chipseq/ctcf_2
input=~/data/chipseq/2/Input_ebv_only


# mapping
bowtie2 -p 8 -x $indexpath/ebv \
 -1 $datapath/trimmed.4C-B_1.fa -2 $datapath/trimmed.4C-B_2.fa \
 -S $chippath/4C-B_ebv.sam --no-unal --local --phred33 -t

bowtie2 -p 8 -x $indexpath/ebv \
 -1 $datapath/trimmed.CTCF_1.fa -2 $datapath/trimmed_CTCF_2.fa \
 -S $chippath/CTCF_ebv_2.sam --no-unal --local --phred33 -t

# making tagdirectory
nice makeTagDirectory $chippath/4C-B_ebv $chippath/4C-B_ebv.sam -format sam \
 -tbp 1 -unique -genome $indexpath/ebv

nice makeTagDirectory $chippath/CTCF_ebv_2 $chippath/CTCF_ebv_2.sam -format sam \
 -tbp 1 -unique -genome $indexpath/ebv

# remove outofbounds
nice removeOutOfBoundsReads.pl $chippath/4C-B_ebv/ $indexpath/ebv.fa
nice removeOutOfBoundsReads.pl $chippath/CTCF_ebv_2/ $indexpath/ebv.fa

# make UCSC files
nice makeUCSCfile $chippath/4C-B_ebv/ -tbp 1 -style chipseq

nice makeUCSCfile $chippath/CTCF_ebv_2/ -tbp 1 -style chipseq

# make bigwig files
nice makeBigWig.pl $chippath/4C-B_ebv/ $indexpath/ebv.fa -normal -force -url $httppath -webdir $linkpath

nice makeBigWig.pl $chippath/CTCF_ebv_2/ $indexpath/ebv.fa -normal -force -url $httppath -webdir $linkpath

# find peaks
nice findPeaks $chippath/4C-B_ebv/ -style factor -o $chippath/4C-B_ebv_input.txt -i $input -tbp 1\
-fragLength 100

nice findPeaks $chippath/CTCF_ebv_2/ -style factor -o $chippath/CTCF_2_ebv_input.txt -i $input -tbp 1\
-fragLength 100

# find motifs
nice findMotifsGenome.pl $chippath/4C-B_ebv_input.txt $indexpath/ebv $chippath/4C-B_ebv_motifs/ \
 -size 200 -mask -p 4

 nice findMotifsGenome.pl $chippath/CTCF_2_ebv_input.txt $indexpath/ebv $chippath/CTCF_2_ebv_motifs/ \
 -size 200 -mask -p 4