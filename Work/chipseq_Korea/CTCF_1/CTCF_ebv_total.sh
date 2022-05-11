#set the directory paths
indexpath=/NextGenSeqData/project-data/sejun/chipseq/2/index/hg19_ebv/
chippath=/NextGenSeqData/project-data/sejun/chipseq/2
datapath=/NextGenSeqData/project-data/sejun/chipseq/2/2/trimmed
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/2
linkpath=/webdata/bric-data/sejun/chipseq/2

# mapping
bowtie2 -p 8 -x $indexpath/hg19_ebv \
 -1 $datapath/CTCF-ab_1_trimmed.fa -2 $datapath/CTCF-ab_2_trimmed.fa \
 -S $chippath/CTCF_ebv.sam --no-unal --local --phred33 -t

bowtie2 -p 8 -x $indexpath/hg19_ebv \
 -1 $datapath/Input_1_trimmed.fa -2 $datapath/Input_2_trimmed.fa \
 -S $chippath/Input_ebv.sam --no-unal --local --phred33 -t

# making tagdirectory
nice makeTagDirectory $chippath/CTCF_ebv $chippath/CTCF_ebv.sam -format sam \
 -tbp 1 -unique -genome $indexpath/hg19_ebv

nice makeTagDirectory $chippath/Input_ebv $chippath/Input_ebv.sam -format sam \
 -tbp 1 -unique -genome $indexpath/hg19_ebv

# remove outofbounds
nice removeOutOfBoundsReads.pl $chippath/CTCF_ebv/ $indexpath/hg19_ebv

nice removeOutOfBoundsReads.pl $chippath/Input_ebv/ $indexpath/hg19_ebv

# make UCSC files
nice makeUCSCfile $chippath/CTCF_ebv/ -tbp 1 -style chipseq

nice makeUCSCfile $chippath/Input_ebv/ -tbp 1 -style chipseq

# make bigwig files
nice makeBigWig.pl $chippath/CTCF_ebv/ $indexpath/hg19_ebv -normal -force -url $httppath -webdir $linkpath

nice makeBigWig.pl $chippath/Input_ebv/ $indexpath/hg19_ebv -normal -force -url $httppath -webdir $linkpath

# find peaks
nice findPeaks $chippath/CTCF_ebv/ -style factor -o $chippath/ctcf_ebv_input.txt -i $chippath/Input_ebv -tbp 1\
-fragLength 100

# find motifs
nice findMotifsGenome.pl $chippath/ctcf_ebv_input.txt $indexpath/hg19_ebv $chippath/ctcf_ebv_motifs/ \
 -size 200 -mask -p 4

