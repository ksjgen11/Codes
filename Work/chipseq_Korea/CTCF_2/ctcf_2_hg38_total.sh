#set the directory paths
chippath=/NextGenSeqData/project-data/sejun/chipseq/ctcf_2
datapath=/NextGenSeqData/project-data/sejun/chipseq/ctcf_2/trimmed
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/ctcf_2
linkpath=/webdata/bric-data/sejun/chipseq/ctcf_2
indexpath=/k/genomes/hg38/index/bowtie2_canonical/hg38
input=~/data/chipseq/2/Input

# mapping
bowtie2 -p 8 -x $indexpath \
 -1 $datapath/trimmed.4C-B_1.fa -2 $datapath/trimmed.4C-B_2.fa \
 -S $chippath/4C-B.sam --no-unal --local --phred33 -t

bowtie2 -p 8 -x $indexpath \
 -1 $datapath/trimmed.CTCF_1.fa -2 $datapath/trimmed_CTCF_2.fa \
 -S $chippath/CTCF_2.sam --no-unal --local --phred33 -t

# making tagdirectory
nice makeTagDirectory $chippath/4C-B $chippath/4C-B.sam -format sam \
 -tbp 1 -unique -genome hg38

nice makeTagDirectory $chippath/CTCF_2 $chippath/CTCF_2.sam -format sam \
 -tbp 1 -unique -genome hg38

# remove outofbounds
nice removeOutOfBoundsReads.pl $chippath/4C-B/ hg38
nice removeOutOfBoundsReads.pl $chippath/CTCF_2/ hg38

# make UCSC files
nice makeUCSCfile $chippath/4C-B/ -tbp 1 -style chipseq

nice makeUCSCfile $chippath/CTCF_2/ -tbp 1 -style chipseq

# make bigwig files
nice makeBigWig.pl $chippath/4C-B/ hg38 -normal -force -url $httppath -webdir $linkpath

nice makeBigWig.pl $chippath/CTCF_2/ hg38 -normal -force -url $httppath -webdir $linkpath

# find peaks
nice findPeaks $chippath/4C-B/ -style factor -o $chippath/4C-B_hg38_input.txt -i $input -tbp 1\
 -fragLength 100

nice findPeaks $chippath/CTCF_2/ -style factor -o $chippath/CTCF_2_hg38_input.txt -i $input -tbp 1\
 -fragLength 100

# find motifs
nice findMotifsGenome.pl $chippath/4C-B_hg38_input.txt hg38 $chippath/4C-B_hg38_motifs/ \
 -size 200 -mask -p 4

 nice findMotifsGenome.pl $chippath/CTCF_2_hg38_input.txt hg38 $chippath/CTCF_2_hg38_motifs/ \
 -size 200 -mask -p 4