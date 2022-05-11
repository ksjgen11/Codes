#set the directory paths

chippath=/NextGenSeqData/project-data/sejun/chipseq/2
datapath=/NextGenSeqData/project-data/sejun/chipseq/2/2/trimmed
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/2
linkpath=/webdata/bric-data/sejun/chipseq/2

indexpath=/NextGenSeqData/project-data/sejun/chipseq/2/index/ebv

bowtie2 -p 8 -x $indexpath/ebv \
 -1 $datapath/CTCF-ab_1_trimmed.fa -2 $datapath/CTCF-ab_2_trimmed.fa \
 -S $chippath/CTCF_ebv_only.sam --no-unal --local --phred33 -t

bowtie2 -p 8 -x $indexpath/ebv \
 -1 $datapath/Input_1_trimmed.fa -2 $datapath/Input_2_trimmed.fa \
 -S $chippath/Input_ebv_only.sam --no-unal --local --phred33 -t

# making tagdirectory
nice makeTagDirectory $chippath/CTCF_ebv_only $chippath/CTCF_ebv_only.sam -format sam \
 -tbp 1 -unique

nice makeTagDirectory $chippath/Input_ebv_only $chippath/Input_ebv_only.sam -format sam \
 -tbp 1 -unique


# make UCSC files
nice makeUCSCfile $chippath/CTCF_ebv_only/ -tbp 1 -style chipseq -o auto

nice makeUCSCfile $chippath/Input_ebv_only/ -tbp 1 -style chipseq -o auto

# make bigwig files
nice makeBigWig.pl $chippath/CTCF_ebv_only/ $indexpath/ebv -normal -force -url $httppath -webdir $linkpath

nice makeBigWig.pl $chippath/Input_ebv_only/ $indexpath/ebv -normal -force -url $httppath -webdir $linkpath

# find peaks
nice findPeaks $chippath/CTCF_ebv_only/ -style factor -o $chippath/ctcf_ebv_only_input.txt \
-i $chippath/Input_ebv_only -tbp 1 -fragLength 100

# find motifs
nice findMotifsGenome.pl $chippath/ctcf_ebv_only_input.txt $indexpath/ebv.fa $chippath/ctcf_ebv_motifs_2/ \
 -size 200 -mask -p 4

