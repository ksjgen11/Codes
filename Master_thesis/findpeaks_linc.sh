rnapath=/NextGenSeqData/project-data/sejun/rnaseq/

nice findPeaks $rnapath/LD005/ -style tss -o $rnapath/LD005_LD001.txt -i $rnapath/LD001 -tbp 1\
 -fragLength 100
nice findPeaks $rnapath/LD006/ -style tss -o $rnapath/LD006_LD002.txt -i $rnapath/LD002 -tbp 1\
 -fragLength 100
nice findPeaks $rnapath/LD010/ -style tss -o $rnapath/LD010_LD009.txt -i $rnapath/LD009 -tbp 1\
 -fragLength 100
nice findPeaks $rnapath/LD007/ -style tss -o $rnapath/LD007_LD001.txt -i $rnapath/LD001 -tbp 1\
 -fragLength 100
nice findPeaks $rnapath/LD008/ -style tss -o $rnapath/LD008_LD002.txt -i $rnapath/LD002 -tbp 1\
 -fragLength 100
nice findPeaks $rnapath/LD011/ -style tss -o $rnapath/LD011_LD009.txt -i $rnapath/LD009 -tbp 1\
 -fragLength 100
