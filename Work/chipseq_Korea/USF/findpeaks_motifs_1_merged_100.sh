chippath=/NextGenSeqData/project-data/sejun/chipseq/

nice findPeaks $chippath/1/VSF_a/ -style factor -o $chippath/1/VSF_a_input.txt -i $chippath/1/IN_a -tbp 1\
 -fragLength 100
nice findPeaks $chippath/1/VSF_b/ -style factor -o $chippath/1/VSF_b_input.txt -i $chippath/1/IN_a -tbp 1\
 -fragLength 100
nice findPeaks $chippath/1/VSF_c/ -style factor -o $chippath/1/VSF_c_input.txt -i $chippath/1/IN_b -tbp 1\
 -fragLength 100


nice findMotifsGenome.pl $chippath/1/VSF_a_input.txt mm10 $chippath/1/motifs/VSF_a_100/ -size -100,100 -mask -p 4
nice findMotifsGenome.pl $chippath/1/VSF_b_input.txt mm10 $chippath/1/motifs/VSF_b_100/ -size -100,100 -mask -p 4
nice findMotifsGenome.pl $chippath/1/VSF_c_input.txt mm10 $chippath/1/motifs/VSF_c_100/ -size -100,100 -mask -p 4


