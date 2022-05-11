#set the directory paths
chippath=/NextGenSeqData/project-data/sejun/chipseq/
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/
linkpath=/webdata/bric-data/sejun/chipseq/

nice findMotifsGenome.pl $chippath/1/VSF_a_input.txt mm10 $chippath/1/motifs/VSF_a/ -size 200 -mask -p 4
nice findMotifsGenome.pl $chippath/1/VSF_b_input.txt mm10 $chippath/1/motifs/VSF_b/ -size 200 -mask -p 4
nice findMotifsGenome.pl $chippath/1/VSF_c_input.txt mm10 $chippath/1/motifs/VSF_c/ -size 200 -mask -p 4
