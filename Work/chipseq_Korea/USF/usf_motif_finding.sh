chippath=/NextGenSeqData/project-data/sejun/chipseq/1/
motifs=$chippath/motifs/VSF_a/homerResults

nice annotatePeaks.pl $chippath/VSF_a_input.txt mm10 -size -500,100 -m /usr/local/src/homer/motifs/usf2.motif \
  > USF_motif_usf2_500_100.txt
