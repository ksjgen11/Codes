chippath=/NextGenSeqData/project-data/sejun/chipseq/1/
motifs=$chippath/motifs/VSF_a/homerResults

nice annotatePeaks.pl $chippath/VSF_a_input.txt mm10 -size 200 -m $motifs/motif1.motif \
 -mbed $chippath/VSF_a_motif.bed > $chippath/VSF_a_motif_annotation.txt
nice annotatePeaks.pl $chippath/VSF_b_input.txt mm10 -size 200 -m $motifs/motif1.motif \
 -mbed $chippath/VSF_b_motif.bed > $chippath/VSF_b_motif_annotation.txt
nice annotatePeaks.pl $chippath/VSF_c_input.txt mm10 -size 200 -m $motifs/motif1.motif \
 -mbed $chippath/VSF_c_motif.bed > $chippath/VSF_c_motif_annotation.txt

nice annotatePeaks.pl $chippath/VSF_7_input.txt mm10 -size 200 -m $motifs/motif1.motif \
 -mbed $chippath/VSF_7_motif.bed > $chippath/VSF_7_motif_annotation.txt
nice annotatePeaks.pl $chippath/VSF_8_input.txt mm10 -size 200 -m $motifs/motif1.motif \
 -mbed $chippath/VSF_8_motif.bed > $chippath/VSF_8_motif_annotation.txt
nice annotatePeaks.pl $chippath/VSF_9_input.txt mm10 -size 200 -m $motifs/motif1.motif \
-mbed $chippath/VSF_9_motif.bed > $chippath/VSF_9_motif_annotation.txt
nice annotatePeaks.pl $chippath/VSF_10_input.txt mm10 -size 200 -m $motifs/motif1.motif \
 -mbed $chippath/VSF_10_motif.bed > $chippath/VSF_10_motif_annotation.txt
nice annotatePeaks.pl $chippath/VSF_11_input.txt mm10 -size 200 -m $motifs/motif1.motif \
 -mbed $chippath/VSF_11_motif.bed > $chippath/VSF_11_motif_annotation.txt
nice annotatePeaks.pl $chippath/VSF_12_input.txt mm10 -size 200 -m $motifs/motif1.motif \
 -mbed $chippath/VSF_12_motif.bed > $chippath/VSF_12_motif_annotation.txt
