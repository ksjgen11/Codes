#!/bin/bash

atac=~/data/atacseq/
blacklist=$atac/hg38.blacklist.bed
peak=`echo "new_result"`
peakpath=$atac/$peak
scales=`echo "scale_factor"`
scalepath=$peakpath/$scales
# mkdir $scalepath
igv=~/bin/igv.sh
igvtools=~/bin/igvtools
web=~/web/atacseq/new/
# mkdir $web

# data normalization

for cond in si673_4h si673_12h si673_24h siNT 
do
    for rep in 1 2 3
    do 
    	fl=$(basename -a ${atac}/${cond}_${rep}_*.bam)
		scale= cat $scalepath/scale_factor_count_${cond}_${rep}_*.txt
		printf "normalizes for %s\n" "$fl"
		echo "$scale"

		bedtools genomecov -ibam $fl -bg -scale $scale > $peakpath/scaled_${cond}_${rep}_new.bedgraph
		$igvtools toTDF -z 9 $peakpath/scaled_${cond}_${rep}_new.bedgraph $peakpath/scaled_${cond}_${rep}_new.tdf hg38
        cp $peakpath/scaled_${cond}_${rep}_new.tdf $web/.
    done
done

# merging bedfiles


bedtools intersect -a $peakpath/blacklist_removed_si673_4h_1_new.bed \
 -b $peakpath/blacklist_removed_si673_4h_2_new.bed $peakpath/blacklist_removed_si673_4h_3_new.bed -sorted > \
 $peakpath/merged_si673_4h_new.bed
bedtools intersect -a $peakpath/blacklist_removed_si673_12h_1_new.bed \
 -b $peakpath/blacklist_removed_si673_12h_2_new.bed $peakpath/blacklist_removed_si673_12h_3_new.bed -sorted > \
 $peakpath/merged_si673_12h_new.bed
bedtools intersect -a $peakpath/blacklist_removed_si673_24h_1_new.bed \
 -b $peakpath/blacklist_removed_si673_24h_2_new.bed $peakpath/blacklist_removed_si673_24h_3_new.bed -sorted > \
 $peakpath/merged_si673_24h_new.bed
bedtools intersect -a $peakpath/blacklist_removed_siNT_1_new.bed \
 -b $peakpath/blacklist_removed_siNT_2_new.bed $peakpath/blacklist_removed_siNT_3_new.bed -sorted > \
 $peakpath/merged_siNT_new.bed


# find motif
for fl in $peakpath/merged_*
do
    samplename=${fl%_lower.bed}
    output="${samplename##merged_}_new"
    findMotifsGenome.pl $fl hg38 $output -S 15 -len 6,8,10,12 -size -250,250 -bg $peakpath/merged_siNT_new.bed -p 4
done
