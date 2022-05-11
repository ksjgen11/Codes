#!/bin/bash

atac=~/data/atacseq/
blacklist=$atac/hg38.blacklist.bed
peak=`echo "new_result"`
# mkdir $atac/$peak
peakpath=$atac/$peak
scales=`echo "scale_factor"`
# mkdir $peakpath/$scales
scalepath=$peakpath/$scales
igv=~/bin/igv.sh
igvtools=~/bin/igvtools
# mkdir ~/web/atacseq/new
web=~/web/atacseq/new/


# functions
# peak calling_threshold 0.005
# blacklist removing

for cond in si673_4h si673_12h si673_24h siNT 
do
    for rep in 1 2 3
    do 
    	fl=$(basename -a ${atac}/${cond}_${rep}_*.bam)
		printf "Calling peaks for %s\n" "$fl"
		
		macs2 callpeak \
		-t $fl \
		-n ${cond}_${rep}_new \
		--outdir $peakpath \
		-f BAMPE \
		-g hs \
		-q 0.005 \
		--nomodel \
		--shift 0 \

		printf "removing blacklist for %s\n" "${cond}_${rep}_new_peaks.xls"		
        sed '/^$/d' $peakpath/${cond}_${rep}_new_peaks.xls > $peakpath/peaks.bed
        sed '/#/d' $peakpath/peaks.bed > $peakpath/peaks1.bed
		sed '1d' $peakpath/peaks1.bed > $peakpath/peaks.bed
        awk '{ print $1"\t"$2"\t"$3"\t"$10"\t"$7}' $peakpath/peaks.bed > $peakpath/peaks1.bed
        mv $peakpath/peaks1.bed $peakpath/peaks.bed
        bedtools intersect -v -a $peakpath/peaks.bed -b $blacklist > $peakpath/blacklist_removed_${cond}_${rep}_new.bed
        rm $peakpath/peaks.bed
		
	done
done


# merge data
cat $peakpath/blacklist_removed_* >> $peakpath/total.bed
bedtools sort -i $peakpath/total.bed > $peakpath/total_sorted.bed
bedtools merge -i $peakpath/total_sorted.bed > $peakpath/merge.bed
awk '{FS="\t";OFS="\t";$4="merged"NR;print}' $peakpath/merge.bed > $peakpath/merge1.bed
awk '{FS="\t";OFS="\t";$6=".";$7="+";$8=".";$9="gene_id \""$4"\"";$4=$2;$5=$3;$2="merged";$3="exon";print}' \
	$peakpath/merge1.bed > $peakpath/merge_new.gff
rm $peakpath/merge.bed

# peak annotation file 

grep 'chr' $peakpath/merge1.bed > $peakpath/merge1_deleted_new.bed


# peak intensity count
for cond in si673_4h si673_12h si673_24h siNT 
do
    for rep in 1 2 3
    do 
    	fl=$(basename -a ${atac}/${cond}_${rep}_*.bam)
		printf "count peak intnesities for %s\n" "$fl"

        htseq-count -f bam -r pos -m intersection-nonempty $fl $peakpath/merge_new.gff > \
        $peakpath/count_${cond}_${rep}_new.txt
		
	done
done

Rscript atac_higher_threshold_normalization.R

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
    samplename=${fl%_new.bed}
    output="${samplename##merged_}_new"
    findMotifsGenome.pl $fl hg38 $output -S 15 -len 6,8,10,12 -size -250,250 -bg $peakpath/merged_siNT_new.bed -p 4
done


bedtools intersect -a $peakpath/blacklist_removed_si673_12h_1_new.bed \
 -b $peakpath/blacklist_removed_si673_12h_2_new.bed -sorted > \
 $peakpath/merged_si673_12h_new_deleted.bed

findMotifsGenome.pl $peakpath/merged_si673_12h_new_deleted.bed hg38 merged_si673_12h_new_deleted \
 -S 15 -len 6,8,10,12 -size -250,250 -bg $peakpath/merged_siNT_new.bed -p 4
