#!/bin/bash

atac=~/data/atacseq/
blacklist=$atac/hg38.blacklist.bed
peak=`echo "new_result"`
mkdir $atac/$peak
peakpath=$atac/$peak
igv=~/bin/igv.sh
igvtools=~/bin/igvtools

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

