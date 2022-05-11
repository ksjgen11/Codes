index_hg19=/k/genomes/hg19/index/bowtie2_canonical/hg19
trimmed=~/data/rnaseq/trimmed/

tophat2 --num-threads 4 -N 4 \
 -o LD001_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD001 --rg-sample LD001 --rg-library LD001 $index_hg19 \
 $trimmed/LD001_trimmed.fa
tophat2 --num-threads 4 -N 4 \
 -o LD002_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD002 --rg-sample LD002 --rg-library LD002 $index_hg19 \
 $trimmed/LD002_trimmed.fa
tophat2 --num-threads 4 -N 4 \
 -o LD005_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD005 --rg-sample LD005 --rg-library LD005 $index_hg19 \
 $trimmed/LD005_trimmed.fa
tophat2 --num-threads 4 -N 4 \
 -o LD006_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD006 --rg-sample LD006 --rg-library LD006 $index_hg19 \
 $trimmed/LD006_trimmed.fa
tophat2 --num-threads 4 -N 4 \
 -o LD007_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD007 --rg-sample LD007 --rg-library LD007 $index_hg19 \
 $trimmed/LD007_trimmed.fa
tophat2 --num-threads 4 -N 4 \
 -o LD008_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD008 --rg-sample LD008 --rg-library LD008 $index_hg19 \
 $trimmed/LD008_trimmed.fa
tophat2 --num-threads 4 -N 4 \
 -o LD009_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD009 --rg-sample LD009 --rg-library LD009 $index_hg19 \
 $trimmed/LD009_trimmed.fa
tophat2 --num-threads 4 -N 4 \
 -o LD010_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD010 --rg-sample LD010 --rg-library LD010 $index_hg19 \
 $trimmed/LD010_trimmed.fa
tophat2 --num-threads 4 -N 4 \
 -o LD011_tophat_mapping --max-multihits 10 --no-coverage-search --read-mismatches 2 \
 --rg-id LD011 --rg-sample LD011 --rg-library LD011 $index_hg19 \
 $trimmed/LD011_trimmed.fa
