# RNAseq analysis pipeline

## Used pipeline

Linux envrionment (Debian 4.9.272-2)
R (ver 4.0.3)
Rstudio server

## Used program

### On Linux
Fastqc # sequencing files quality check
Trimmomatic # trimming adatper sequences
STAR # mapping reads on the transcriptome

### On R
Rsubread # count reads
DEseq # normalize data, DEseq analysis
ggplot2 # draw plots

# Analysis order

## Fastqc
At first, you need to check the sequencing quality, and remove adapter sequences from reads. 
You can use fastqc by GUI format or command line format.
If you want to use GUI format, just use `fastqc` on the command line.
```{bash}
fastqc
```

You can see the report through the program. 

If you want to use command line format, use `fastqc` with other options.
```{bash}
fastqc [-o outdir] file1 file2 file3....
```

-o option set the directory where the report would save. Before use this option, make the folder first.

### Advantages
GUI format : intuitively working with it. 
Command line format : command for many files in one line. 

After run, check the qc report and decide how you trim the sequences. 

## trimmomatic

There are many programs to trim reads (fastx-toolkit, cutadopt, ...). Here, I generally use trimmomatic because it controls
paired-end data once, removing adapter sequence easier. 
Unfortunately, trimmomatic is not a basic program on tycho. If you want to use, download the program following the guidance. 
Or, use this path.
```{bash}
trimmomaticpath=~/data/Trimmomatic-0.39
```

To remove adapter sequence, adapter sequence-reference file is needed. It is in genome folder on tycho.
```{bash}
adapter=/k/genomes/adapters/fa/adapters.fa
```

Trimmomatic is based on java. If you want to use your local computer, please check the java is available.

```{bash}
java -jar $trimmomaticpath/trimmomatic-0.39.jar PE -threads 4 -phred33 \
read1.fastq.gz read2.fastq.gz \
trimmed_read1.fa untrimmed_read1.fa \
trimmed_read2.fa untrimmed_read2.fa \
CROP:{tail} HEADCROP:{head} ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:{min_length}
```
Please get .fa file to next step. STAR cannot accept fastq.gz file.

PE means paired-end file. If you want to use single-end file, use SE.
SE option is only available trimmed. untrimmed is not available.
-thread : multithreads to make your work faster.
-phred33 : only accept Q > 33. 
CROP : trimming tail
HEADCROP : trimming head
ILLUMINACLIP : trimming adapter
MINLEN : minimum length criteria of read. 

After trimming, do fastqc one more to check the trimming result.

## mapping

STAR is very fast transciptome mapping program. It needs STAR_index for it. the indices are in the genome folder on tycho. Please check the tycho instruction. 

```{bash}
STAR --genomeDir $indexpath --runThreadN 24 --readFilesIn read1.fa read2.fa \
    --outFileNamePrefix output_prefix --outFilterIntronMotifs RemoveNoncanonicalUnannotated
```

--genomeDir : the index path
--runThreadN : multithreads
--readFilesIn : reads for mapping
--outFileNamePrefix : prefix for output sample.
--outFilterIntronMotifs : remove reads if it is mapped intron range.

The output is shown {prefix}AlignedOut.sam. The output folder should be exist before running.

To handle conveniently, convert sam file to bam file by samtools.

```{bash}
samtools view --threads 16 -Sb {output_prefix}Aligned.out.sam > output.bam
```

## count reads, normalization, DEseq --> Please check rmd file.