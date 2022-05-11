library(airway)
setwd('~/data/rnaseq/bamfiles')
dir <- '~/data/rnaseq/bamfiles'
file_list <- list.files(path = dir, pattern = '^LD0..')
filenames <- file.path(path = dir, file_list)
file.exists(file_list)

library(Rsamtools)
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
seqinfo(bamfiles[1])

library(GenomicFeatures)

gtffile <- file.path('~/data/rnaseq', 'Homo_sapiens.GRCh38.104.gtf.gz')
txdb <- makeTxDbFromGFF(gtffile, format = 'gtf')

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

ebg <- exonsBy(txdb, by = 'gene')
ebg

ebg2 <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = 'gene')

library("GenomicAlignments")
library("BiocParallel")
register(MulticoreParam(8))
se2 <- summarizeOverlaps(features=ebg2, 
                        reads=bamfiles,
                        mode="Union",
                        singleEnd=T,
                        ignore.strand=TRUE)
se2
tail(assay(se))
library(tidyverse)
assay(se2) %>%
  filter(LD001_human.bam != 0)
assay(se2)[assay(se2) != 0]
assay(se)[assay(se) != 0]

se3 <- summarizeOverlaps(features=ebg, 
                         reads=bamfiles,
                         mode="Union",
                         singleEnd=T,
                         ignore.strand=F)


