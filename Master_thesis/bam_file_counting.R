library(tidyverse)
library(airway)
library(Rsamtools)
library(GenomicFeatures)
library("GenomicAlignments")
library("BiocParallel")
library(Rsubread)

setwd('d:/Projects/UCPH/Work/rnaseq/')
dir <- 'bam files/'
file_list <- list.files(path = dir, pattern = '^LD0..')
filenames <- file.path(path = dir, file_list)
file.exists(filenames)
register(SnowParam(4))
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
seqinfo(bamfiles[1])



gtffile <- file.path('Homo_sapiens.GRCh37.87.gtf.gz')
txdb <- makeTxDbFromGFF(gtffile, format = 'gtf')



ebg <- exonsBy(txdb, by = 'gene')
ebg


se <- summarizeOverlaps(features=ebg, 
                         reads=bamfiles,
                         mode="Union",
                         singleEnd=T,
                         ignore.strand=F)


se
assay(se)


fc <- featureCounts(files=filenames, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=F)
fc$counts
fc$annotation
fc$targets
fc$stat

design <- data.frame(condition = c('NT', 'NT', 'LNC1', 'LNC1', 'LNC2', 'LNC2', 'NT', 'LNC1', 'LNC2'))
rownames(design) <- colnames(fc$counts)
identical(assay(se), fc$counts)
count1 <- assay(se)
count2 <- fc$counts

library(DESeq2)
dds1 <- DESeqDataSetFromMatrix(countData = count1, 
                               colData = design, 
                               design = ~condition)
deseq1 <- DESeq(dds1)

# sort result with logFC threshold 0.25, padj < 0.05
result1 <- results(deseq1, lfcThreshold = 0.25, alpha = 0.05)

# confirm result
summary(result1)

# make tibble from DESeq result
result_tb1 <- result1 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

# make a plot
ggplot(result_tb1, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = padj<0.05), alpha = 0.5) +
  scale_x_log10() +
  geom_smooth() +
  geom_hline(yintercept = 0, alpha = 0.75, 
             color = 'darkgreen', linetype = 'dashed') +
  theme_bw()


ggplot(result_tb1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.75, 
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = 1) +
  theme_bw()





# 
# dds2 <- DESeqDataSetFromMatrix(countData = count2, 
#                                colData = design, 
#                                design = ~condition)
# deseq2 <- DESeq(dds2)
# 
# # sort result with logFC threshold 0.25, padj < 0.05
# result2 <- results(deseq2, lfcThreshold = 0.25, alpha = 0.05)
# 
# # confirm result
# summary(result2)
# 
# # make tibble from DESeq result
# result_tb2 <- result2 %>%
#   as.data.frame() %>%
#   rownames_to_column('gene') %>%
#   as_tibble()
# 
# # make a plot
# ggplot(result_tb2, aes(x = baseMean, y = log2FoldChange)) +
#   geom_point(aes(color = padj<0.05), alpha = 0.5) +
#   scale_x_log10() +
#   geom_smooth() +
#   geom_hline(yintercept = 0, alpha = 0.75, 
#              color = 'darkgreen', linetype = 'dashed') +
#   theme_bw()
# 
# 
# ggplot(result_tb2, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(color = padj < 0.05), alpha = 0.5) +
#   geom_vline(xintercept = 0, alpha = 0.75, 
#              color = 'black', linetype = 'dashed') +
#   geom_hline(yintercept = 1) +
#   theme_bw()

vsd <- vst(dds1)
class(vsd)
colData(vsd)
plotPCA(vsd, 'condition')

