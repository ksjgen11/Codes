library(tidyverse)
library(Rsamtools)
library(Rsubread)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(factoextra)
library(ggfortify)
library(ggrepel)
library(Mus.musculus)
library(VennDiagram)

figure <- '~/data/rnaseq_arnes/results/figures/'
dataset <- '~/data/rnaseq_arnes/results/datasets/'
figures <- '~/data/rnaseq_arnes/results/figures_wo_legend/'


# DESeq2 object
# DESeq2
design <- data.frame(condition = c('Group1', 'Group1', 
                                   'Group2', 'Group2', 'Group2', 
                                   'Group3', 'Group3', 'Group4', 
                                   'Group4', 'Group5', 'Group5', 
                                   'Group3', 'Group4', 'Group6', 
                                   'Group6', 'Group6', 'Group7', 
                                   'Group7', 'Group1')) # KL003 is the last column

count <- read.table(paste0(dataset, 'raw_count.txt'))

colnames(count) <- c('KL001', 'KL002', 'KL004', 'KL005', 'KL006', 
                     'KL007', 'KL008', 'KL009', 'KL010', 'KL011', 'KL012', 
                     'KL013', 'KL014', 'KL015', 'KL016', 'KL017', 'KL018', 
                     'KL019', 'KL003') # KL003 is the last column
rownames(design) <- colnames(count)


dds <- DESeqDataSetFromMatrix(countData = count,
                               colData = design,
                               design = ~condition)
deseq <- DESeq(dds)

# sort result with logFC threshold 0.25, padj < 0.05
result <- results(deseq, lfcThreshold = 0.25, alpha = 0.05)


# confirm result
summary(result)

# make tibble from DESeq result
result_tb <- result %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

# annotation result


result_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_tb$entrez <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

result_tb$genename <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

symbol <- result_tb$symbol

# Data transfrom

vsd <- vst(dds) # variance stabilizing transformation

# heatmap sample dist

sampleDists <- dist( t( assay(vsd) ) )
sampleDistMatrix <- as.matrix( sampleDists ) 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
distheatmap <- pheatmap(sampleDistMatrix, 
                        # main = "Sample distance clustering",
                        # annotation_legend = F,
                        annotation_col = design
                        )
distheatmap
ggsave(distheatmap, file = paste0(figure, 'rnaseq_before_sample_dist_heatmap.png'))

# sample correlation plot

# calculate correlation between pairs of samples for all possible pairs.
sampleCor <- cor(assay(vsd))
# generate heatmap
samplecorheatmap <- pheatmap(
  as.matrix(sampleCor),
  # main = "Correlation clustering",
  width = 12,
  height = 10, 
  # annotation_legend = F,
  annotation_col = design
  )
samplecorheatmap
ggsave(samplecorheatmap, file = paste0(figure, 'rnaseq_before_sample_correlation_heatmap.png'))


# PCA

pca <- prcomp(t(assay(vsd)), center=TRUE, scale.=FALSE)
fviz_eig(pca) # elbow method

pca.x <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column('condition') %>%
  as_tibble()

pca.var <- summary(pca)$importance['Proportion of Variance',] * 100


pca_2 <- ggplot(pca.x, aes(x = PC1, y = PC2, color = design$condition, 
                           label = condition)) +
  geom_point() + 
  xlab(label = paste('PC1', 'with variance', pca.var[1])) +
  ylab(label = paste('PC2', 'with variance', pca.var[2])) +
  geom_label_repel() +
  # scale_color_manual(values = colors) +
  # labs(color = 'Condition') +
  theme_bw() 
# +
#   theme(legend.position = 'none')

pca_2

ggsave(pca_2, file = paste0(figure, 'rnaseq_before_pca_plot.png'))

# heatmap of genes for higher padj values


x<- assay(vsd)
row.names(x) <- symbol
mat <- x[head(order(result_tb$padj),20), ]

mat <- mat - rowMeans(mat)
adjheatmap <- pheatmap(mat, annotation_col=design, 
                       # annotation_legend = F,
                       cellheight = 10
                       )
adjheatmap
ggsave(adjheatmap, file = paste0(figure, 'rnaseq_before_gene_significant_heatmap.png'))

# pairwise comparison based on PCA

result_1_vs_3 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group3', 'Group1'))

result_1_vs_6 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group6', 'Group1'))

result_1_vs_4 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group4', 'Group1'))

result_1_vs_7 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group7', 'Group1'))

result_2_vs_3 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group3', 'Group2'))

result_2_vs_6 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group6', 'Group2'))

result_2_vs_4 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group4', 'Group2'))

result_2_vs_7 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group7', 'Group2'))

result_5_vs_3 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group3', 'Group5'))

result_5_vs_6 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group6', 'Group5'))

result_5_vs_4 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group4', 'Group5'))

result_5_vs_7 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group7', 'Group5'))

result_3_vs_4 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group4', 'Group3'))

result_3_vs_7 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group7', 'Group3'))

result_6_vs_4 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group4', 'Group6'))

result_6_vs_7 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group7', 'Group6'))


summary(result_1_vs_3)
summary(result_1_vs_6)
summary(result_1_vs_4)
summary(result_1_vs_7)
summary(result_2_vs_3)
summary(result_2_vs_6)
summary(result_2_vs_4)
summary(result_2_vs_7)
summary(result_5_vs_3)
summary(result_5_vs_6)
summary(result_5_vs_4)
summary(result_5_vs_7)
summary(result_3_vs_4)
summary(result_3_vs_7)
summary(result_6_vs_4)
summary(result_6_vs_7)

result_1_vs_3_tb <- result_1_vs_3 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_1_vs_6_tb <- result_1_vs_6 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_1_vs_4_tb <- result_1_vs_4 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_1_vs_7_tb <- result_1_vs_7 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()


result_2_vs_3_tb <- result_2_vs_3 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_2_vs_6_tb <- result_2_vs_6 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_2_vs_4_tb <- result_2_vs_4 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_2_vs_7_tb <- result_2_vs_7 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_5_vs_3_tb <- result_5_vs_3 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_5_vs_6_tb <- result_5_vs_6 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_5_vs_4_tb <- result_5_vs_4 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_5_vs_7_tb <- result_5_vs_7 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_3_vs_4_tb <- result_3_vs_4 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_3_vs_7_tb <- result_3_vs_7 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_6_vs_4_tb <- result_6_vs_4 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_6_vs_7_tb <- result_6_vs_7 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

# annotate genes each result

result_1_vs_3_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_1_vs_6_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_1_vs_4_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_1_vs_7_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")


result_2_vs_3_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_2_vs_6_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_2_vs_4_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_2_vs_7_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_5_vs_3_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_5_vs_6_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_5_vs_4_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_5_vs_7_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_3_vs_4_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_3_vs_7_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_6_vs_4_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

result_6_vs_7_tb$symbol <- mapIds(Mus.musculus,
                     keys=result_tb$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")


write.csv(result_1_vs_3_tb, file = paste0(dataset, 'deseq_result_1_vs_3.csv'), row.names = FALSE)
write.csv(result_1_vs_6_tb, file = paste0(dataset, 'deseq_result_1_vs_6.csv'), row.names = FALSE)
write.csv(result_1_vs_4_tb, file = paste0(dataset, 'deseq_result_1_vs_4.csv'), row.names = FALSE)
write.csv(result_1_vs_7_tb, file = paste0(dataset, 'deseq_result_1_vs_7.csv'), row.names = FALSE)

write.csv(result_2_vs_3_tb, file = paste0(dataset, 'deseq_result_2_vs_3.csv'), row.names = FALSE)
write.csv(result_2_vs_6_tb, file = paste0(dataset, 'deseq_result_2_vs_6.csv'), row.names = FALSE)
write.csv(result_2_vs_4_tb, file = paste0(dataset, 'deseq_result_2_vs_4.csv'), row.names = FALSE)
write.csv(result_2_vs_7_tb, file = paste0(dataset, 'deseq_result_2_vs_7.csv'), row.names = FALSE)

write.csv(result_5_vs_3_tb, file = paste0(dataset, 'deseq_result_5_vs_3.csv'), row.names = FALSE)
write.csv(result_5_vs_6_tb, file = paste0(dataset, 'deseq_result_5_vs_6.csv'), row.names = FALSE)
write.csv(result_5_vs_4_tb, file = paste0(dataset, 'deseq_result_5_vs_4.csv'), row.names = FALSE)
write.csv(result_5_vs_7_tb, file = paste0(dataset, 'deseq_result_5_vs_7.csv'), row.names = FALSE)

write.csv(result_3_vs_4_tb, file = paste0(dataset, 'deseq_result_3_vs_4.csv'), row.names = FALSE)
write.csv(result_3_vs_7_tb, file = paste0(dataset, 'deseq_result_3_vs_7.csv'), row.names = FALSE)

write.csv(result_6_vs_4_tb, file = paste0(dataset, 'deseq_result_6_vs_4.csv'), row.names = FALSE)
write.csv(result_6_vs_7_tb, file = paste0(dataset, 'deseq_result_6_vs_7.csv'), row.names = FALSE)



# volcano plot

result_1_vs_3_tb <- read.csv(paste0(dataset, 'deseq_result_1_vs_3.csv')) 
result_1_vs_6_tb <- read.csv(paste0(dataset, 'deseq_result_1_vs_6.csv'))
result_1_vs_4_tb <- read.csv(paste0(dataset, 'deseq_result_1_vs_4.csv'))
result_1_vs_7_tb <- read.csv(paste0(dataset, 'deseq_result_1_vs_7.csv'))


result_2_vs_3_tb <- read.csv(paste0(dataset, 'deseq_result_2_vs_3.csv'))
result_2_vs_6_tb <- read.csv(paste0(dataset, 'deseq_result_2_vs_6.csv'))
result_2_vs_4_tb <- read.csv(paste0(dataset, 'deseq_result_2_vs_4.csv'))
result_2_vs_7_tb <- read.csv(paste0(dataset, 'deseq_result_2_vs_7.csv'))

result_5_vs_3_tb <- read.csv(paste0(dataset, 'deseq_result_5_vs_3.csv'))
result_5_vs_6_tb <- read.csv(paste0(dataset, 'deseq_result_5_vs_6.csv'))
result_5_vs_4_tb <- read.csv(paste0(dataset, 'deseq_result_5_vs_4.csv'))
result_5_vs_7_tb <- read.csv(paste0(dataset, 'deseq_result_5_vs_7.csv'))

result_3_vs_4_tb <- read.csv(paste0(dataset, 'deseq_result_3_vs_4.csv'))
result_3_vs_7_tb <- read.csv(paste0(dataset, 'deseq_result_3_vs_7.csv'))
result_6_vs_4_tb <- read.csv(paste0(dataset, 'deseq_result_6_vs_4.csv'))
result_6_vs_7_tb <- read.csv(paste0(dataset, 'deseq_result_6_vs_7.csv'))



# make a plot - volcano plot


# subset <- result_padj %>%
#   filter(padj < 0.05, abs(log2FoldChange) > 3.0)

volcano_1_vs_3 <- ggplot(result_1_vs_3_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 1 vs 3') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_1_vs_3

volcano_1_vs_6 <- ggplot(result_1_vs_6_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 1 vs 6') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_1_vs_6

volcano_1_vs_4 <- ggplot(result_1_vs_4_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 1 vs 4') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_1_vs_4

volcano_1_vs_7 <- ggplot(result_1_vs_7_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 1 vs 7') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_1_vs_7

volcano_2_vs_3 <- ggplot(result_2_vs_3_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 2 vs 3') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_2_vs_3

volcano_2_vs_6 <- ggplot(result_2_vs_6_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 2 vs 6') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_2_vs_6

volcano_2_vs_4 <- ggplot(result_2_vs_4_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 2 vs 4') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_2_vs_4

volcano_2_vs_7 <- ggplot(result_2_vs_7_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 2 vs 7') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_2_vs_7

volcano_5_vs_3 <- ggplot(result_5_vs_3_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 5 vs 3') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_5_vs_3

volcano_5_vs_6 <- ggplot(result_5_vs_6_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 5 vs 6') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_5_vs_6

volcano_5_vs_4 <- ggplot(result_5_vs_4_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 5 vs 4') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_5_vs_4

volcano_5_vs_7 <- ggplot(result_5_vs_7_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 5 vs 7') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_5_vs_7

volcano_3_vs_4 <- ggplot(result_3_vs_4_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 3 vs 4') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_3_vs_4

volcano_3_vs_7 <- ggplot(result_3_vs_7_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 3 vs 7') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_3_vs_7

volcano_6_vs_4 <- ggplot(result_6_vs_4_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 6 vs 4') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_6_vs_4

volcano_6_vs_7 <- ggplot(result_6_vs_7_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 6 vs 7') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_6_vs_7

ggsave(volcano_1_vs_3, file = paste0(figure, 'rnaseq_volcano_1_vs_3.png'))
ggsave(volcano_1_vs_6, file = paste0(figure, 'rnaseq_volcano_1_vs_6.png'))
ggsave(volcano_1_vs_4, file = paste0(figure, 'rnaseq_volcano_1_vs_4.png'))
ggsave(volcano_1_vs_7, file = paste0(figure, 'rnaseq_volcano_1_vs_7.png'))

ggsave(volcano_2_vs_3, file = paste0(figure, 'rnaseq_volcano_2_vs_3.png'))
ggsave(volcano_2_vs_6, file = paste0(figure, 'rnaseq_volcano_2_vs_6.png'))
ggsave(volcano_2_vs_4, file = paste0(figure, 'rnaseq_volcano_2_vs_4.png'))
ggsave(volcano_2_vs_7, file = paste0(figure, 'rnaseq_volcano_2_vs_7.png'))

ggsave(volcano_5_vs_3, file = paste0(figure, 'rnaseq_volcano_5_vs_3.png'))
ggsave(volcano_5_vs_6, file = paste0(figure, 'rnaseq_volcano_5_vs_6.png'))
ggsave(volcano_5_vs_4, file = paste0(figure, 'rnaseq_volcano_5_vs_4.png'))
ggsave(volcano_5_vs_7, file = paste0(figure, 'rnaseq_volcano_5_vs_7.png'))

ggsave(volcano_3_vs_4, file = paste0(figure, 'rnaseq_volcano_3_vs_4.png'))
ggsave(volcano_3_vs_7, file = paste0(figure, 'rnaseq_volcano_3_vs_7.png'))

ggsave(volcano_6_vs_4, file = paste0(figure, 'rnaseq_volcano_6_vs_4.png'))
ggsave(volcano_6_vs_7, file = paste0(figure, 'rnaseq_volcano_6_vs_7.png'))

result_1_vs_2 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group2', 'Group1'))

result_3_vs_6 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group6', 'Group3'))

result_4_vs_7 <- results(deseq, lfcThreshold = 0.25, alpha = 0.05, 
                         contrast = c('condition', 'Group7', 'Group4'))

summary(result_1_vs_2)
summary(result_3_vs_6)
summary(result_4_vs_7)

result_1_vs_2_tb <- result_1_vs_2 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_3_vs_6_tb <- result_3_vs_6 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_4_vs_7_tb <- result_4_vs_7 %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

result_1_vs_2_tb$symbol <- mapIds(Mus.musculus,
                                  keys=result_tb$gene,
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")

result_3_vs_6_tb$symbol <- mapIds(Mus.musculus,
                                  keys=result_tb$gene,
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")

result_4_vs_7_tb$symbol <- mapIds(Mus.musculus,
                                  keys=result_tb$gene,
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")

write.csv(result_1_vs_2_tb, file = paste0(dataset, 'deseq_result_1_vs_2.csv'), row.names = FALSE)

write.csv(result_3_vs_6_tb, file = paste0(dataset, 'deseq_result_3_vs_6.csv'), row.names = FALSE)

write.csv(result_4_vs_7_tb, file = paste0(dataset, 'deseq_result_4_vs_7.csv'), row.names = FALSE)

result_1_vs_2_tb <- read.csv(paste0(dataset, 'deseq_result_1_vs_2.csv')) 
result_3_vs_6_tb <- read.csv(paste0(dataset, 'deseq_result_3_vs_6.csv')) 
result_4_vs_7_tb <- read.csv(paste0(dataset, 'deseq_result_4_vs_7.csv')) 

volcano_1_vs_2 <- ggplot(result_1_vs_2_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 1 vs 2') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_1_vs_2

volcano_3_vs_6 <- ggplot(result_3_vs_6_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 3 vs 6') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_3_vs_6

volcano_4_vs_7 <- ggplot(result_4_vs_7_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = c(-0.25, 0.25), alpha = 0.75,
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05)) +
  xlab('log2 Fold Change') +
  ylab('-log10 adjusted P value') +
  labs(color = 'Adjusted P value < 0.05') +
  ggtitle('group 4 vs 7') +
  # geom_text_repel(aes(label = symbol), data = subset, max.overlaps = 1000) +
  theme_bw()

# +
#   theme(legend.position = 'none')
volcano_4_vs_7

ggsave(volcano_1_vs_2, file = paste0(figure, 'rnaseq_volcano_1_vs_2.png'))
ggsave(volcano_3_vs_6, file = paste0(figure, 'rnaseq_volcano_3_vs_6.png'))
ggsave(volcano_4_vs_7, file = paste0(figure, 'rnaseq_volcano_4_vs_7.png'))
