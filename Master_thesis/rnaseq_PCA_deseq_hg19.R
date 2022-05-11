library(tidyverse)
library(tximport)
library(pheatmap)
library(ggrepel)
setwd('d:/Projects/UCPH/Work/rnaseq/')

my.files <- 
  c('LD001_salmon_hg19', 'LD002_salmon_hg19', 'LD005_salmon_hg19', 'LD006_salmon_hg19', 
    'LD007_salmon_hg19', 'LD008_salmon_hg19', 'LD009_salmon_hg19', 'LD010_salmon_hg19', 'LD011_salmon_hg19') %>%
  file.path(., 'quant.sf') %>%
  set_names(c('NT1', 'NT2', 'LNC1-1', 'LNC1-2', 'LNC2-1', 'LNC2-2', 'NT3', 'LNC1-3', 'LNC2-3'))
my.files

file.exists(my.files)

txi <- tximport(my.files, type = 'salmon', txOut = T)

knitr::kable(head(txi$abundance))


knitr::kable(head(txi$counts))

tpm_df <- as_data_frame(txi$abundance + 0.01)

ggplot(tpm_df, aes(x = NT1, y = NT2)) + geom_point(alpha = 0.05) + scale_x_log10() + scale_y_log10() + theme_bw()


pca <- prcomp(t(txi$abundance))

pca.x <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column('condition') %>%
  as_tibble()

pca.var <- summary(pca)$importance['Proportion of Variance',] * 100

pca_1 <- ggplot(pca.x, aes(x = PC1, y = PC2, color = condition, label = condition)) +
  geom_point() + 
  xlab(label = paste('PC1', pca.var[1])) +
  ylab(label = paste('PC2', pca.var[2])) +
  geom_label_repel() +
  theme_bw()
pca_1
# #ggsave(pca_1, file = 'pca_first.png')
# dist <- dist(t(txi$abundance), method = 'euclidian')
# tree <- hclust(dist, method = 'complete')
# plot(tree)
# 
# ggsave(plot(tree), file = 'hierachycal_plot_first.png')
# 
# abun <- txi$abundance %>%
#   as.data.frame() %>%
#   select(-c('LNC1-3', 'LNC2-3', 'NT3')) %>%
#   as.matrix()
# 
# pca_2 <- prcomp(t(abun))
# 
# pca.x2 <- pca_2$x %>%
#   as.data.frame() %>%
#   rownames_to_column('condition') %>%
#   as_tibble()
# 
# pca.var2 <- summary(pca)$importance['Proportion of Variance',] * 100
# 
# pca_2 <- ggplot(pca.x2, aes(x = PC1, y = PC2, color = condition, label = condition)) +
#   geom_point() + 
#   xlab(label = paste('PC1', pca.var2[1])) +
#   ylab(label = paste('PC2', pca.var2[2])) +
#   geom_label_repel() +
#   theme_bw()
# 
# ggsave(pca_2, file = 'pca_second.png')
# 
# dist2 <- dist(t(abun), method = 'euclidian')
# tree2 <- hclust(dist2, method = 'complete')
# plot(tree2)
# dist1 <- dist(txi$abundance, method = 'euclidian')
# tree1 <- hclust(dist1, method = 'complete')
# plot(tree1)
# 
# library(cluster)
# 
# a <- agnes(txi$abundance, method = 'average')
# 
# pheatmap(txi$abundance, cluster_rows = F, cluster_cols = F, show_rownames = T)

design <- data.frame(condition = c('NT', 'NT', 'LNC1', 'LNC1', 'LNC2', 'LNC2', 'NT', 'LNC1', 'LNC2'))
rownames(design) <- colnames(txi$abundance)

# design_2 <- design[1:6, ] %>%
#   as.data.frame()
# rownames(design_2) <- c('NT1', 'NT2', 'LNC1-1', 'LNC1-2', 'LNC2-1', 'LNC2-2')

count <- txi$counts

# # delete an outlier sample
# count_2 <- count %>%
#   as.data.frame() %>%
#   rownames_to_column('gene') %>%
#   as_tibble() %>%
#   select(-c(NT3, 'LNC1-3', 'LNC2-3'))
# 
# # convert to matrix
# count_2_mat <- as.matrix(count_2[2:7])
# rownames(count_2_mat) <- count_2$gene
# 
# # delete an outlier sample
# q3_normal_correct_tb <- q3_normal %>%
#   as.data.frame() %>%
#   rownames_to_column('gene') %>%
#   as_tibble() %>%
#   select(-Sample7)
# 
# # convert to matrix
# q3_normal_correct <- as.matrix(q3_normal_correct_tb[2:12])
# rownames(q3_normal_correct) <- q3_normal_correct_tb$gene
# 
# # correct sample annotations
# q3_design_correct <- q3_design[-7,]
# rownames(q3_design_correct) <- NULL
# 
# # confirm the result by hierarchical clustering
# q3_dist_confirm <- dist(t(q3_normal_correct), method = 'euclidian')
# q3_tree_confirm <- hclust(q3_dist_confirm, method = 'average')
# plot(q3_tree_confirm)



library(DESeq2)
deseq_matrix <- DESeqDataSetFromMatrix(countData = round(count), 
                                       colData = design, 
                                       design = ~condition)
deseq <- DESeq(deseq_matrix)

# sort result with logFC threshold 0.25, padj < 0.05
result <- results(deseq, lfcThreshold = 0.25, alpha = 0.05)

# confirm result
summary(result)

# make tibble from DESeq result
result_tb <- result %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  as_tibble()

# make a plot
ggplot(result_tb, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = padj<0.05), alpha = 0.5) +
  scale_x_log10() +
  geom_smooth() +
  geom_hline(yintercept = 0, alpha = 0.75, 
             color = 'darkgreen', linetype = 'dashed') +
  theme_bw()


ggplot(result_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.75, 
             color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = 1) +
  theme_bw()


# lowest padj
result_padj <- result_tb %>%
  arrange(padj) %>%
  filter(padj < 0.05)


# deseq_matrix2 <- DESeqDataSetFromMatrix(countData = round(count_2_mat), 
#                                         colData = design_2, 
#                                         design = ~condition)
# deseq2 <- DESeq(deseq_matrix2)
# 
# # sort result with logFC threshold 0.25, padj < 0.05
# result <- results(deseq, lfcThreshold = 0.25, alpha = 0.05)
# 
# # confirm result
# summary(result)
# 
# # make tibble from DESeq result
# result_tb <- result %>%
#   as.data.frame() %>%
#   rownames_to_column('gene') %>%
#   as_tibble()
# 
# # make a plot
# ggplot(result_tb, aes(x = baseMean, y = log2FoldChange)) +
#   geom_point(aes(color = padj<0.05), alpha = 0.5) +
#   scale_x_log10() +
#   geom_smooth() +
#   geom_hline(yintercept = 0, alpha = 0.75, 
#              color = 'darkgreen', linetype = 'dashed') +
#   theme_bw()
# 
# 
# ggplot(result_tb, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(color = padj < 0.05), alpha = 0.5) +
#   geom_vline(xintercept = 0, alpha = 0.75, 
#              color = 'black', linetype = 'dashed') +
#   geom_hline(yintercept = 1) +
#   theme_bw()
# 
# 
# # lowest padj
# result_padj <- result_tb %>%
#   arrange(padj) %>%
#   filter(padj < 0.05)
# 
# ggsave(result_padj, file = 'result_padj.png')
# 

library(limma)
rld <- rlog(deseq_matrix)
vsd <- vst(deseq_matrix, blind=FALSE)
mat <- assay(vsd)
mat <- removeBatchEffect(mat, vsd$batch)
assay(vsd) <- mat
plotPCA(vsd)

# + geom_text_repel(aes(label = colData(rld)$precise))
