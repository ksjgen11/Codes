library(tidyverse)
library(pheatmap)
setwd('d:')

C1 <- read_tsv('reusf2/ClusterSelectGene.red.1.txt', col_names = 'gene') %>%
  mutate(cluster = 'Cluster 1')
C2 <- read_tsv('reusf2/ClusterSelectGene.purple.2.txt', col_names = 'gene') %>%
  mutate(cluster = 'Cluster 2')
C3 <- read_tsv('reusf2/ClusterSelectGene.blue.3.txt', col_names = 'gene') %>%
  mutate(cluster = 'Cluster 3')
C4 <- read_tsv('reusf2/ClusterSelectGene.yellow.4.txt', col_names = 'gene') %>%
  mutate(cluster = 'Cluster 4')
C5 <- read_tsv('reusf2/ClusterSelectGene.green.5.txt', col_names = 'gene') %>%
  mutate(cluster = 'Cluster 5')
C6 <- read_tsv('reusf2/ClusterSelectGene.orange.6.txt', col_names = 'gene') %>%
  mutate(cluster = 'Cluster 6')
clusters <- bind_rows(C1, C2, C3, C4, C5, C6)

motif <- read_tsv('motif_usf2.txt')



gene.norm <- read.table('reusf2/Gene.norm.txt') %>%
  rownames_to_column('gene')
  


gene.norm.m <- as.matrix(gene.norm[2:9])
rownames(gene.norm.m) <- gene.norm$gene




genes<- gene.norm$gene
annotation <- clusters %>%
  mutate(USF2 = ifelse(gene %in% motif$`Gene Name`, 'Included', 'Not Included'))

annotation2 <- data.frame(cluster = factor(annotation$cluster),
                          USF2 = factor(annotation$USF2))
rownames(annotation2) <- annotation$gene

pheatmap(gene.norm.m, annotation_row = annotation2, cutree_cols = 8, cutree_rows = 6,
         border_color = 'white')

heatmap <- pheatmap(gene.norm.m, annotation_row = annotation2, cutree_cols = 8, cutree_rows = 6, 
                    border_color = 'white', show_rownames = F)

save(gene.norm, gene.norm.m, motif, clusters, annotation, annotation2, 
     heatmap, file = 'USF_heatmap.RData')
