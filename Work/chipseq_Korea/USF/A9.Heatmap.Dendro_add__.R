#!/opt/local/stow/R-3.4.0/bin/Rscript

library(gplots)
#library(ggplot2)
library(reshape2)
library(pvclust)
library(RColorBrewer);
library(stringr)
library(tidyverse)
library(pheatmap)


flagArrangeColumn = 0


Files = c("Gene.norm.txt")

colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
for (ff in 1:length(Files)){
	fname = Files[ff]
	E = read.delim(fname,row.names=1,header=T) %>%
	  t() %>%
	  as.data.frame() %>%
	  rownames_to_column('sample') %>%
	  arrange(sample) 
	
	row_E <- E$sample
	E2 <- E %>%
	  select(-sample)
	rownames(E2) <- row_E
	
	E3 <- E2 %>%
	  t() %>%
	  as.matrix()
	

	mExp= E3
	maxval = 2
	minval = -2
	my = mExp
	my = my[apply(my,1,function(x) sd(x)!=0),]
	my <- t(scale(t(my)))
	
	

# 	  cluster=vector(mode="character",length=nrow(mExp))
#
# 	d=as.dist(1-cor(t(my)))
# 	h=hclust(d, method="ward.D2")
# 	dend = as.dendrogram(h)
# 	#dend <- reorder(dend,1:4998)
#
# 	d2 =as.dist(1-cor((my)))
# 	h2 = hclust(d2,method="ward.D2")
# 	dend2 = as.dendrogram(h2)
#
  	bk <- seq(-2, 2, by=0.1)
# 	data.mat=as.matrix(my)
#
# 	pngnamedist = paste("Results/",fname,".dist.png",sep="")
# 	png(pngnamedist,width=500,height=500)
# 	par(cex=1.2)
# 	plot(h2,main="transcripts",cex.main=1)
# 	dev.off()
#
# 	N.cluster=7; title="7"
#
# 	rowColor=vector(mode="character",length=nrow(my))
# 	cluster=vector(mode="character",length=nrow(my))
# 	unlink("Results/Cluster*.txt")
# 	for(j in 1:N.cluster){
# 		sub.tf=cutree(h,k=N.cluster)==j
# 		rowColor[sub.tf]=colorL[j]
# 		clustername = paste("Results/Cluster",str_replace(fname,".txt",""),".",colorL[j],".",j,".txt",sep="")
# 		print(clustername)
# 		write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
# 		#write.table(names(sub.tf[sub.tf==TRUE]),clustername,col.names=F,row.names=F)
# 		cluster[sub.tf]=j
# }
	
	C1 <- read_tsv('Results/ClusterGene.norm.red.1.txt', col_names = 'gene') %>%
	  mutate(cluster = 'Cluster 1')
	C2 <- read_tsv('Results/ClusterGene.norm.purple.2.txt', col_names = 'gene') %>%
	  mutate(cluster = 'Cluster 2')
	C3 <- read_tsv('Results/ClusterGene.norm.blue.3.txt', col_names = 'gene') %>%
	  mutate(cluster = 'Cluster 3')
	C4 <- read_tsv('Results/ClusterGene.norm.yellow.4.txt', col_names = 'gene') %>%
	  mutate(cluster = 'Cluster 4')
	C5 <- read_tsv('Results/ClusterGene.norm.green.5.txt', col_names = 'gene') %>%
	  mutate(cluster = 'Cluster 5')
	C6 <- read_tsv('Results/ClusterGene.norm.orange.6.txt', col_names = 'gene') %>%
	  mutate(cluster = 'Cluster 6')
	C7 <- read_tsv('Results/ClusterGene.norm.brown.7.txt', col_names = 'gene') %>%
	  mutate(cluster = 'Cluster 7')
	clustering <- bind_rows(C1, C2, C3, C4, C5, C6, C7)
	clustering2 <- data.frame(cluster = clustering$cluster, row.names = clustering$gene)
	motif <- read_tsv('USF_motif_usf2_500_100.txt')
	motif2 <- motif %>%
	  filter(!is.na(`Usf2(bHLH)/C2C12-Usf2-ChIP-Seq(GSE36030)/Homer Distance From Peak(sequence,strand,conservation)`))
	
	gene <- as_tibble(rownames(clustering2))
	colnames(gene) <- 'genes'
	USF2 <- gene %>%
	  mutate(USF2 = ifelse(genes %in% motif2$`Gene Name`, 'Included', 'Not_Included'))
	
	annotation <- data.frame(Cluster = clustering$cluster, 
	                         USF2 = factor(USF2$USF2), row.names = USF2$genes)
	
	gene.norm <- my %>%
	  as.data.frame() %>%
	  rownames_to_column('gene') %>%
	  mutate(cluster=ifelse(gene %in% C1$gene, 'Cluster 1', 
	                        ifelse(gene %in% C2$gene, 'Cluster 2', 
	                               ifelse(gene %in% C3$gene, 'Cluster 3', 
	                                      ifelse(gene %in% C4$gene, 'Cluster 4', 
	                                             ifelse(gene %in% C5$gene, 'Cluster 5', 
	                                                    ifelse(gene %in% C6$gene, 'Cluster 6', 'Cluster 7'))))))) %>%
	  arrange(cluster)
	

	
	my2 <- as.matrix(gene.norm[2:9])
	rownames(my2) <- gene.norm$gene
	annotation_colors = list(
	  Cluster = c('Cluster 1' = "red",
	              'Cluster 2' = "purple",
	              'Cluster 3' ="blue",
	              'Cluster 4' ="yellow",
	              'Cluster 5' ="green",
	              'Cluster 6' ="orange",
	              'Cluster 7' ="brown"),
	  USF2 = c(Included="black", Not_Included="white"))
	
	png(filename = 'heatmap_100bps_reordered_500_100.png', width = 1100, height = 1400)
	pheatmap(as.matrix(my2), cutree_row = 7, show_rownames = F, annotation_row = annotation, 
	         cluster_cols = F, cluster_rows = F, annotation_colors = annotation_colors,  
	         color=colorpanel(length(bk)-1,"blue","white","red"), fontsize = 20, scale = 'row')
	dev.off()
	
	
	includegenes <- annotation %>%
	  filter(USF2 == 'Included', Cluster == 'Cluster 1') %>%
	  write.table(file = 'usf2_cluster1_500_100.txt', quote=F,col.names=T,row.names=T)
	
	includegenes2 <- annotation %>%
	  filter(USF2 == 'Included', Cluster == 'Cluster 4') %>%
	  write.table(file = 'usf2_cluster4_500_100.txt', quote=F,col.names=T,row.names=T)
	  
	# cluster <- as.matrix(cluster)
	# rownames(cluster)<-rownames(my)
	# imgDat = t(my[h$labels[h$order],])
	# if (flagArrangeColumn ==1) imgDat = imgDat[h2$labels[h2$order],]	# cluster column as well
	# imgDat[imgDat>maxval]=maxval
	# imgDat[imgDat<minval]=minval
	# 
	# mycol = colorL[1:N.cluster]
	# cluster = cluster[h$labels[h$order],]
	# 
	# pngname=paste("Results/",fname,".Dendro.png",sep="") 
	# 
	# png(pngname,width=1200,height=1400)
	# 
	# par(mar=c(13,1,1,1))						    
	# 
	#     heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
	# 	 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, RowSideColors=rowColor,
	# 	 cexCol=1.5,margins=c(15,5) )
	# 
	# dev.off()
	# 
	# png("Results/Gene.norm.Dendro2.png",width=1000,height=1400)
	#     heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
	# 	 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, #RowSideColors=F,
	# 	 cexCol=1.5,margins=c(16,5) )
	# dev.off()

}
print("Done")


