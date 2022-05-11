#!/opt/local/stow/R-3.4.0/bin/Rscript

library(gplots)
#library(ggplot2)
library(reshape2)
library(pvclust)
library(RColorBrewer);
library(stringr)
library(tidyverse)
library(pheatmap)
setwd('d:/reusf2')

flagArrangeColumn = 0
Files = c("Gene.norm.txt")
motif <- read_tsv('d:/motif_usf2.txt')
colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
for (ff in 1:length(Files)){
	fname = Files[ff]
	E = read.delim(fname,row.names=1,header=T)

	mExp= E
	maxval = 2
	minval = -2
	my = mExp
	my = my[apply(my,1,function(x) sd(x)!=0),]
	my <- t(scale(t(my)))
	gene <- as_tibble(rownames(E))
	colnames(gene) <- 'genes'
	USF2 <- gene %>%
	  mutate(USF2 = ifelse(genes %in% motif$`Gene Name`, 'Included', 'Not Included'))
	  
	annotation <- data.frame(USF2 = factor(USF2$USF2), row.names = USF2$genes)
	

	cluster=vector(mode="character",length=nrow(mExp))

	d=as.dist(1-cor(t(my)))
	h=hclust(d, method="ward.D2")
	dend = as.dendrogram(h)
	#dend <- reorder(dend,1:4998)

	d2 =as.dist(1-cor((my)))
	h2 = hclust(d2,method="ward.D2")
	dend2 = as.dendrogram(h2)
	
	bk <- seq(-2, 2, by=0.1)
	data.mat=as.matrix(my)

	pngnamedist = paste("Results/",fname,".dist.png",sep="")
	png(pngnamedist,width=500,height=500)
	par(cex=1.2)
	plot(h2,main="transcripts",cex.main=1)
	dev.off()

	pngname=paste("Results/",fname,".Dendro.png",sep="") ; N.cluster=7; title="7"
	png(pngname,width=1200,height=1400)
	par(mar=c(13,1,1,1))

	rowColor=vector(mode="character",length=nrow(my))
	cluster=vector(mode="character",length=nrow(my))
	unlink("Results/Cluster*.txt")
	for(j in 1:N.cluster){
		sub.tf=cutree(h,k=N.cluster)==j
		rowColor[sub.tf]=colorL[j]
		clustername = paste("Results/Cluster",str_replace(fname,".txt",""),".",colorL[j],".",j,".txt",sep="")
		print(clustername)
		write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
		#write.table(names(sub.tf[sub.tf==TRUE]),clustername,col.names=F,row.names=F)
		cluster[sub.tf]=j
	}
	cluster <- as.matrix(cluster)
	rownames(cluster)<-rownames(my)
	imgDat = t(my[h$labels[h$order],])
	if (flagArrangeColumn ==1) imgDat = imgDat[h2$labels[h2$order],]	# cluster column as well
	imgDat[imgDat>maxval]=maxval
	imgDat[imgDat<minval]=minval

	mycol = colorL[1:N.cluster]
	cluster = cluster[h$labels[h$order],]
							    
	    heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
		 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, RowSideColors=rowColor,
		 cexCol=1.5,margins=c(15,5) )
	dev.off()

	png("Results/Gene.norm.Dendro2.png",width=1000,height=1400)
	    heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
		 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, #RowSideColors=F,
		 cexCol=1.5,margins=c(16,5) )
	dev.off()

}
print("Done")


