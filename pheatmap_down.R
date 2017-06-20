library(pheatmap)
library(RColorBrewer)
library(dendsort)

down<-read.table(file="aa",row.names=1,header=T,fill=T,quote="",sep = "\t")
row_annotation<-data.frame(annotation=down$annotation)
rownames(row_annotation) = rownames(down)
down<-down[,43:66]
col1=colorRampPalette(c("blue","white","red"))(100)

set.seed(4)
annotation_col = data.frame(
	condition = factor(rep(c("control", "treatment"), 4, each=3)),
	phenotype = c(rep("WT", 6), rep("mutant", 18))
	)
rownames(annotation_col) = colnames(down)
down_p<-pheatmap(log2(down+1),border=F,col=col1,cluster_rows=T,cluster_cols=T,scale="row",legend=T,show_colnames=T,show_rownames=F,annotation_col = annotation_col,annotation_row=row_annotation,cutree_rows=5,annotation_colors=list(annotation=c("jasmonic acid metabolic process"="red","response to alcohol"="blue","response to wounding"="green","no_annotation"="grey","response to hormone stimulus"="purple")),cellwidth=20)
down_cluster<-cbind(down, cluster = cutree(down_p$tree_row, k = 5))
row_annotation=cbind(row_annotation,cluster = cutree(down_p$tree_row, k = 5))
row_annotation$cluster=as.character(row_annotation$cluster)
pheatmap(log2(down+1),border=F,col=col1,cluster_rows=T,cluster_cols=T,scale="row",legend=T,show_colnames=T,show_rownames=F,annotation_col = annotation_col,annotation_row=row_annotation,cutree_rows=5,annotation_colors=list(annotation=c("jasmonic acid metabolic process"="red","response to alcohol"="blue","response to wounding"="green","no_annotation"="grey","response to hormone stimulus"="purple"),cluster=c("1"="red","2"="green","3"="yellow","4"="black","5"="orange")),filename="pheatmap_down.pdf",cellwidth=20)
write.table(down_cluster[down_cluster$cluster==1,],file="down_cluster1",quote=FALSE,sep="\t")
write.table(down_cluster[down_cluster$cluster==2,],file="down_cluster2",quote=FALSE,sep="\t")
write.table(down_cluster[down_cluster$cluster==3,],file="down_cluster3",quote=FALSE,sep="\t")
write.table(down_cluster[down_cluster$cluster==4,],file="down_cluster4",quote=FALSE,sep="\t")
write.table(down_cluster[down_cluster$cluster==5,],file="down_cluster5",quote=FALSE,sep="\t")
