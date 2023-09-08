## Author: Mahdi

## Step1
# Set the working directory in the folder containing the cellranger files ("./")
setwd("./")
# load the cellranger matrix
library(cellrangerRkit)
genome<-"mm10.e87"
gbm <- load_cellranger_matrix("./", genome=genome)
# normalization
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
exp<-exprs(gbm_log)
expdf<-as.data.frame(as.matrix(exp))
texpdf<-t(expdf)
k<-kmeans(texpdf,center=9,nstart=10)
write.csv(k$cluster,"AllCells_KM9.csv")
write.csv(k$tot.withinss,"AllCells_KM9_totwithin.csv")

## Step2
# drawing heatmep for 9 clusters
# Drawing Heat map graph for marker gene experssion in different clusters
gene<-fData(gbm)
cellM<-read.csv("./data/CellMarkers.csv")

tg<-gene[gene$symbol %in% cellM$Gene,] #finding target gene ID
tge<-expdf[rownames(expdf) %in% tg$id,]
colnames(tg)[2]<-"Gene"
tg1<-merge(tg,cellM, by="Gene")
tg1<-tg1[order(tg1$CellType),]

tge<-tge[match(tg1$id,rownames(tge)),]
km9<-read.csv("AllCells_KM9.csv")
km9s<-km9[order(km9[,2]),]
tge<-tge[,match(km9s[,1],colnames(tge))]

anr<-data.frame(CellType=tg1$CellType)
rownames(anr)<-tg1$Gene
anc<-data.frame(Cluster=as.character(km9s$x))
rownames(anc)<-km9s$X
rownames(tge)<-tg1$Gene

maxtge<-apply(tge,1,max)
tge1<-tge
for (i in 1:nrow(tge1)){
  tge1[i,]<-tge[i,]/maxtge[i]
}

myColor <- colorRampPalette(c("white", "red"))(50)
library(pheatmap)
pheatmap(tge1,cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_row = anr,
         annotation_col = anc, annotation_names_row = TRUE, 
         annotation_names_col = TRUE,fontsize=6, color=myColor)

## Step3
# seperating the cells which are in the 7 (Enteroendocrine and Tuft)
cl7<-km9[km9$x==7,]
expdf7<-expdf[,colnames(expdf) %in% cl7$X]
k<-kmeans(t(expdf7),center=2,nstart=10)
write.csv(k$cluster,"cl7_KM2.csv")

# seperating tuft and Enteroendocrines 
cl7<-read.csv("cl7_KM2.csv")
km10<-km9
km10[km10$X %in% cl7[cl7$x==2,"X"],"x"]<-10

## Step4
#heat map for km10
tg<-gene[gene$symbol %in% cellM$Gene,] #finding target gene ID
tge<-expdf[rownames(expdf) %in% tg$id,]
colnames(tg)[2]<-"Gene"
tg1<-merge(tg,cellM, by="Gene")
tg1<-tg1[order(tg1$CellType),]

tge<-tge[match(tg1$id,rownames(tge)),]
km10s<-km10[order(km10[,2]),]
tge<-tge[,match(km10s[,1],colnames(tge))]

anr<-data.frame(CellType=tg1$CellType)
rownames(anr)<-tg1$Gene
anc<-data.frame(Cluster=as.character(km10s$x))
rownames(anc)<-km10s$X
rownames(tge)<-tg1$Gene

maxtge<-apply(tge,1,max)
tge1<-tge
for (i in 1:nrow(tge1)){
  tge1[i,]<-tge[i,]/maxtge[i]
}

myColor <- colorRampPalette(c("white", "red"))(50)
library(pheatmap)
pheatmap(tge1,cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, annotation_row = anr,
         annotation_col = anc, annotation_names_row = TRUE, 
         annotation_names_col = TRUE,fontsize=6, color=myColor)
write.csv(km10,"KM10.csv")

## Step5
# replacing the cluster numbers with names
colnames(km10)<-c("Index","Cluster")
km10n<-km10
name<-c("Early Enterocytes Progenitors","Secretory Progenitors","Enterocyte",
        "Goblet","Stem","Paneth","Enteroendocrine","TA","Late Enterocytes Progenitors","Tuft")
km10n$Cluster<-name[km10$Cluster]
write.csv(km10n,"Cluster10_withname.csv")

# PCA
t.pca <- prcomp(t(expdf),center = TRUE,scale. = TRUE) 
save.image("pca.RData")

library(ggbiplot)
g <- ggbiplot(t.pca, groups = as.character(km10n$Cluster),obs.scale = 1, var.scale = 1, circle = TRUE,var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

# Finding number of young and old cells in each cluster
sam<-substr(as.character(km10n$Index),18,18)
g<-c("Old","Young")
sam<-g[as.numeric(sam)]
km10<-cbind(km10n,Group=sam)
countR<-data.frame(ClusterName=NULL,Group=NULL,Percent=NULL)
yt<-nrow(km10[km10$Group=="Young",])
ot<-nrow(km10[km10$Group=="Old",])
km10$Cluster<-as.factor(km10$Cluster)
for (i in levels(km10$Cluster)){
  yn<-nrow(km10[km10$Cluster==i & km10$Group=="Young",])
  on<-nrow(km10[km10$Cluster==i & km10$Group=="Old",])
  
  countR<-rbind(countR,data.frame(ClusterName=i,Group="Young",Percent=100*yn/yt))
  countR<-rbind(countR,data.frame(ClusterName=i,Group="Old",Percent=100*on/ot))
}
write.csv(countR,"NumberOfcellsinClusters10.csv")

ggplot(countR,aes(x=as.character(countR$ClusterName),y=Percent))+
  geom_bar(aes(fill = Group), position = "dodge", stat="identity")+
  xlab("Cluster")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))

# fold change in young and old
cY<-countR[countR$Group=="Young",]
cO<-countR[countR$Group=="Old",]
lfc<-log2(cO$Percent/cY$Percent)
names(lfc)<-cY$ClusterName
barplot(lfc,las=2,col=rainbow(10),ylab="Log2 fold change")

# Finding top markers
library(cellrangerRkit)
genome<-"mm10.e87"
gbm <- load_cellranger_matrix("./", genome=genome)
cluster_result <-read.csv("KM10.csv")[,c(2,3)]
colnames(cluster_result)<-c("Barcode","Cluster")

# sort the cells by the cluster labels
cells_to_plot <- order_cell_by_clusters(gbm, cluster_result$Cluster)
# order the genes from most up-regulated to most down-regulated in each cluster
prioritized_genes <- prioritize_top_genes(gbm, cluster_result$Cluster, "sseq", min_mean=0.1,logscale =FALSE)
save.image("topgene.RData")
for(i in 1:10){
  tmp<-as.data.frame(prioritized_genes[i])
  tmp<-tmp[tmp[,10]>0,]
  write.csv(tmp,paste("Marker_all/TopGene_Cluster",i,".csv",sep=""))
}
dir.create("gene_set")
gene<-write_cluster_specific_genes(prioritized_genes, "gene_set", n_genes=70)
gbm_pheatmap(log_gene_bc_matrix(gbm), prioritized_genes, cells_to_plot,
             n_genes=10, limits=c(0,2),colour = rainbow(10))
#finding best markers for each cluster
allg<-data.frame()
for (i in 1:10){
  p<-prioritized_genes[[i]]
  p<-p[p$significant,]
  if(nrow(p)==0){next}
  allg<-rbind(allg,data.frame(Cluster=i,GeneSymbol=p$symbol,Ensemble=p$id,logfc=p$log2fc,p_adj=p$p_adj))
}

table(allg$Cluster)
write.csv(allg,"AllClusters_topGene_logfc1_qvalue0.05mean0.1.csv")
name<-c("Early Enterocytes Progenitors","Secretory Progenitors","Enterocyte",
        "Goblet","Stem","Paneth","Enteroendocrine","TA","Late Enterocytes Progenitors","Tuft")
allg$Cluster<-name[allg$Cluster]
write.csv(allg,"AllClusters_withclustername_topGene_logfc1_qvalue0.05mean0.1.csv")

#drawing heatmap for my top markers
gene<-fData(gbm)
km10n<-read.csv("Cluster10_withname.csv",row.names=1)
allg<-read.csv("AllClusters_withclustername_topGene_logfc1_qvalue0.05.csv",row.names=1)
cellM<-allg[,c(1,2)]
colnames(cellM)<-c("CellType",	"Gene")
tg<-gene[gene$symbol %in% cellM$Gene,] #finding target gene ID
tge<-expdf[rownames(expdf) %in% tg$id,]
colnames(tg)[2]<-"Gene"
tg1<-merge(tg,cellM, by="Gene")
tg1<-tg1[order(tg1$CellType),]

tge<-tge[match(tg1$id,rownames(tge)),]
km10ns<-km10n[order(km10n[,2]),]
tge<-tge[,match(km10ns[,1],colnames(tge))]

anr<-data.frame(CellType=tg1$CellType)
rownames(anr)<-make.names(tg1$Gene,unique=TRUE)
anc<-data.frame(Cluster=as.character(km10ns$Cluster))
rownames(anc)<-km10ns$Gene
rownames(tge)<-make.names(tg1$Gene,unique=TRUE)

tge1<-t(apply(tge, 1, function(x)(x-min(x))/(max(x)-min(x))))

myColor <- colorRampPalette(c("white", "red"))(50)
library(pheatmap)
pheatmap(tge1,cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE,show_rownames = FALSE, annotation_row = anr,
         annotation_col = anc, annotation_names_row = TRUE, 
         annotation_names_col = TRUE,fontsize=6, color=myColor)

## Step6
#Distribution of cells in PC1 and PC2
load("pca.RData")
library(cellrangerRkit)
library(sm)
genome<-"mm10.e87"
gbm <- load_cellranger_matrix("./", genome=genome)
## normalization
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
exp<-exprs(gbm_log)
expdf<-as.data.frame(as.matrix(exp))
gene<-fData(gbm)

pdf("PC1AndPC2.pdf")
km10<-read.csv("Cluster10_withname_Age.csv",row.names=1)
pcasel<-data.frame(t.pca$x[,1:2])

pcasel<-cbind(pcasel,Cluster=km10$Cluster,Age=km10$Age)
pcasel<-pcasel[pcasel$Cluster %in% c("Enteroendocrine","Secretory Progenitors","Tuft","Goblet","Paneth","Stem"),]
library(sm)
sm.density.compare(pcasel$PC1,pcasel$Age,col=c("red","Blue"))
legend("topright", levels(pcasel$Age), fill=c("red","Blue"))

pcasel$Cluster<-droplevels(pcasel$Cluster)
# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 8, 4))
sm.density.compare(pcasel$PC1,pcasel$Cluster,col=rainbow(6))
legend("top", cex=0.5,levels(pcasel$Cluster), fill=rainbow(6))

# with average of scaled markers
mar<-read.csv("./AllClusters_withclustername_topGene_logfc1_qvalue0.05mean0.1.csv",stringsAsFactors = FALSE)
mar<-mar[mar$Cluster %in% c("Enteroendocrine","Tuft","Goblet","Paneth"),]
km<-km10[km10$Cluster %in% c("Enteroendocrine","Secretory Progenitors","Tuft","Goblet","Paneth","Stem"),]
expdfsel<-expdf[,colnames(expdf) %in% rownames(km)]
expmar<-expdfsel[rownames(expdfsel) %in% mar$Ensemble,]
expmar<-t(scale(t(expmar)))
expmar<-apply(expmar,2,mean)
res<-data.frame(PC1=pcasel$PC1,SMarker=expmar)
library(ggplot2)
ggplot(res,aes(x=PC1,y=SMarker))+
  geom_smooth()+
  ylab("Secretory Markers")

#For stem markers
mar<-read.csv("AllClusters_withclustername_topGene_logfc1_qvalue0.05mean0.1.csv",stringsAsFactors = FALSE)
mar<-mar[mar$Cluster =="Stem",]
expmar<-expdfsel[rownames(expdfsel) %in% mar$Ensemble,]
expmar<-t(scale(t(expmar)))
expmar<-apply(expmar,2,mean)
res<-data.frame(PC1=pcasel$PC1,SMarker=expmar)
library(ggplot2)
ggplot(res,aes(x=PC1,y=expmar))+
  geom_smooth()+
  ylab("Stem Markers")



#for PC2
pcasel<-data.frame(t.pca$x[,1:2])

pcasel<-cbind(pcasel,Cluster=km10$Cluster,Age=km10$Age)
pcasel<-pcasel[pcasel$Cluster %in% c("TA","Early Enterocytes Progenitors","Late Enterocytes Progenitors","Enterocyte","Stem"),]
library(sm)
sm.density.compare(pcasel$PC2,pcasel$Age,col=c("red","Blue"))
legend("topright", levels(pcasel$Age), fill=c("red","Blue"))

pcasel$Cluster<-droplevels(pcasel$Cluster)
# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 8, 4))
sm.density.compare(pcasel$PC2,pcasel$Cluster,col=rainbow(5))

legend("topright", cex=0.5,levels(pcasel$Cluster), fill=rainbow(5))

#with average of scaled markers
mar<-read.csv("AllClusters_withclustername_topGene_logfc1_qvalue0.05mean0.1.csv",stringsAsFactors = FALSE)
mar<-mar[mar$Cluster %in% c("Enterocyte"),]
km<-km10[km10$Cluster %in% c("TA","Early Enterocytes Progenitors","Late Enterocytes Progenitors","Enterocyte","Stem"),]
expdfsel<-expdf[,colnames(expdf) %in% rownames(km)]
expmar<-expdfsel[rownames(expdfsel) %in% mar$Ensemble,]
expmar<-t(scale(t(expmar)))
expmar<-apply(expmar,2,mean)
res<-data.frame(PC2=pcasel$PC2,SMarker=expmar)
library(ggplot2)
ggplot(res,aes(x=PC2,y=SMarker))+
  geom_smooth()+
  ylab("Entrocyte Markers")

#For stem markers
mar<-read.csv("AllClusters_withclustername_topGene_logfc1_qvalue0.05mean0.1.csv",stringsAsFactors = FALSE)
mar<-mar[mar$Cluster =="Stem",]
expmar<-expdfsel[rownames(expdfsel) %in% mar$Ensemble,]
expmar<-t(scale(t(expmar)))
expmar<-apply(expmar,2,mean)
res<-data.frame(PC2=pcasel$PC2,SMarker=expmar)
library(ggplot2)
ggplot(res,aes(x=PC2,y=expmar))+
  geom_smooth()+
  ylab("Stem Markers")
dev.off()
