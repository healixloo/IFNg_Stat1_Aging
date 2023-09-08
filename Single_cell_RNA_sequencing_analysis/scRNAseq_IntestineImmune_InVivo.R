## Author: Jing

setwd("./")
myData <- readRDS(file = "./data/immuneF_single_F3.rds")
# subset small intestine proximal
myData_SIP<-subset(myData, subset = pos == "SI_P")
unique(myData_SIP$pos)
myData_SIP$seurat_clusters<-NULL
myData_SIP$RNA_snn_res.0.7<-NULL
myData_SIP$RNA_snn_res.0.3<-NULL
myData_SIP$clusters_age<-NULL
# clustering
library(Seurat)
myData_SIP <- NormalizeData(object = myData_SIP, normalization.method = "LogNormalize", scale.factor = 1e4)
myData_SIP <- FindVariableFeatures(object = myData_SIP,selection.method = 'mvp',mean.cutoff = c(0.0125, 2), dispersion.cutoff = c(0.9, Inf))
all.genes <- rownames(x = myData_SIP)
myData_SIP <- ScaleData(object = myData_SIP)
myData_SIP <- RunPCA(object = myData_SIP, features = VariableFeatures(object = myData_SIP))
DimPlot(object = myData_SIP, reduction = 'pca')
myData_SIP <- FindNeighbors(object = myData_SIP, dims = 1:20)
# res=0.3
myData_SIP <- FindClusters(object = myData_SIP, resolution = 0.3)
myData_SIP <- RunTSNE(object = myData_SIP, dims = 1:20)
DimPlot(object = myData_SIP, reduction = 'tsne')
# res=0.7
myData_SIP <- FindClusters(object = myData_SIP, resolution = 0.7)
myData_SIP <- RunTSNE(object = myData_SIP, dims = 1:20)
DimPlot(object = myData_SIP, reduction = 'tsne')
# res=0.4
myData_SIP <- FindClusters(object = myData_SIP, resolution = 0.4)
myData_SIP <- RunTSNE(object = myData_SIP, dims = 1:20)
DimPlot(object = myData_SIP, reduction = 'tsne')
DimPlot(object = myData_SIP, reduction = 'tsne',label = T)
DimPlot(object = myData_SIP, reduction = 'tsne',split.by = "age")+coord_fixed(ratio = 1)
# Calculate markers
myData_SIP.markers <- FindAllMarkers(object = myData_SIP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(myData_SIP.markers,"./data/myData_SIP.markers.csv")
# compare with previous markers top20
markers_all<-read.csv("./data/immuneF_single_F3.markers.csv",header = T)
library(dplyr)
myData_SIP.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) -> top20_target
markers_all %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) -> top20_reference
top20_target$cluster<-paste("t",top20_target$cluster,sep="_")
top20_reference$cluster<-paste("r",top20_reference$cluster,sep="_")
clusters_celltypes_top20<-merge(top20_target,top20_reference,by.x='gene',by.y='gene',all=T)
library(RColorBrewer)
nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
ggplot(clusters_celltypes_top20, aes(cluster.x))+geom_bar(aes(fill=cluster.y))+scale_fill_manual(values = mycolors)
# compare with previous markers top50
myData_SIP.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50_target
markers_all %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC) -> top50_reference
top50_target$cluster<-paste("t",top50_target$cluster,sep="_")
top50_reference$cluster<-paste("r",top50_reference$cluster,sep="_")
clusters_celltypes_top50<-merge(top50_target,top50_reference,by.x='gene',by.y='gene',all=T)
ggplot(clusters_celltypes_top50, aes(cluster.x))+geom_bar(aes(fill=cluster.y))+scale_fill_manual(values = mycolors)
write.csv(top20_target,"./data/markers_top20.csv")
write.csv(top50_target,"./data/markers_top50.csv")

## assign celltypes
celltypes<-c("B-cell","CD4-naive","UNKNOWN-Th17_CD4-TEM_CD4-Gata3","CD4-Cytotoxic","Plasma-cell","Monocyte-Neutrophils-Macrophages","B-early","CD4-Gata3","CD8_NK","UNKNOWN")
# check
head(top20_target[top20_target$cluster=="t_9",])
celltypes[10]

myData_SIP$celltypes_new<-celltypes[myData_SIP$seurat_clusters]
myData_SIP$celltypes_old<-myData_SIP$celltypes
myData_SIP$celltypes<-NULL
# transform seurat object to normal dataframe to draw ggplot
my.data<-FetchData(myData_SIP,c("orig.ident","sampleID","Index","nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters","type","TYPE","TYPES","age","pos","celltypes_old","celltypes_new"))
pca<-Embeddings(myData_SIP, reduction = "pca")[, 1:20]
tsne<-Embeddings(myData_SIP, reduction = "tsne")[, 1:2]
my.data<-merge(my.data,tsne,by="row.names",all.x=T)
rownames(my.data)<-my.data$Row.names
my.data<-my.data[,-1]
my.data<-merge(my.data,pca,by="row.names",all.x=T)
rownames(my.data)<-my.data$Row.names
my.data<-my.data[,-1]
write.csv(my.data,"./data/my.data.csv")

#ggplot(my.data, aes(celltypes_new))+geom_bar(aes(fill=celltypes_old))+scale_fill_manual(values = mycolors)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=10),legend.text = element_text(size = 10),legend.title = element_text())
ggplot(my.data, aes(celltypes_new,fill=celltypes_old))+geom_bar()+scale_fill_manual(values = mycolors)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=10),legend.text = element_text(size = 10),legend.title = element_text())

#my.data<-read.csv("./data/my.data.csv",header = T)
my.data$celltypes_new_2<-as.character(my.data$celltypes_new)
my.data[my.data$celltypes_new_2=="UNKNOWN","celltypes_new_2"]<-"B-cell"
my.data[my.data$celltypes_new_2=="UNKNOWN-Th17_CD4-TEM_CD4-Gata3","celltypes_new_2"]<-"Th17"
my.data[my.data$celltypes_new_2=="Monocyte-Neutrophils-Macrophages","celltypes_new_2"]<-"Myeloblast"

library(stringr)
library(Rmisc)

my.data2<-my.data[,c("celltypes_new_2", "age","Index")]
my.data2$var<-paste(my.data2$celltypes_new_2,my.data2$age,my.data2$Index,sep="#")
my.data2<-as.data.frame(table(my.data2$var))
my.data2$celltypes_new_2<-str_split_fixed(my.data2$Var1,"#",3)[,1]
my.data2$age<-str_split_fixed(my.data2$Var1,"#",3)[,2]
my.data2$Index<-str_split_fixed(my.data2$Var1,"#",3)[,3]

# calculate celltype ratio in each sample according to celltype counts in my.data2

library(dplyr)
my.data2<-as.data.frame(my.data2)
# group_by will be useless /invalid if Rmisc(use plyr as default) was load, so deteach it before loading dplyr, or assign delyr in font of function 
detach(package:Rmisc)
detach(package:plyr)
#tg<-my.data2 %>% group_by(age,Index)  %>% mutate(p = Freq/sum(Freq))
tg<-my.data2 %>% dplyr::group_by(age,Index)  %>% dplyr::mutate(ratio = Freq/sum(Freq))
tg$Var2<-paste(tg$Index,tg$age,sep="#")
# check 
sum(tg[tg$Var2=="5#old","Freq"])
tg[tg$Var2=="5#old",]
121/674

# /method1: calculate sd/se for ggplot
library(Rmisc)
tgc <- summarySE(tg, measurevar="ratio", groupvars=c("celltypes_new_2","age"))
# ggplot sd
ggplot(tgc, aes(x=celltypes_new_2, y=ratio, fill=age)) +
geom_bar(position=position_dodge(), stat="identity") +
geom_errorbar(aes(ymin=ratio-sd, ymax=ratio+sd),
width=.2,                    # Width of the error bars
position=position_dodge(.9))
# ggplot se
#tgc <- summarySE(tg, measurevar="ratio", groupvars=c("celltypes_new_2","age"))
tgc$age<-factor(tgc$age,levels = c("young","old"))
ggplot(tgc, aes(x=celltypes_new_2, y=ratio, fill=age)) +
geom_bar(position=position_dodge(), stat="identity") +
geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
width=.2,                    # Width of the error bars
position=position_dodge(.9))

# /method2: ggpubr for sd/se and p value
library(ggpubr)
tg$age<-factor(tg$age,levels = c("young","old"))

ggbarplot(tg, x = "celltypes_new_2", y = "ratio", add = "mean_se",
color = "age", palette = "jco",
position = position_dodge(0.8))+
stat_compare_means(aes(group = age), label = "p.signif", label.y = 0.75)

ggbarplot(tg, x = "celltypes_new_2", y = "ratio", add = "mean_se",
color = "age", palette = "jco",
position = position_dodge(0.8))+
stat_compare_means(aes(group = age), label = "p.signif")

ggbarplot(tg, x = "celltypes_new_2", y = "ratio", add = "mean_se",
color = "age", palette = "jco",
position = position_dodge(0.8))+
stat_compare_means(aes(group = age),method = "wilcox.test", label = "p.signif", label.y =  c(0.75, 0.19, 0.2,0.15,0.2,0.15,0.2,0.2,0.2))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggbarplot(tg, x = "celltypes_new_2", y = "ratio", add = "mean_se",
color = "age", palette = "jco",
position = position_dodge(0.8))+
stat_compare_means(aes(group = age),method = "t.test", label = "p.signif", label.y =  c(0.75, 0.19, 0.2,0.15,0.2,0.15,0.2,0.2,0.2))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


## 2
## SingleImmune_ggbarplot_se_p_ttest_2.pdf
p<-ggbarplot(tg, x = "celltypes_new_2", y = "ratio", add = "mean_se",
           color = "age", palette = "jco",
           position = position_dodge(0.8))
p+stat_compare_means(aes(group = age),method = "t.test", label = "p.format", label.y =  c(0.75, 0.19, 0.2,0.15,0.2,0.15,0.2,0.2,0.2))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ggtitle("t.test_pairedF")+scale_color_manual(values = c("grey","grey0"))

p<-ggbarplot(tg, x = "celltypes_new_2", y = "ratio", add = "mean_sd",
           color = "age", palette = "jco",
           position = position_dodge(0.8))
p+stat_compare_means(aes(group = age),method = "t.test", label = "p.format", label.y =  c(0.75, 0.19, 0.2,0.15,0.2,0.15,0.2,0.2,0.2))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ggtitle("t.test_pairedF")+scale_color_manual(values = c("grey","grey0"))

write.csv(tg,"./data/tg.csv")
write.csv(tgc,"./data/tgc.csv")
write.csv(my.data,"./data/my.data.csv")
write.csv(my.data2,"./data/my.data2.csv")
saveRDS(myData_SIP,"./data/myData_SIP.RDS")

myData_SIP$celltypes_new_2<-(my.data[match(myData_SIP$sampleID,my.data$sampleID),]$celltypes_new_2)
install.packages("patchwork")
# display cytokine expression in all the celltypes

plots <- VlnPlot(myData_SIP, features = c("Ifng", "Il4", "Il13"), split.by = "age", group.by = "celltypes_new_2", pt.size = 1, combine = FALSE)
patchwork::wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(myData_SIP, features = c("Ifng", "Il4", "Il13"), split.by = "age", group.by = "celltypes_new_2", pt.size = 1, combine = FALSE,log=TRUE)
patchwork::wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(myData_SIP, features = c("Ifng", "Il4", "Il13"),  group.by = "celltypes_new_2", pt.size = 1, combine = FALSE,cols=rep("black",length(unique(myData_SIP$celltypes_new_2))))+ theme(legend.position = 'none')
patchwork::wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(myData_SIP, features = c("Ifng", "Il4", "Il13"),  group.by = "celltypes_new_2", pt.size = 1, cols=rep("black",length(unique(myData_SIP$celltypes_new_2))))+ theme(legend.position = 'none')
patchwork::wrap_plots(plots = plots, ncol = 1)

FeaturePlot(myData_SIP, features = "Ifng",split.by = "age")
FeaturePlot(myData_SIP, features = "Il4",split.by = "age")
FeaturePlot(myData_SIP, features = "Il13",split.by = "age")
# postive ratio test
myData_SIP_pos<-subset(myData_SIP, subset = Ifng > 1)
library(dplyr)
library(ggpubr)
myData_SIP_pos$Var1<-paste(myData_SIP_pos$celltypes_new_2,myData_SIP_pos$age,myData_SIP_pos$Index,sep="#")
md<-data.frame("Var1"=myData_SIP_pos$Var1)
data.ratio <- dplyr::group_by(md,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratio$celltype<-str_split_fixed(data.ratio$Var1,"#",3)[,1]
data.ratio$age<-str_split_fixed(data.ratio$Var1,"#",3)[,2]
data.ratio$index<-str_split_fixed(data.ratio$Var1,"#",3)[,3]
data.ratio2 <- data.ratio %>% dplyr::group_by(age,index) %>% dplyr::mutate(ratio2=count/sum(count))
# check 
sum(data.ratio2$count)
8/303
sum(data.ratio2[data.ratio2$age=="old"&data.ratio2$index=="1","count"])
8/106
# merge pos cell number and all cell number to calculate posivte ratio in each cell type 
my.data2.ratio2<-merge(my.data2,data.ratio2,by="Var1")
my.data2.ratio2$ratio_pos<-my.data2.ratio2$count/my.data2.ratio2$Freq

## draw figure
my.data2.ratio2$age.y<-factor(my.data2.ratio2$age.y,levels=c("young","old"))
ggbarplot(my.data2.ratio2, x = "celltype", y = "ratio_pos", add = "mean_se", color = "age.y", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = age.y), label = "p.signif")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## check other genes in positve ratio
# form function
pos_calculate <- function(myData_SIP_pos) {
myData_SIP_pos$Var1<-paste(myData_SIP_pos$celltypes_new_2,myData_SIP_pos$age,myData_SIP_pos$Index,sep="#")
md<-data.frame("Var1"=myData_SIP_pos$Var1)
data.ratio <- dplyr::group_by(md,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratio$celltype<-str_split_fixed(data.ratio$Var1,"#",3)[,1]
data.ratio$age<-str_split_fixed(data.ratio$Var1,"#",3)[,2]
data.ratio$index<-str_split_fixed(data.ratio$Var1,"#",3)[,3]
data.ratio2 <- data.ratio %>% dplyr::group_by(age,index) %>% dplyr::mutate(ratio2=count/sum(count))
my.data2.ratio2<-merge(my.data2,data.ratio2,by="Var1",all.x = T)
my.data2.ratio2$ratio_pos<-my.data2.ratio2$count/my.data2.ratio2$Freq
my.data2.ratio2[is.na(my.data2.ratio2$ratio_pos),"ratio_pos"]<-0
my.data2.ratio2$age.x<-factor(my.data2.ratio2$age.x,levels=c("young","old"))
ggbarplot(my.data2.ratio2, x = "celltypes_new_2", y = "ratio_pos", add = "mean_sd", color = "age.x", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = age.x), method="wilcox.test",paired=F,label = "p.signif")+scale_color_manual(values = c("grey","grey0"))+ggtitle("wilcox.test_pairedF")+xlab("")+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
}
myData_SIP_pos<-subset(myData_SIP, subset = Ifng > 1)
pos_calculate(myData_SIP_pos)
myData_SIP_pos<-subset(myData_SIP, subset = Il4 > 1)
pos_calculate(myData_SIP_pos)
myData_SIP_pos<-subset(myData_SIP, subset = Il13 > 1)
pos_calculate(myData_SIP_pos)
saveRDS(myData_SIP,"./data/myData_SIP.RDS")

## deg
myData_SIP$celltype_age<-paste(myData_SIP$celltypes_new_2,myData_SIP$age,sep=".")
celltypes<-unique(myData_SIP$celltypes_new_2)
for(i in celltypes){assign(paste("age.DEG.wilcox",i,sep = "_") , FindMarkers(myData_SIP, ident.1 = paste(i,"old",sep="."), ident.2 = paste(i,"young",sep="."), verbose = FALSE,group.by='celltype_age'))}
for(i in celltypes){write.csv(get(paste("age.DEG.wilcox",i,sep = "_")),paste(i,"age.DEG.wilcox.csv",sep = "_"))}

## immune markers
immune_markers<-read.csv("./data/immune_markers.csv",head=T,stringsAsFactors = F)
immune_markers$scm<-paste(immune_markers$source,immune_markers$cell_type,immune_markers$markers,sep="#")
top20_immune<-merge(top20_target,immune_markers,by.x='gene',by.y='markers',all =T)
top20_immune2<-top20_immune[(!is.na(top20_immune$cluster))&(!is.na(top20_immune$scm)),]
for(i in unique(top20_immune2$gene)){pdf(paste(i,".pdf",sep=""),width=3,height = 3.5);plots<-VlnPlot(myData_SIP, features = i,  group.by = "celltypes_new_2", pt.size = 0,cols=rep("grey",length(unique(myData_SIP$celltypes_new_2))))+NoLegend();print(plots);dev.off()}
write.csv(top20_immune2,"top20_immune2.csv")
write.csv(top20_immune,"top20_immune.csv")
write.csv(top20_target,"top20_target.csv")
##cytokine
cytokine<-read.csv("./data/cytokine_gene",header = T,stringsAsFactors = F)
cytokines<-cytokine[cytokine$gene!="Il17",]
cytokines<-cytokines[!cytokines %in% c("Il1","Ifnb1","Il25","Il9","Il8","Il3","Il12","Ifna2")]
for(i in unique(cytokines)){pdf(paste(i,".pdf",sep=""),width=3,height = 3.5);plots<-VlnPlot(myData_SIP, features = i,  group.by = "celltypes_new_2", pt.size = 0,cols=rep("grey",length(unique(myData_SIP$celltypes_new_2))))+NoLegend();print(plots);dev.off()}

## heatmap for immune
myData_SIP.immune<-subset(myData_SIP,features = top20_immune2$gene)
top20_immune3<-top20_immune2[!duplicated(top20_immune2$gene),]
myData_SIP.immune <- ScaleData(object = myData_SIP.immune, features = top20_immune3$gene)
DoHeatmap(myData_SIP.immune, assay = "RNA", features = top20_immune3$gene, size = 2, group.by='celltypes_new_2')+theme(axis.text.y = element_text(size = 4))+ scale_y_discrete(limit=top20_immune3$gene,labels=top20_immune3$scm)+ theme(
legend.title = element_text(size = 4),
legend.text = element_text(size = 4)
) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
# reorder
top20_immune3$num<-str_split_fixed(top20_immune3$cluster,"_",2)[,2]
top20_immune3<-top20_immune3[order(top20_immune3$num,decreasing = T),]


DoHeatmap(myData_SIP.immune, assay = "RNA", features = top20_immune3$gene, size = 2, group.by='celltypes_new_2')+theme(axis.text.y = element_text(size = 8))+ scale_y_discrete(limit=top20_immune3$gene,labels=top20_immune3$scm)+ theme(
legend.title = element_text(size = 4),
legend.text = element_text(size = 8)
) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))

# simplify heatmap to form heatmap2
top20_immune3_2<-top20_immune3[!(top20_immune3$scm %in% c("abcam#LMPP#Il7r","abcam#T_Cell#Cd3g","Ido_nonBT#Lymphocytes#Tcf7","Ido_nonBT#Lymphocytes#Trbc2","Alon_CD4#aTregs#Cd74","Ido_nonBT#cDC2s_Res#H2-DMb2","Ido_nonBT#pDCs#Ly6d")),]
ct<-c("B_cells","B_cells_early","Cytotoxic","T_cells_CD4","Naive","Naive_T_cell_NaiveIsg15","T_cells_CD8","NK_Cell","NK_cells","NK","NK_DC","Monocyte","Macrophages","Plasma_cells")
cts<-data.frame(ct=as.character(ct),orders=1:14)
for(i in 1:nrow(top20_immune3_2)){top20_immune3_2$num[i]<-cts[cts$ct==top20_immune3_2$cell_type[i],"orders"]}
top20_immune3_2<-top20_immune3_2[order(as.numeric(top20_immune3_2$num),decreasing = T),]
DoHeatmap(myData_SIP.immune, assay = "RNA", features = top20_immune3_2$gene, size = 2, group.by='celltypes_new_2')+theme(axis.text.y = element_text(size = 8))+ scale_y_discrete(limit=top20_immune3_2$gene,labels=top20_immune3_2$scm)+ theme(
legend.title = element_text(size = 4),
legend.text = element_text(size = 8)
) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))

## redo postive cell ratio
pos_calculate2 <- function(myData_SIP_pos,yp=NULL) {
    myData_SIP_pos$Var1<-paste(myData_SIP_pos$celltypes_new_2,myData_SIP_pos$age,myData_SIP_pos$Index,sep="#")
    md<-data.frame("Var1"=myData_SIP_pos$Var1)
    data.ratio <- dplyr::group_by(md,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
    data.ratio$celltype<-str_split_fixed(data.ratio$Var1,"#",3)[,1]
    data.ratio$age<-str_split_fixed(data.ratio$Var1,"#",3)[,2]
    data.ratio$index<-str_split_fixed(data.ratio$Var1,"#",3)[,3]
    data.ratio2 <- data.ratio %>% dplyr::group_by(age,index) %>% dplyr::mutate(ratio2=count/sum(count))
    my.data2.ratio2<-merge(my.data2,data.ratio2,by="Var1",all.x = T)
    my.data2.ratio2[is.na(my.data2.ratio2$ratio_pos),"ratio_pos"]<-0
    my.data2.ratio2$ratio_pos<-my.data2.ratio2$count/my.data2.ratio2$Freq
    my.data2.ratio2$age.x<-factor(my.data2.ratio2$age.x,levels=c("young","old"))
    ggbarplot(my.data2.ratio2, x = "celltypes_new_2", y = "ratio_pos", add = "mean_sd", color = "age.x", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = age.x), method="wilcox.test",paired=F,label = "p.format",label.y = yp)+scale_color_manual(values = c("grey","grey0"))+ggtitle("wilcox.test_pairedF")+xlab("")+theme(axis.text.x=element_text(angle=45,hjust=0.5,vjust=0.5))+ylab("Posive-cell ratio")
}
myData_SIP_pos<-subset(myData_SIP, subset = Ifng > 1)
pos_calculate2(myData_SIP_pos)
myData_SIP_pos<-subset(myData_SIP, subset = Il4 > 1)
pos_calculate2(myData_SIP_pos)
myData_SIP_pos<-subset(myData_SIP, subset = Il13 > 1)
pos_calculate2(myData_SIP_pos)

# version 2
load("./data/Point4.RData")
library(Seurat)
DimPlot(object = myData_SIP, reduction = 'tsne',split.by = "age")+coord_fixed(ratio = 1)
### odds ratio
tg2<-tg %>% dplyr::group_by(celltypes_new_2,age)  %>% dplyr::mutate(celltypes_sum = sum(Freq))
# check
sum(tg[tg$age=="old"&tg$celltypes_new_2=="B-cell","Freq"])
sum(tg[tg$age=="young"&tg$celltypes_new_2=="B-cell","Freq"])

tg2<-unique(tg2[,c(3,4,8)])
tg2<-tg2 %>% dplyr::group_by(age)  %>% dplyr::mutate(celltypes_sum_ratio = celltypes_sum/sum(celltypes_sum))
# check
sum(tg2[tg2$age=="old","celltypes_sum"])
2151/4574
63/4574

library(reshape2)
library(tidyr)
#tg2_wide<-spread(tg2,age,celltypes_sum_ratio)
tg2_wide<-dcast(tg2,celltypes_new_2~tg2$age,value.var = 'celltypes_sum_ratio')
tg2_wide$odds_ratio<-(tg2_wide$old/(1-tg2_wide$old))/(tg2_wide$young/(1-tg2_wide$young))
tg2<-tg2 %>% dplyr::group_by(age)  %>% dplyr::mutate(age_sum = sum(celltypes_sum))
tg2_wide2<-dcast(tg2,celltypes_new_2~tg2$age,value.var = c("celltypes_sum"))
tg2_wide2_2<-dcast(tg2,celltypes_new_2~tg2$age,value.var = c("age_sum"))
tg2_wide2<-merge(tg2_wide2,tg2_wide2_2,by="celltypes_new_2",suffixes=c(".specific",".total"))
tg2_wide2$k<-(tg2_wide2$young.specific+tg2_wide2$old.specific)
colnames(tg2_wide2)<-c("celltypes_new_2","young.specific","q_old.specific","n_young.total","m_old.total","k")
## check
phyper(2151,4574,4423,4910)
phyper(2151,4910,4574,4423)

tg2_wide2$hypergeometric_p<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k)
colnames(tg2_wide2)[2]<-"x_young.specific"

fisher_conf<-function(q,m,n,x){tmp<-fisher.test(matrix(c(q,(m-q),x,(n-x)),nrow=2))
return(tmp$conf.int[c(1,2)])}
fisher_conf(2151,4574,4423,2759)

for (i in 1:nrow(tg2_wide2)){tg2_wide2$conf_1[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[1];tg2_wide2$conf_2[i]<-fisher_conf(tg2_wide2[i,"q_old.specific"],tg2_wide2[i,"m_old.total"],tg2_wide2[i,"n_young.total"],tg2_wide2[i,"x_young.specific"])[2]}
tg2_wide2_m<-merge(tg2_wide,tg2_wide2,by = "celltypes_new_2",all = T)

library(ggplot2)
tg2_wide2$hypergeometric_p_upper<-phyper(tg2_wide2$q_old.specific,tg2_wide2$m_old.total,tg2_wide2$n_young.total,tg2_wide2$k,lower.tail = F)
tg2_wide2_m$hypergeometric_p_upper<-tg2_wide2$hypergeometric_p_upper
tg2_wide2_m$hypergeometric_p_tails<-pmin(tg2_wide2_m$hypergeometric_p,tg2_wide2_m$hypergeometric_p_upper)

ggplot(tg2_wide2_m, aes(x=celltypes_new_2, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + coord_flip()+scale_x_discrete(limits = rev(levels(factor(tg2_wide2_m$celltypes_new_2))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)

## cytokines version2
# pdf
for(i in unique(cytokines)){pdf(paste(i,".pdf",sep=""),width=3,height = 3.5);plots<-VlnPlot(myData_SIP, features = i, split.by = "age", group.by = "celltypes_new_2", pt.size = 0)+ggplot2::scale_fill_manual(values = c('gray','gray0')) + geom_violin(colour = "NA")+ theme(legend.position = 'top');print(plots);dev.off()}
# eps
for(i in unique(cytokines)){plots<-VlnPlot(myData_SIP, features = i, split.by = "age", group.by = "celltypes_new_2", pt.size = 0)+ggplot2::scale_fill_manual(values = c('gray','gray0')) + geom_violin(colour = "NA")+ theme(legend.position = 'top');ggsave(file=paste(i,".eps",sep=""), plot=plots, width=3, height=3.5)}
# tiff
for(i in unique(cytokines)){plots<-VlnPlot(myData_SIP, features = i, split.by = "age", group.by = "celltypes_new_2", pt.size = 0)+ggplot2::scale_fill_manual(values = c('gray','gray0')) + geom_violin(colour = "NA")+ theme(legend.position = 'top');ggsave(file=paste(i,".tiff",sep=""), plot=plots, width=3, height=3.5)}

# heatmap
myData_SIP.cytokines<-subset(myData_SIP,features = unique(cytokines))
myData_SIP.cytokines <- ScaleData(object = myData_SIP.cytokines,features = unique(cytokines))
DoHeatmap(myData_SIP.cytokines, assay = "RNA", features = unique(cytokines), size = 2, group.by='celltypes_new_2')+theme(axis.text.y = element_text(size = 8))+ scale_y_discrete(limit=unique(cytokines),labels=unique(cytokines))+ theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 8)
) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
myData_SIP.cytokines$celltype_age<-factor(myData_SIP.cytokines$celltype_age,levels=c("B-cell.young","B-cell.old","B-early.young","B-early.old","CD4-Cytotoxic.young","CD4-Cytotoxic.old","CD4-Gata3.young","CD4-Gata3.old","CD4-naive.young","CD4-naive.old","CD8_NK.young","CD8_NK.old","Myeloblast.young","Myeloblast.old","Plasma-cell.young","Plasma-cell.old","Th17.young","Th17.old"))
DoHeatmap(myData_SIP.cytokines, assay = "RNA", features = unique(cytokines), size = 2, group.by='celltype_age')+theme(axis.text.y = element_text(size = 8))+ scale_y_discrete(limit=unique(cytokines),labels=unique(cytokines))+ theme(
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 8)
) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))

## deg
age.DEG.wilcox_all<-FindMarkers(myData_SIP, ident.1 = "old", ident.2 = "young", verbose = FALSE,group.by='age')
write.csv(age.DEG.wilcox_all,"./data_output/age.DEG.wilcox_all.csv")
min(abs(age.DEG.wilcox_all$avg_logFC))
## label p_val_adj
deg_bar<-function(deg){
p<-ggplot(data=na.omit(deg[unique(cytokines),]), aes(x=row.names(na.omit(deg[unique(cytokines),])), y=avg_logFC)) + geom_bar(stat="identity") + coord_flip()+theme_classic(base_size = 20)+geom_text(aes(label=signif(p_val_adj,digits = 3), hjust=ifelse(avg_logFC < 0 , -0.1, 1)))+xlab("")
print(p)
}
for(i in celltypes){
pdf(paste(i,"age.DEG.cytokines.pdf",sep = "_"),width=4,height = 3)
deg_bar(get(paste("age.DEG.wilcox",i,sep = "_")))
dev.off()
}
## label p_val
deg_bar<-function(deg){
p<-ggplot(data=na.omit(deg[unique(cytokines),]), aes(x=row.names(na.omit(deg[unique(cytokines),])), y=avg_logFC)) + geom_bar(stat="identity") + coord_flip()+theme_classic(base_size = 20)+geom_text(aes(label=signif(p_val,digits = 3), hjust=ifelse(avg_logFC < 0 , -0.1, 1)))+xlab("")
print(p)
}
for(i in celltypes){
pdf(paste(i,"age.DEG.cytokines.pdf",sep = "_"),width=4,height = 3)
deg_bar(get(paste("age.DEG.wilcox",i,sep = "_")))
dev.off()
}
#na.omit(age.DEG.wilcox_Th17[row.names(age.DEG.wilcox_Th17)%in%unique(cytokines),])
## bulk immune cells
ggplot(data=na.omit(age.DEG.wilcox_all[row.names(age.DEG.wilcox_all)%in%unique(cytokines),]), aes(x=row.names(na.omit(age.DEG.wilcox_all[row.names(age.DEG.wilcox_all)%in%unique(cytokines),])), y=avg_logFC)) + geom_bar(stat="identity") + coord_flip()+theme_classic(base_size = 20)+geom_text(aes(label=signif(p_val,digits = 3), hjust=ifelse(avg_logFC < 0 , -0.1, 1)))+xlab("")
VlnPlot(myData_SIP, features = "Ccl5",  group.by = "age", pt.size = 0,cols=rep("grey",length(unique(myData_SIP$age))))+NoLegend()+ggplot2::scale_fill_manual(values = c('gray','gray0'))
VlnPlot(myData_SIP, features = "Il22",  group.by = "age", pt.size = 0,cols=rep("grey",length(unique(myData_SIP$age))))+NoLegend()+ggplot2::scale_fill_manual(values = c('gray','gray0'))
save.image(file = "./data/Point4.RData")

## check positive cell ratio for all cytokines
for (i in cytokines){tryCatch({print(i);expr <- FetchData(object = myData_SIP, vars = i);p<-pos_calculate2(myData_SIP[, which(x = expr > 1)]);ggsave(paste(i,"_PositiveCellRatio.pdf",sep=""),width = 7.05, height = 6.03, units = "in")},error=function(e){})}

## Average expression in each cluster by sample
myData_SIP.cytokines$celltypes_new_2_Index_age<-paste(myData_SIP.cytokines$celltypes_new_2,myData_SIP.cytokines$Index_age,sep="#")
myData_SIP.cytokines_Avg<-AverageExpression(myData_SIP.cytokines,add.ident = "celltypes_new_2_Index_age")
myData_SIP.cytokines_Avg2<-as.data.frame(t(myData_SIP.cytokines_Avg$RNA))
myData_SIP.cytokines_Avg2$celltypes_new_2<-str_split_fixed(row.names(myData_SIP.cytokines_Avg2),"#",2)[,1]
myData_SIP.cytokines_Avg2$Index<-str_split_fixed(row.names(myData_SIP.cytokines_Avg2),"#",2)[,2]
myData_SIP.cytokines_Avg2$age<-str_split_fixed(myData_SIP.cytokines_Avg2$Index,"_",2)[,2]
ggbarplot(myData_SIP.cytokines_Avg2, x = "celltypes_new_2", y = "Il17a", add = "mean_sd", color = "age", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = age), method="wilcox.test",paired=F,label = "p.format",label.y = NULL)+scale_color_manual(values = c("grey","grey0"))+ggtitle("wilcox.test_pairedF")+xlab("")+theme(axis.text.x=element_text(angle=45,hjust=0.5,vjust=0.5))+ylab("average expression")
myData_SIP.cytokines_Avg2$age<-factor(myData_SIP.cytokines_Avg2$age,levels=c("young","old"))
for (i in cytokines){tryCatch({print(i);p<-ggbarplot(myData_SIP.cytokines_Avg2, x = "celltypes_new_2", y = i, add = "mean_sd", color = "age", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = age), method="wilcox.test",paired=F,label = "p.format",label.y = NULL)+scale_color_manual(values = c("grey","grey0"))+ggtitle("wilcox.test_pairedF")+xlab("")+theme(axis.text.x=element_text(angle=45,hjust=0.5,vjust=0.5))+ylab("average expression");ggsave(paste(i,"_AverageExpression.pdf",sep=""),width = 7.05, height = 6.03, units = "in")},error=function(e){})}
write.csv(myData_SIP.cytokines_Avg2,"./data_output/myData_SIP.cytokines_Avg2.csv")
## Average expression in each cluster by age
myData_SIP.cytokines_AvgS<-AverageExpression(myData_SIP.cytokines,add.ident = "celltype_age")
myData_SIP.cytokines_AvgS2<-as.data.frame(t(myData_SIP.cytokines_AvgS$RNA))
myData_SIP.cytokines_AvgS2$celltypes_new_2<-str_split_fixed(row.names(myData_SIP.cytokines_AvgS2),"\\.",2)[,1]
myData_SIP.cytokines_AvgS2$age<-str_split_fixed(row.names(myData_SIP.cytokines_AvgS2),"\\.",2)[,2]
write.csv(myData_SIP.cytokines_AvgS2,"./data_output/myData_SIP.cytokines_AvgS2.csv")

############################################################
############################################################
############################################################
# Modification 1
############################################################
############################################################
############################################################
load("./data/Point4.RData")
library(dplyr)
library(stringr)
library(ggpubr)
library(Seurat)

volcano_immune<-function(data,clname,cf=0.5,mark="sig") {
data$geneID<-row.names(data)
require(ggrepel)
data$log2FoldChange<-(data$avg_logFC)
data$sig[(data$p_val_adj > 0.05|data$p_val_adj=="NA")|(data$log2FoldChange < cf)& data$log2FoldChange > (-cf)] <- "no"
data$sig[data$p_val_adj <= 0.05 & data$log2FoldChange >= cf] <- "up"
data$sig[data$p_val_adj <= 0.05 & data$log2FoldChange <= (-cf)] <- "down"
data$sign<-NA
if (mark=="sig") {
data[data$sig %in% c("up","down"),"sign"] <- as.character(data[data$sig %in% c("up","down"),"geneID"])
} else {
data[data$geneID %in% cytokines,"sign"] <- as.character(data[data$geneID %in% cytokines,"geneID"])
}
p<-ggplot(data[!is.na(data$sig),],aes(log2FoldChange,-1*log10(p_val_adj),color = sig))+geom_point(size=1)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(deparse(quote(data)))+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.5, "lines"), show.legend = F, size = 3)
print(p)
}
# bubble plot in all celltypes
pdf("All_age.DEG.wilcox.pdf")
volcano_immune(age.DEG.wilcox_all,cf=0.5,mark = "sig")
dev.off()
pdf("All_age.DEG.wilcox_2.pdf")
volcano_immune(age.DEG.wilcox_all,cf=0.5,mark = "sig2")
dev.off()
## bubble plot in specific celltypes
for(i in celltypes){pdf(paste(i,"age.DEG.wilcox.pdf",sep = "_"));volcano_immune(get(paste("age.DEG.wilcox",i,sep = "_")),cf=0.5,mark="sig");dev.off()}
for(i in celltypes){pdf(paste(i,"age.DEG.wilcox_2.pdf",sep = "_"));volcano_immune(get(paste("age.DEG.wilcox",i,sep = "_")),cf=0.5,mark="sig2");dev.off()}

## MeanExpression by sample
myData_SIP$celltypes_new_2_Index_age<-paste(myData_SIP$celltypes_new_2,myData_SIP$Index_age,sep="#")
myData_SIP.allGenes_Avg_bySample<-AverageExpression(myData_SIP,add.ident = "celltypes_new_2_Index_age")
myData_SIP.allGenes_Avg2_bySample<-as.data.frame(t(myData_SIP.allGenes_Avg_bySample$RNA))
write.csv(myData_SIP.allGenes_Avg2_bySample,"./data_output/myData_SIP.allGenes_Avg2_bySample.csv")
## MeanExpression by age
myData_SIP$celltypes_new_2_age<-paste(myData_SIP$celltypes_new_2,myData_SIP$age,sep="#")
myData_SIP.allGenes_Avg_byAge<-AverageExpression(myData_SIP,add.ident = "celltypes_new_2_age")
myData_SIP.allGenes_Avg2_byAge<-as.data.frame(t(myData_SIP.allGenes_Avg_byAge$RNA))
write.csv(myData_SIP.allGenes_Avg2_byAge,"./data_output/myData_SIP.allGenes_Avg2_byAge.csv")

###### heatmap
load("./data/Point4.RData")
library(Seurat)
library(dplyr)
myData_SIP.markers_new <- FindAllMarkers(object = myData_SIP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
myData_SIP.markers_new %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) -> top20_target_new
#myData_SIP.top20<-subset(myData_SIP,features = top20_target$gene)
myData_SIP.top20<-subset(myData_SIP,features = top20_target_new$gene)
myData_SIP.top20.counts<-myData_SIP.top20[["RNA"]]@counts
myData_SIP.celltypes_new_2<-as.data.frame(myData_SIP$celltypes_new_2)
colnames(myData_SIP.celltypes_new_2)<-"celltypes_new_2"
myData_SIP.celltypes_new_2$cell_id<-row.names(myData_SIP.celltypes_new_2)
colnames(myData_SIP.top20.counts)<-myData_SIP.celltypes_new_2[match(colnames(myData_SIP.top20.counts),myData_SIP.celltypes_new_2$cell_id),"celltypes_new_2"]
OrderCol<-data.frame(a=c("B-cell","CD4-naive","Th17","CD4-Cytotoxic","Plasma-cell","Myeloblast","B-early","CD4-Gata3","CD8_NK"),b=c("A#B-cell","B#CD4-naive","C#Th17","D#CD4-Cytotoxic","E#Plasma-cell","F#Myeloblast","G#B-early","H#CD4-Gata3","I#CD8_NK"))
colnames(myData_SIP.top20.counts)<-OrderCol[match(colnames(myData_SIP.top20.counts),OrderCol$a),"b"]
myData_SIP.top20.counts<-myData_SIP.top20.counts[,order(colnames(myData_SIP.top20.counts))]
my_sample_col <- data.frame(sample = colnames(myData_SIP.top20.counts))
row.names(my_sample_col)<-make.names(colnames(myData_SIP.top20.counts),unique = T)
myData_SIP.top20.counts.c<-myData_SIP.top20.counts
colnames(myData_SIP.top20.counts)<-make.names(colnames(myData_SIP.top20.counts),unique = T)

top20_target_new$cluster2<-OrderCol[match(top20_target_new$cluster,OrderCol$a),"b"]
top20_target_new<-top20_target_new[order(top20_target_new$cluster2),]
#my_gene_col <- data.frame(gene = top20_target$cluster)
#myData_SIP.top20.counts<-myData_SIP.top20.counts[match(top20_target$gene,row.names(myData_SIP.top20.counts)),]
my_gene_col <- data.frame(gene = top20_target_new$cluster)
myData_SIP.top20.counts<-myData_SIP.top20.counts[match(top20_target_new$gene,row.names(myData_SIP.top20.counts)),]
row.names(my_gene_col)<-make.names(row.names(myData_SIP.top20.counts),unique=T)
row.names(myData_SIP.top20.counts)<-make.names(row.names(myData_SIP.top20.counts),unique=T)
saveRDS(myData_SIP.top20.counts,"./data//myData_SIP.top20.counts.RDS")
library(pheatmap)
myData_SIP.top20.counts.log2<-log2(as.data.frame(myData_SIP.top20.counts)+0.01)
pheatmap(myData_SIP.top20.counts.log2, annotation_row = my_gene_col, annotation_col = my_sample_col,cluster_rows = F,cluster_cols = F,scale = "row",color = colorRampPalette(c( "white", "red"))(10),show_colnames = FALSE,breaks=seq(0,1,by=0.1))
write.csv(myData_SIP.markers_new,"./data/myData_SIP.markers_new.csv")
write.csv(top20_target_new,"./data/top20_target_new.csv")
save.image("./data/Point4.RData")

####################################
#### Modification 2          #######
####################################

## celltype order in change
# SingleImmune_OddsRatio_NewOrder.pdf

tg2_wide2_m$celltypes_new_2 <-factor(tg2_wide2_m$celltypes_new_2,levels = c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8_NK","Myeloblast","CD4-Gata3","Th17"))
ggplot(tg2_wide2_m, aes(x=celltypes_new_2, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + coord_flip()+scale_x_discrete(limits = rev(levels(factor(tg2_wide2_m$celltypes_new_2))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)

# Cytokines_AverageExpression.pdf
Idents(myData_SIP.cytokines)<-myData_SIP.cytokines$celltypes_new_2
myData_SIP.cytokines_Avg<-AverageExpression(myData_SIP.cytokines,add.ident = "Index_age")
myData_SIP.cytokines_Avg2<-as.data.frame(t(myData_SIP.cytokines_Avg$RNA))
myData_SIP.cytokines_Avg2$celltypes_new_2<-str_split_fixed(myData_SIP.cytokines_Avg2$celltypes_new_2,"_",3)[,1]
myData_SIP.cytokines_Avg2[myData_SIP.cytokines_Avg2$celltypes_new_2=="CD8","celltypes_new_2"]<-"CD8_NK"
myData_SIP.cytokines_Avg2$celltypes_new_2 <-factor(myData_SIP.cytokines_Avg2$celltypes_new_2,levels = c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8_NK","Myeloblast","CD4-Gata3","Th17"))
myData_SIP.cytokines_Avg2$age<-gsub("^.*_","",row.names(myData_SIP.cytokines_Avg2))
myData_SIP.cytokines_Avg2$age<-factor(myData_SIP.cytokines_Avg2$age,levels = c("young","old"))
for (i in cytokines){tryCatch({print(i);p<-ggbarplot(myData_SIP.cytokines_Avg2, x = "celltypes_new_2", y = i, add = "mean_sd", color = "age", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = age), method="wilcox.test",paired=F,label = "p.format",label.y = NULL)+scale_color_manual(values = c("grey","grey0"))+ggtitle("wilcox.test_pairedF")+xlab("")+theme(axis.text.x=element_text(angle=45,hjust=0.5,vjust=0.5))+ylab(paste(i,"average expression"));ggsave(paste(i,"_AverageExpression.pdf",sep=""),width = 7.05, height = 6.03, units = "in")},error=function(e){})}
write.csv(myData_SIP.cytokines_Avg2,"./data_output/myData_SIP.cytokines_Avg2_new.csv")
# AveExp_Cytotoxic.pdf AveExp_ILC2.pdf
myData_SIP.cytokines_Avg_byAge<-AverageExpression(myData_SIP.cytokines,add.ident = "age")
myData_SIP.cytokines_Avg2_byAge<-as.data.frame(t(myData_SIP.cytokines_Avg_byAge$RNA))
write.csv(myData_SIP.cytokines_Avg2_byAge,"./data_output/myData_SIP.cytokines_AvgS2_byAge.csv")

## improved heatmap
# SingleImmune_heatmap_top20_3.pdf
library(pheatmap)
my_gene_col_2 <- data.frame(gene = top20_target_new$cluster2)
row.names(my_gene_col_2)<-make.names(row.names(myData_SIP.top20.counts),unique=T)
my_sample_col_2<-my_sample_col
my_sample_col_2$sample<-as.character(my_sample_col_2$sample)
my_sample_col_2[my_sample_col_2$sample=="C#Th17","sample"]<-"C#ILC3"
my_sample_col_2[my_sample_col_2$sample=="H#CD4-Gata3","sample"]<-"H#ILC2"
my_gene_col_2$gene<-as.character(my_gene_col_2$gene)
my_gene_col_2[my_gene_col_2$gene=="C#Th17","gene"]<-"C#ILC3"
my_gene_col_2[my_gene_col_2$gene=="H#CD4-Gata3","gene"]<-"H#ILC2"
ann_colors = list(gene = c('A#B-cell'="firebrick2",'B#CD4-naive'="orange2",'C#ILC3'="khaki4",'D#CD4-Cytotoxic'="limegreen",'E#Plasma-cell'="mediumseagreen",'F#Myeloblast'="lightseagreen",'G#B-early'="steelblue2",'H#ILC2'="slateblue2",'I#CD8_NK'="plum4"), sample = c('A#B-cell'="firebrick2",'B#CD4-naive'="orange2",'C#ILC3'="khaki4",'D#CD4-Cytotoxic'="limegreen",'E#Plasma-cell'="mediumseagreen",'F#Myeloblast'="lightseagreen",'G#B-early'="steelblue2",'H#ILC2'="slateblue2",'I#CD8_NK'="plum4"))
pheatmap(myData_SIP.top20.counts.log2,annotation_colors=ann_colors, annotation_row = my_gene_col_2, annotation_col = my_sample_col_2,cluster_rows = F,cluster_cols = F,scale = "row",color = colorRampPalette(c( "white", "red"))(10),show_colnames = FALSE,breaks=seq(0,1,by=0.1))

## bubble plot with mean expression and log2fc
# bubble_sig_0.4.pdf bubble_cytokines_0.4.pdf
########

myData_SIP.allGenes_Avg_byCellType<-AverageExpression(myData_SIP)
age.DEG.wilcox.CD4.Cytotoxic.meanExp<-merge(tmp,myData_SIP.allGenes_Avg_byCellType,all.x = T,by=0)
myData_SIP.allGenes_Avg_byCellType<-as.data.frame(myData_SIP.allGenes_Avg_byCellType$RNA)
volcano_immune_2<-function(data,clname,cf=0.5,mark="sig",ct="RNA.CD4.Cytotoxic") {
tmp<-merge(get(data),myData_SIP.allGenes_Avg_byCellType,all.x = T,by=0)
tmp$geneID<-tmp$Row.names
require(ggrepel)
tmp$log2FoldChange<-(tmp$avg_logFC)
tmp$sig[(tmp$log2FoldChange < cf)& tmp$log2FoldChange > (-cf)] <- "no"
tmp$sig[tmp$log2FoldChange >= cf] <- "up"
tmp$sig[tmp$log2FoldChange <= (-cf)] <- "down"
tmp$sign<-NA
if (mark=="sig") {
tmp[tmp$sig %in% c("up","down"),"sign"] <- as.character(tmp[tmp$sig %in% c("up","down"),"geneID"])
} else {
tmp[tmp$geneID %in% cytokines,"sign"] <- as.character(tmp[tmp$geneID %in% cytokines,"geneID"])
}
p<-ggplot(tmp[!is.na(tmp$sig),],aes(get(ct),log2FoldChange,color = sig))+geom_point(size=1)+labs(x="mean_expression",y="log2(FoldChange)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+ggtitle(ct)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.5, "lines"), show.legend = F, size = 3)
print(p)
}

pdf("bubble_cytokines.pdf")
for(i in celltypes){
volcano_immune_2(paste("age.DEG.wilcox",i,sep = "_"),cf=0.4,ct=i,mark="cytokines")
}
dev.off()


pdf("bubble_sig.pdf")
for(i in celltypes){
volcano_immune_2(paste("age.DEG.wilcox",i,sep = "_"),cf=0.4,ct=i)
}
dev.off()
# SingleImmune_bubble_ILC2.pdf SingleImmune_bubble_CD4Cytotoxic.pdf
volcano_immune_2('age.DEG.wilcox_CD4-Gata3',cf=0.4,ct="CD4-Gata3",mark="cytokines")
volcano_immune_2('age.DEG.wilcox_CD4-Cytotoxic',cf=0.4,ct="CD4-Cytotoxic",mark="cytokines")


#### remind: Seurat use natural log (e) not log2

### DESeq2
for(i in celltypes){assign(paste("age.DEG.DESeq2",i,sep = "_") , FindMarkers(myData_SIP, ident.1 = paste(i,"old",sep="."), ident.2 = paste(i,"young",sep="."),test.use = "DESeq2",verbose = FALSE,group.by='celltype_age'))}
for(i in celltypes){write.csv(get(paste("age.DEG.DESeq2",i,sep = "_")),paste(i,"age.DEG.DESeq2.csv",sep = "_"))}

# volcano plot with logFC with mean expression
pdf("bubble_sig_DESeq2.pdf")
for(i in celltypes){
volcano_immune_2(paste("age.DEG.DESeq2",i,sep = "_"),cf=0.4,ct=i)
}
dev.off()
pdf("bubble_cytokines_DESeq2.pdf")
for(i in celltypes){
volcano_immune_2(paste("age.DEG.DESeq2",i,sep = "_"),cf=0.4,ct=i,mark="cytokines")
}
dev.off()

## volcano plot with logFC with p value and logFC
volcano_immune<-function(data,clname,cf=0.5,mark="sig") {
data$geneID<-row.names(data)
require(ggrepel)
data$log2FoldChange<-(data$avg_logFC)
data$sig[(data$p_val_adj > 0.05|data$p_val_adj=="NA")|(data$log2FoldChange < cf)& data$log2FoldChange > (-cf)] <- "no"
data$sig[data$p_val_adj <= 0.05 & data$log2FoldChange >= cf] <- "up"
data$sig[data$p_val_adj <= 0.05 & data$log2FoldChange <= (-cf)] <- "down"
data$sign<-NA
if (mark=="sig") {
data[data$sig %in% c("up","down"),"sign"] <- as.character(data[data$sig %in% c("up","down"),"geneID"])
} else {
data[data$geneID %in% cytokines,"sign"] <- as.character(data[data$geneID %in% cytokines,"geneID"])
}
p<-ggplot(data[!is.na(data$sig),],aes(log2FoldChange,-1*log10(p_val_adj),color = sig))+geom_point(size=1)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(deparse(quote(data)))+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.5, "lines"), show.legend = F, size = 3)
print(p)
}

for(i in celltypes){pdf(paste(i,"age.DEG.DESeq2.pdf",sep = "_"));volcano_immune(get(paste("age.DEG.DESeq2",i,sep = "_")),cf=0.5,mark="sig");dev.off()}
for(i in celltypes){pdf(paste(i,"age.DEG.DESeq2_2.pdf",sep = "_"));volcano_immune(get(paste("age.DEG.DESeq2",i,sep = "_")),cf=0.5,mark="sig2");dev.off()}

save.image("./data/Point4.RData")

##############################
#### Cytokine content level
##############################
load("./data/Point4.RData")
library(Seurat)
# https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
myData_SIP.RawCounts<-as.data.frame(myData_SIP[["RNA"]]@counts) #  the count matrix is stored in pbmc[["RNA"]]@counts
# myData_SIP.RawCounts<-as.data.frame(myData_SIP[["RNA"]]@data) #  Normalized values are stored in pbmc[["RNA"]]@data
unique(colnames(myData_SIP.RawCounts)==myData_SIP$sampleID)
colnames(myData_SIP.RawCounts)<-paste(myData_SIP$Index_age,myData_SIP$celltypes_new_2,sep="#")
rowSums(myData_SIP.RawCounts["Ifng",colnames(myData_SIP.RawCounts)=="1_old#CD4-Cytotoxic"])
library(stringr)
# sum by sample_celltype
myData_SIP.RawCounts2<-t(rowsum(t(myData_SIP.RawCounts), colnames(myData_SIP.RawCounts)))
myData_SIP.RawCounts<-myData_SIP.RawCounts2
rm(myData_SIP.RawCounts2)
sum(myData_SIP.RawCounts["Ifng",str_split_fixed(colnames(myData_SIP.RawCounts),"#",2)[,1]=="1_old"])
# sum by sample
myData_SIP.RawCounts2<-t(rowsum(t(myData_SIP.RawCounts), str_split_fixed(colnames(myData_SIP.RawCounts),"#",2)[,1]))

> colSums(myData_SIP.RawCounts2)
    1_old     2_old   3_young   4_young     5_old   5_young     6_old   6_young 
  3486557.6 1956733.6 2389204.5 3174582.2 1089180.8  955679.7 1051032.6  522844.7 

> median(colSums(myData_SIP.RawCounts2))
1522957
# normalize raw counts by total counts of sample, and scale it with median total counts os 8 samples 
myData_SIP.RawCounts3<-sapply(colnames(myData_SIP.RawCounts),function(x){(myData_SIP.RawCounts[,x]/colSums(myData_SIP.RawCounts2)[str_split_fixed(x,"#",2)[,1]])*median(colSums(myData_SIP.RawCounts2))})
sum(myData_SIP.RawCounts["Ifng",colnames(myData_SIP.RawCounts)=="5_old#CD4-Cytotoxic"])
# group by age
myData_SIP.RawCounts4<-t(rowsum(t(myData_SIP.RawCounts3), str_split_fixed(colnames(myData_SIP.RawCounts3),"#",2)[,1]))
sum(myData_SIP.RawCounts3["Ifng",str_split_fixed(colnames(myData_SIP.RawCounts3),"#",2)[,1]=="6_old"])
myData_SIP.RawCounts4<-as.data.frame(myData_SIP.RawCounts4)
library(ggpubr)
# test
#test<-as.data.frame(t(myData_SIP.RawCounts4["Ifng",]))
#test$age<-str_split_fixed(row.names(test),"_",2)[,2]
#test$age<-factor(test$age,levels = c("young","old"))
#ggbarplot(test, x = "age", y = "Ifng", add = "mean_sd",color = "age", palette = c("grey","grey0"))+ stat_compare_means(method = "t.test")



## cytokine content between age
# SingleImmune_cytokines_content_aging_@data.pdf 
cytokine_bar<-function(myData_SIP.RawCounts4,cytokine){
test<-as.data.frame(t(myData_SIP.RawCounts4[cytokine,]))
test$age<-str_split_fixed(row.names(test),"_",2)[,2]
test$age<-factor(test$age,levels = c("young","old"))
p<-ggbarplot(test, x = "age", y = cytokine, add = "mean_sd",color = "age", palette = c("grey","grey0"))+ stat_compare_means(method = "t.test")+ylab(paste(cytokine," normalized counts (by sample),sd"))
print(p)
}
pdf("SingleImmune_cytokines_content_aging.pdf",width = 4)
for (i in unique(cytokines)) {cytokine_bar(myData_SIP.RawCounts4,i)}
dev.off()


#test<-as.data.frame(t(myData_SIP.RawCounts3["Ifng",]))
#test$age<-str_split_fixed(row.names(test),"_",2)[,2]
#test$age<-str_split_fixed(test$age,"#",2)[,1]
#test$celltype<-str_split_fixed(row.names(test),"#",2)[,2]
#test$celltype<-factor(test$celltype,levels = c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8_NK","Myeloblast","CD4-Gata3","Th17"))
#test$age<-factor(test$age,levels = c("young","old"))
#ggbarplot(test, x = "celltype", y = "Ifng", add = "mean_sd",color = "age", palette = c("grey","grey0"),position = position_dodge(0.8))+ stat_compare_means(aes(group = age),method = "t.test",size=2)+ylab(paste("IFng"," normalized counts (by sample),sd"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## cytokine content between age_celltype
# SingleImmune_cytokines_content_aging_celltypes_@data.pdf
cytokine_bar_2<-function(myData_SIP.RawCounts3,cytokine){
test<-as.data.frame(t(myData_SIP.RawCounts3[cytokine,]))
test$age<-str_split_fixed(row.names(test),"_",2)[,2]
test$age<-str_split_fixed(test$age,"#",2)[,1]
test$celltype<-str_split_fixed(row.names(test),"#",2)[,2]
test$celltype<-factor(test$celltype,levels = c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8_NK","Myeloblast","CD4-Gata3","Th17"))
test$age<-factor(test$age,levels = c("young","old"))
p<-ggbarplot(test, x = "celltype", y = cytokine, add = "mean_sd",color = "age", palette = c("grey","grey0"),position = position_dodge(0.8))+ stat_compare_means(aes(group = age),method = "t.test",size=2)+ylab(paste(cytokine," normalized counts (by sample),sd"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
}
myData_SIP.RawCounts3<-as.data.frame(myData_SIP.RawCounts3)
pdf("SingleImmune_cytokines_content_aging_celltypes.pdf")
for (i in unique(cytokines)) {cytokine_bar_2(myData_SIP.RawCounts3,i)}
dev.off()

save.image("./data/Point4.RData")

## add dot plot and tsne plot: 
# SingleImmune_Cytokines.pdf
load("./data/Point4.RData")
library(stringr)
library(ggpubr)
library(Seurat)
library(reshape2)
myData_SIP$celltypes_new_2<-factor(myData_SIP$celltypes_new_2,levels = rev(levels(myData_SIP$celltypes_new_2)))
DotPlot(myData_SIP, features = unique(cytokines),group.by = "celltypes_new_2") + RotatedAxis()

FeaturePlot(myData_SIP, features = "Ifng",split.by = "age",coord.fixed = 1)
for(i in unique(cytokines)){pdf(paste(i,"pdf",sep="."),width=9,height=6);print(FeaturePlot(myData_SIP, features = i,split.by = "age",coord.fixed = 1));dev.off()}

test<-as.data.frame(t(myData_SIP.RawCounts4))
test<-test[,unique(cytokines)]
test$age<-str_split_fixed(row.names(test),"_",2)[,2]
test$age<-factor(test$age,levels = c("young","old"))
test_long1<-melt(test,
id.vars = c('age'),
variable.name='cytokines',
value.name='norm_exp')

ggbarplot(test_long1, x = "cytokines", y = "norm_exp", add = "mean_sd",color = "age", palette = c("grey","grey0"),position = position_dodge(0.8))+ stat_compare_means(aes(group = age),method = "t.test",label = "p.signif",size=4)+ylab(paste("cytokines"," normalized counts (by sample),sd"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggbarplot(test_long1, x = "cytokines", y = "norm_exp", add = "mean_sd",color = "age", palette = c("grey","grey0"),position = position_dodge(0.8))+ stat_compare_means(aes(group = age),method = "t.test",size=1)+ylab(paste("cytokines"," normalized counts (by sample),sd"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
cytokine_bar_2(myData_SIP.RawCounts3,"Ccr5")

====================================================================
====================================================================
#### Modify1 
load("./data/Point4.RData")
library(Seurat)
## PCA
## omid's order
Idents(myData_SIP)<-factor(Idents(myData_SIP),levels = c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8_NK","Myeloblast","CD4-Gata3","Th17"))
mycolors2<-mycolors
mycolors2[7]<-"purple"
DimPlot(object = myData_SIP, reduction = 'pca')+scale_color_manual(values = mycolors2)+coord_fixed(ratio = 1)
DimPlot(object = myData_SIP, reduction = 'tsne')+scale_color_manual(values = mycolors2)+coord_fixed(ratio=1)
DimPlot(object = myData_SIP, reduction = 'tsne',split.by = "age")+scale_color_manual(values = mycolors2)+coord_fixed(ratio=1)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Ccl5"))+scale_fill_manual(values = mycolors2)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Klrg1"))+scale_fill_manual(values = mycolors2)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Cd117"))+scale_fill_manual(values = mycolors2)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Kit"))+scale_fill_manual(values = mycolors2)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Rorc"))+scale_fill_manual(values = mycolors2)
## change cell type order: new order
Idents(myData_SIP)<-factor(Idents(myData_SIP),levels = c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8_NK","Th17","CD4-Gata3","Myeloblast"))
mycolors3<-mycolors2
mycolors3<-mycolors3[-7]
mycolors3[9]<-"purple"
## check Ccl5 and Klrg1 markers: SingleImmune_Ccl5_cytotoxic.pdf SingleImmune_Klrg1_ILC2.pdf
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Ccl5"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Klrg1"))+scale_fill_manual(values = mycolors3)
## tsne with new color: SingleImmune_tsne_newColor_age.pdf SingleImmune_tsne_newColor.pdf 
## SingleImmune_tsne_splitAge_v2.pdf
DimPlot(object = myData_SIP, reduction = 'tsne',split.by = "age")+scale_color_manual(values = mycolors3)+coord_fixed(ratio=1)
## SingleImmune_tsne_res_0.4.pdf
DimPlot(object = myData_SIP, reduction = 'tsne')+scale_color_manual(values = mycolors3)+coord_fixed(ratio=1)
## PCA with new color: SingleImmune_PCA.pdf
DimPlot(object = myData_SIP, reduction = 'pca')+scale_color_manual(values = mycolors3)+coord_fixed(ratio=1)


## check ILC3 or Th17: all files in /Users/jlu/Desktop/Pro_Mine/Paper_Omid/Point4_1/figure/ILC3
## SingleImmune_Thy1-CD90_ILC3.pdf
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Rorc"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Kit"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Il22"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Il8"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Il17"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Csf2"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Stat3"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Ifng"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Ptprc"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Il7r"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Thy1"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Cd4"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Cxcr3"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Cd6"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Klrb1"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Ccr6"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Ifng"))+scale_fill_manual(values = mycolors3)
VlnPlot(myData_SIP,pt.size = 0.2, features = c("Ifngr1"))+scale_fill_manual(values = mycolors3)
expr <- FetchData(object = myData_SIP, vars = "Cd4")
## check ILC3 markers: CD254(RANKL), IL23r, IL1rb, CD122, CD314(NKG2D), CD49d, Ahr.
pos_calculate2(myData_SIP[, which(x = expr > 1)])
FeaturePlot(myData_SIP, features = c("Cd4"))
FeaturePlot(myData_SIP, features = c("Tnfsf11"))
FeaturePlot(myData_SIP, features = c("Il23r"))
FeaturePlot(myData_SIP, features = c("Il1r2"))
FeaturePlot(myData_SIP, features = c("Il2rb"))
FeaturePlot(myData_SIP, features = c("Klrk1"))
FeaturePlot(myData_SIP, features = c("Itga4"))
FeaturePlot(myData_SIP, features = c("Ahr"))



==============================================================
# Attached
# Heatmap for markers in each cluster
==============================================================

load("./data/Point4.RData")

#######################################################
##### Heatmap for top5 markers in each celltype #######
#######################################################
# SingleImmune_heatmap_top5.pdf
library(pheatmap)
pheatmap(myData_SIP.top20.counts.log2,annotation_colors=ann_colors, annotation_row = my_gene_col_2, annotation_col = my_sample_col_2,cluster_rows = F,cluster_cols = F,scale = "row",color = colorRampPalette(c( "white", "red"))(10),show_colnames = FALSE,breaks=seq(0,1,by=0.1))
library(Seurat)
library(dplyr)
myData_SIP.markers_new %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) -> top5_target_new
myData_SIP.top5<-subset(myData_SIP,features = top5_target_new$gene)
myData_SIP.top5.counts<-myData_SIP.top5[["RNA"]]@counts
#myData_SIP.celltypes_new_2<-as.data.frame(myData_SIP$celltypes_new_2)
#colnames(myData_SIP.celltypes_new_2)<-"celltypes_new_2"
#myData_SIP.celltypes_new_2$cell_id<-row.names(myData_SIP.celltypes_new_2)
colnames(myData_SIP.top5.counts)<-myData_SIP.celltypes_new_2[match(colnames(myData_SIP.top5.counts),myData_SIP.celltypes_new_2$cell_id),"celltypes_new_2"]
OrderCol<-data.frame(a=c("B-cell","CD4-naive","Th17","CD4-Cytotoxic","Plasma-cell","Myeloblast","B-early","CD4-Gata3","CD8_NK"),b=c("A#B-cell","B#CD4-naive","C#Th17","D#CD4-Cytotoxic","E#Plasma-cell","F#Myeloblast","G#B-early","H#CD4-Gata3","I#CD8_NK"))
colnames(myData_SIP.top5.counts)<-OrderCol[match(colnames(myData_SIP.top5.counts),OrderCol$a),"b"]
myData_SIP.top5.counts<-myData_SIP.top5.counts[,order(colnames(myData_SIP.top5.counts))]
my_sample_col <- data.frame(sample = colnames(myData_SIP.top5.counts))
row.names(my_sample_col)<-make.names(colnames(myData_SIP.top5.counts),unique = T)
#myData_SIP.top5.counts.c<-myData_SIP.top5.counts
colnames(myData_SIP.top5.counts)<-make.names(colnames(myData_SIP.top5.counts),unique = T)
top5_target_new$cluster2<-OrderCol[match(top5_target_new$cluster,OrderCol$a),"b"]
top5_target_new<-top5_target_new[order(top5_target_new$cluster2),]
#my_gene_col <- data.frame(gene = top20_target$cluster)
#myData_SIP.top5.counts<-myData_SIP.top5.counts[match(top20_target$gene,row.names(myData_SIP.top5.counts)),]
my_gene_col <- data.frame(gene = top5_target_new$cluster)
myData_SIP.top5.counts<-myData_SIP.top5.counts[match(top5_target_new$gene,row.names(myData_SIP.top5.counts)),]
row.names(my_gene_col)<-make.names(row.names(myData_SIP.top5.counts),unique=T)
row.names(myData_SIP.top5.counts)<-make.names(row.names(myData_SIP.top5.counts),unique=T)
library(pheatmap)
myData_SIP.top5.counts.log2<-log2(as.data.frame(myData_SIP.top5.counts)+0.01)
pheatmap(myData_SIP.top5.counts.log2, annotation_row = my_gene_col, annotation_col = my_sample_col,cluster_rows = F,cluster_cols = F,scale = "row",color = colorRampPalette(c( "white", "red"))(10),show_colnames = FALSE,breaks=seq(0,1,by=0.1))
## improved
my_gene_col_2 <- data.frame(gene = top5_target_new$cluster2)
row.names(my_gene_col_2)<-make.names(row.names(myData_SIP.top5.counts),unique=T)
my_sample_col_2<-my_sample_col
my_sample_col_2$sample<-as.character(my_sample_col_2$sample)
my_sample_col_2[my_sample_col_2$sample=="C#Th17","sample"]<-"C#ILC3"
my_sample_col_2[my_sample_col_2$sample=="H#CD4-Gata3","sample"]<-"H#ILC2"
my_gene_col_2$gene<-as.character(my_gene_col_2$gene)
my_gene_col_2[my_gene_col_2$gene=="C#Th17","gene"]<-"C#ILC3"
my_gene_col_2[my_gene_col_2$gene=="H#CD4-Gata3","gene"]<-"H#ILC2"
ann_colors = list(gene = c('A#B-cell'="firebrick2",'B#CD4-naive'="orange2",'C#ILC3'="khaki4",'D#CD4-Cytotoxic'="limegreen",'E#Plasma-cell'="mediumseagreen",'F#Myeloblast'="lightseagreen",'G#B-early'="steelblue2",'H#ILC2'="slateblue2",'I#CD8_NK'="plum4"), sample = c('A#B-cell'="firebrick2",'B#CD4-naive'="orange2",'C#ILC3'="khaki4",'D#CD4-Cytotoxic'="limegreen",'E#Plasma-cell'="mediumseagreen",'F#Myeloblast'="lightseagreen",'G#B-early'="steelblue2",'H#ILC2'="slateblue2",'I#CD8_NK'="plum4"))
pheatmap(myData_SIP.top5.counts.log2,annotation_colors=ann_colors, annotation_row = my_gene_col_2, annotation_col = my_sample_col_2,cluster_rows = F,cluster_cols = F,scale = "row",color = colorRampPalette(c( "white", "red"))(10),show_colnames = FALSE,breaks=seq(0,1,by=0.1))


#######################################################
##### Heatmap for cytokines in each celltype ##########
#######################################################
# values were z-score of normalized counts by sample in each cytokine


# SingleImmune_heatmap_cytokines_cytotoxic.pdf
heatmap_Cytotoxic<-myData_SIP.RawCounts3[,grepl("Cytotoxic",colnames(myData_SIP.RawCounts3))]
pheatmap(na.omit(heatmap_Cytotoxic[unique(cytokines),c("3_young#CD4-Cytotoxic","4_young#CD4-Cytotoxic","5_young#CD4-Cytotoxic","6_young#CD4-Cytotoxic","1_old#CD4-Cytotoxic","2_old#CD4-Cytotoxic","5_old#CD4-Cytotoxic","6_old#CD4-Cytotoxic")]),cluster_cols = F,scale = "row",cluster_rows = F)
hinput<-heatmap_Cytotoxic[unique(cytokines),c("3_young#CD4-Cytotoxic","4_young#CD4-Cytotoxic","5_young#CD4-Cytotoxic","6_young#CD4-Cytotoxic","1_old#CD4-Cytotoxic","2_old#CD4-Cytotoxic","5_old#CD4-Cytotoxic","6_old#CD4-Cytotoxic")]
pheatmap(hinput[rowSums(hinput)!=0,],cluster_cols = F,scale = "row",cluster_rows = T)
# SingleImmune_heatmap_cytokines_ILC2.pdf
celltype<-"CD4-Gata3"
heatmap_celltype<-myData_SIP.RawCounts3[,grepl(celltype,colnames(myData_SIP.RawCounts3))]
hinput<-heatmap_celltype[unique(cytokines),c("3_young#CD4-Gata3","4_young#CD4-Gata3","5_young#CD4-Gata3","6_young#CD4-Gata3","1_old#CD4-Gata3","2_old#CD4-Gata3","5_old#CD4-Gata3","6_old#CD4-Gata3")]
table(myData_SIP$celltypes_new_2,myData_SIP$Index_age)
hinput<-heatmap_celltype[unique(cytokines),c("3_young#CD4-Gata3","4_young#CD4-Gata3","5_young#CD4-Gata3","1_old#CD4-Gata3","2_old#CD4-Gata3","5_old#CD4-Gata3","6_old#CD4-Gata3")]
pheatmap(hinput[rowSums(hinput)!=0,],cluster_cols = F,scale = "row",cluster_rows = T)

# SingleImmune_heatmap_cytokines.pdf
heatmap_cytokines<-function(celltype){
heatmap_celltype<-myData_SIP.RawCounts3[,grepl(celltype,colnames(myData_SIP.RawCounts3))]
cn<-c(paste("3_young",celltype,sep="#"),paste("4_young",celltype,sep="#"),paste("5_young",celltype,sep="#"),paste("6_young",celltype,sep="#"),paste("1_old",celltype,sep="#"),paste("2_old",celltype,sep="#"),paste("5_old",celltype,sep="#"),paste("6_old",celltype,sep="#"))
cn<-intersect(cn,colnames(heatmap_celltype))
hinput<-heatmap_celltype[unique(cytokines),cn]
print(pheatmap(hinput[rowSums(hinput)!=0,],cluster_cols = F,scale = "row",cluster_rows = T))
}

pdf("heatmap_cytokines.pdf",width = 3.63,height = 4.81)
for(i in unique(myData_SIP$celltypes_new_2)){heatmap_cytokines(i)}
dev.off()


##########################################
################ Modify ##################
##########################################
## check cytokines in DEG
deg_wilcox<-read.csv("/Users/jlu/Desktop/Pro_Mine/Paper_Omid/Point4_1/data_output/DEG_wilcox/all.celltypes.DEG.wilcox.csv",head=T,stringsAsFactors = F)
deg_wilcox$CellTypes<-sapply(strsplit(deg_wilcox$CellTypes,"_"),'[',1)
load("./data/Point4.RData")
cytokines_pure<-unique(cytokines)[c(3,4,5,6,7,8,10:15,17:23)]
deg_wilcox_cytokines<-deg_wilcox[deg_wilcox$gene %in% cytokines_pure,]
unique(deg_wilcox_cytokines$gene)
###


myData_SIP.cytokines_Avg2_byAge_cytokinesNeri<-myData_SIP.cytokines_Avg2_byAge[,colnames(myData_SIP.cytokines_Avg2_byAge) %in% cytokines_pure]
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri<-t(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri<-myData_SIP.cytokines_Avg2_byAge_cytokinesNeri+0.01
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri<-myData_SIP.cytokines_Avg2_byAge_cytokinesNeri[,c(1:2,4,3,5:ncol(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri))]

myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc<- as.data.frame(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`B-cell`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,1]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,2])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`CD4-naive`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,3]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,4])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`Th17`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,5]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,6])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`CD4-Cytotoxic`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,7]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,8])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`Plasma-cell`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,9]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,10])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`Myeloblast`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,11]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,12])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`B-early`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,13]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,14])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`CD4-Gata3`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,15]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,16])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc$`CD8_NK`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,17]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,18])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc<-myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc[,-c(1:18)]
library(pheatmap)
pheatmap(t(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc),angle_col = 90,scale = "none")
## Heat map_avgExp_scaleCol.pdf
pheatmap(t(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_logfc),angle_col = 90,scale = "column")
#### psedo normalized counts
myData_SIP.RawCounts3_dots<-myData_SIP.RawCounts3
colnames(myData_SIP.RawCounts3_dots)<-sapply(strsplit(colnames(myData_SIP.RawCounts3_dots),"_"),'[',2)

myData_SIP.RawCounts3_dots<-as.data.frame(sapply(unique(colnames(myData_SIP.RawCounts3_dots)), function(col) rowMeans(myData_SIP.RawCounts3_dots[names(myData_SIP.RawCounts3_dots)==col])))
myData_SIP.RawCounts3_dots<-as.data.frame(myData_SIP.RawCounts3_dots)

myData_SIP.RawCounts3_dots$`logfc#B-cell`<-log2((myData_SIP.RawCounts3_dots$`old#B-cell`+0.01)/(myData_SIP.RawCounts3_dots$`young#B-early`+0.01))
myData_SIP.RawCounts3_dots$`logfc#B-cell`<-log2((myData_SIP.RawCounts3_dots$`old#B-cell`+0.01)/(myData_SIP.RawCounts3_dots$`young#B-cell`+0.01))
myData_SIP.RawCounts3_dots$`logfc#B-early`<-log2((myData_SIP.RawCounts3_dots$`old#B-early`+0.01)/(myData_SIP.RawCounts3_dots$`young#B-early`+0.01))
myData_SIP.RawCounts3_dots$`logfc#CD4-Cytotoxic`<-log2((myData_SIP.RawCounts3_dots$`old#CD4-Cytotoxic`+0.01)/(myData_SIP.RawCounts3_dots$`young#CD4-Cytotoxic`+0.01))
myData_SIP.RawCounts3_dots$`logfc#CD4-Gata3`<-log2((myData_SIP.RawCounts3_dots$`old#CD4-Gata3`+0.01)/(myData_SIP.RawCounts3_dots$`young#CD4-Gata3`+0.01))
myData_SIP.RawCounts3_dots$`logfc#CD4-naive`<-log2((myData_SIP.RawCounts3_dots$`old#CD4-naive`+0.01)/(myData_SIP.RawCounts3_dots$`young#CD4-naive`+0.01))
myData_SIP.RawCounts3_dots$`logfc#CD8`<-log2((myData_SIP.RawCounts3_dots$`old#CD8`+0.01)/(myData_SIP.RawCounts3_dots$`young#CD8`+0.01))
myData_SIP.RawCounts3_dots$`logfc#Myeloblast`<-log2((myData_SIP.RawCounts3_dots$`old#Myeloblast`+0.01)/(myData_SIP.RawCounts3_dots$`young#Myeloblast`+0.01))
myData_SIP.RawCounts3_dots$`logfc#Plasma-cell`<-log2((myData_SIP.RawCounts3_dots$`old#Plasma-cell`+0.01)/(myData_SIP.RawCounts3_dots$`young#Plasma-cell`+0.01))
myData_SIP.RawCounts3_dots$`logfc#Th17`<-log2((myData_SIP.RawCounts3_dots$`old#Th17`+0.01)/(myData_SIP.RawCounts3_dots$`young#Th17`+0.01))

myData_SIP.RawCounts3_dots_logfc<-myData_SIP.RawCounts3_dots[,grepl("logfc",colnames(myData_SIP.RawCounts3_dots))]
myData_SIP.RawCounts3_dots_expYoung<-myData_SIP.RawCounts3_dots[,grepl("young",colnames(myData_SIP.RawCounts3_dots))]
myData_SIP.RawCounts3_dots_logfc<-reshape2::melt(as.matrix(myData_SIP.RawCounts3_dots_logfc))
myData_SIP.RawCounts3_dots_expYoung<-reshape2::melt(as.matrix(myData_SIP.RawCounts3_dots_expYoung))

library(stringr)
myData_SIP.RawCounts3_dots_logfc_both$CellTypes<-str_split_fixed(myData_SIP.RawCounts3_dots_logfc_both$Var2.x,"#",2)[,2]
unique(myData_SIP.RawCounts3_dots_expYoung$Var1==myData_SIP.RawCounts3_dots_logfc$Var1)
unique(str_split_fixed(myData_SIP.RawCounts3_dots_expYoung$Var2,"#",2)[,2]==str_split_fixed(myData_SIP.RawCounts3_dots_logfc$Var2,"#",2)[,2])
myData_SIP.RawCounts3_dots_logfc_both<-cbind(myData_SIP.RawCounts3_dots_logfc,myData_SIP.RawCounts3_dots_expYoung)
myData_SIP.RawCounts3_dots_logfc_both$CellTypes<-str_split_fixed(myData_SIP.RawCounts3_dots_logfc_both$Var2,"#",2)[,2]

colnames(myData_SIP.RawCounts3_dots_logfc_both)[3]<-"logfc"
colnames(myData_SIP.RawCounts3_dots_logfc_both)[6]<-"expYoung"
colnames(myData_SIP.RawCounts3_dots_logfc_both)[4]<-"Var1.exp"
colnames(myData_SIP.RawCounts3_dots_logfc_both)[5]<-"Var2.exp"
colnames(myData_SIP.RawCounts3_dots_logfc_both)[2]<-"Var2.logfc"
colnames(myData_SIP.RawCounts3_dots_logfc_both)[1]<-"Var1.logfc"

myData_SIP.RawCounts3_dots_logfc_both[myData_SIP.RawCounts3_dots_logfc_both$CellTypes=="CD8","CellTypes"]<-"CD8-NK"
myData_SIP.RawCounts3_dots_logfc_both$CellTypes<-factor(myData_SIP.RawCounts3_dots_logfc_both$CellTypes,levels = rev(c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8-NK","Myeloblast","CD4-Gata3","Th17")))
## Points_psedoCounts_SizeColor_Cytokines_Exp_AgingChange.pdf
ggplot(myData_SIP.RawCounts3_dots_logfc_both[myData_SIP.RawCounts3_dots_logfc_both$Var1.logfc %in% cytokines_pure,], aes(Var1.logfc, CellTypes)) +
geom_point(aes(size = log2(expYoung+1), colour=logfc)) +scale_color_gradient2(midpoint=0, low="blue", mid="white",
high="red", space ="Lab" )

saveRDS(myData_SIP.RawCounts3_dots_logfc_both,"./data/myData_SIP.RawCounts3_dots_logfc_both.RDS")

######################################
######################################
## Refinedment

## Point plot remove low expression contidions
## Points_psedoCounts_SizeColor_Cytokines_Exp_AgingChange_removeSmall.pdf
myData_SIP.RawCounts3_dots_logfc_both_cytokine_bigExp<-myData_SIP.RawCounts3_dots_logfc_both[(myData_SIP.RawCounts3_dots_logfc_both$Var1.logfc %in% cytokines_pure) & (log2(myData_SIP.RawCounts3_dots_logfc_both$expYoung+1)>1),]
ggplot(myData_SIP.RawCounts3_dots_logfc_both_cytokine_bigExp, aes(Var1.logfc, CellTypes)) +
geom_point(aes(size = log2(expYoung+1), colour=logfc)) +scale_color_gradient2(midpoint=0, low="blue", mid="white",
high="red", space ="Lab" )

##### Do as Francesco suggested
##### heatmap by average expression
##### Heatmap_AverageExpression.pdf
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-myData_SIP.cytokines_Avg2_byAge_cytokinesNeri
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-as.data.frame(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter)
## Remove conditions where expression level < 1 in both young and old: make log2fc = 0
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`B-cell_young`[`B-cell_old`<1 & `B-cell_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`B-cell_old`[`B-cell_old`<1 & `B-cell_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`CD4-naive_young`[`CD4-naive_old`<1 & `CD4-naive_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`CD4-naive_old`[`CD4-naive_old`<1 & `CD4-naive_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`Th17_young`[`Th17_old`<1 & `Th17_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`Th17_old`[`Th17_old`<1 & `Th17_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`CD4-Cytotoxic_young`[`CD4-Cytotoxic_old`<1 & `CD4-Cytotoxic_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`CD4-Cytotoxic_old`[`CD4-Cytotoxic_old`<1 & `CD4-Cytotoxic_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`Plasma-cell_young`[`Plasma-cell_old`<1 & `Plasma-cell_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`Plasma-cell_old`[`Plasma-cell_old`<1 & `Plasma-cell_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`Myeloblast_young`[`Myeloblast_old`<1 & `Myeloblast_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`Myeloblast_old`[`Myeloblast_old`<1 & `Myeloblast_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`B-early_young`[`B-early_old`<1 & `B-early_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`B-early_old`[`B-early_old`<1 & `B-early_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`CD4-Gata3_young`[`CD4-Gata3_old`<1 & `CD4-Gata3_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`CD4-Gata3_old`[`CD4-Gata3_old`<1 & `CD4-Gata3_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`CD8_NK_young`[`CD8_NK_old`<1 & `CD8_NK_young`<1]<-0)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter<-within(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter,`CD8_NK_old`[`CD8_NK_old`<1 & `CD8_NK_young`<1]<-0)
##
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc<- as.data.frame(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter)
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`B-cell`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,1]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,2])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`CD4-naive`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,3]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,4])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`Th17`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,5]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,6])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`CD4-Cytotoxic`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,7]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,8])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`Plasma-cell`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,9]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,10])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`Myeloblast`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,11]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,12])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`B-early`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,13]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,14])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`CD4-Gata3`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,15]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,16])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc$`CD8_NK`<-log2(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,17]/myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,18])
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc<-myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,-c(1:18)]


myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[is.na(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc)]<-0
paletteLength <- 50
myBreaks <- c(seq(min(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc), 0, length.out=ceiling(paletteLength/2) + 1),
seq(max(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc)/paletteLength, max(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc), length.out=floor(paletteLength/2)))
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
## change cell type order
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc<-myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[,c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8_NK","Myeloblast","CD4-Gata3","Th17")]
## Heatmap_AverageExpression.pdf
pheatmap(t(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[rowSums(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc)!=0,]),angle_col = 90,scale = "none",color = myColor, breaks=myBreaks,cluster_rows = F)
###

#### Refinement on pseudo bulk RNA counts for cytokines 
myData_SIP.RawCounts3_dots_expOld<-myData_SIP.RawCounts3_dots[,grepl("old",colnames(myData_SIP.RawCounts3_dots))]
myData_SIP.RawCounts3_dots_expOld<-reshape2::melt(as.matrix(myData_SIP.RawCounts3_dots_expOld))
myData_SIP.RawCounts3_dots_logfc_both_2<-cbind(myData_SIP.RawCounts3_dots_logfc,myData_SIP.RawCounts3_dots_expYoung,myData_SIP.RawCounts3_dots_expOld)
myData_SIP.RawCounts3_dots_logfc_both_2$CellTypes<-str_split_fixed(myData_SIP.RawCounts3_dots_logfc_both_2$Var2,"#",2)[,2]
colnames(myData_SIP.RawCounts3_dots_logfc_both_2)<-c("gene","name_logfc","log2fc","gene_young","young","exp_pseudo_counts_young","gene_old","old","exp_pseudo_counts_old","CellTypes")
myData_SIP.RawCounts3_dots_logfc_both_2$mean_pseudo_counts<-rowMeans(myData_SIP.RawCounts3_dots_logfc_both_2[,c("exp_pseudo_counts_young","exp_pseudo_counts_old")])
>mean(myData_SIP.RawCounts3_dots_logfc_both_2$mean_pseudo_counts)
4.800734
## Points_psedoCounts_SizeColor_Cytokines_Exp_AgingChange_sizeMeanCounts.pdf
myData_SIP.RawCounts3_dots_logfc_both_2[myData_SIP.RawCounts3_dots_logfc_both_2$CellTypes=="CD8","CellTypes"]<-"CD8-NK"
myData_SIP.RawCounts3_dots_logfc_both_2$CellTypes<-factor(myData_SIP.RawCounts3_dots_logfc_both_2$CellTypes,levels = rev(c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8-NK","Myeloblast","CD4-Gata3","Th17")))
ggplot(myData_SIP.RawCounts3_dots_logfc_both_2[myData_SIP.RawCounts3_dots_logfc_both_2$gene %in% cytokines_pure,], aes(gene, CellTypes)) +
geom_point(aes(size = log2(mean_pseudo_counts+1), colour=log2fc)) +scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )

## filter by meanPseudoCounts by cutoff of mean of all genes, show them in heatmap
myData_SIP.RawCounts3_dots_logfc_both_2_filter<-myData_SIP.RawCounts3_dots_logfc_both_2
myData_SIP.RawCounts3_dots_logfc_both_2_filter[myData_SIP.RawCounts3_dots_logfc_both_2_filter$mean_pseudo_counts<4.800734,"log2fc"]<-0
myData_SIP.RawCounts3_dots_logfc_both_2_filter<-myData_SIP.RawCounts3_dots_logfc_both_2_filter[myData_SIP.RawCounts3_dots_logfc_both_2_filter$gene %in% cytokines_pure,]
myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide <- reshape2::dcast(myData_SIP.RawCounts3_dots_logfc_both_2_filter, gene ~ CellTypes, value.var="log2fc")
row.names(myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide)<-myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide$gene
myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide<-myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide[,-1]
myBreaks2 <- c(seq(min(myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide), 0, length.out=ceiling(paletteLength/2) + 1),
seq(max(myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide)/paletteLength, max(myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide), length.out=floor(paletteLength/2)))
myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide<-myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide[,c("B-early","B-cell","Plasma-cell","CD4-naive","CD4-Cytotoxic","CD8-NK","Myeloblast","CD4-Gata3","Th17")]
## Heatmap_psedoCounts_filteredBy4.8_log2fc.pdf
pheatmap(t(myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide[rowSums(myData_SIP.RawCounts3_dots_logfc_both_2_filter_wide)!=0,]),angle_col = 90,scale = "none",color = myColor, breaks=myBreaks2,cluster_rows = F)
save.image("./data/heatmap_Cytokines.RData")



## Heatmap_AverageExp_Cytokines_Modified_white0.31.pdf
pheatmap(t(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[rowSums(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc)!=0,]),angle_col = 90,scale = "row", breaks=c(seq(-3,-0.31,length.out = paletteLength/2),seq(0.31,3,length.out = paletteLength/2)),cluster_rows = F,color = colorRampPalette(c("blue", "white", "red"))(paletteLength-1))

## Heatmap_AverageExp_Cytokines_Modified_white0.5.pdf
pheatmap(t(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[rowSums(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc)!=0,]),angle_col = 90,scale = "row", breaks=c(seq(-3,-0.5,length.out = paletteLength/2),seq(0.5,3,length.out = paletteLength/2)),cluster_rows = F,color = colorRampPalette(c("blue", "white", "red"))(paletteLength-1))

## check conditions where cytokines are not detected
myData_SIP.cytokines_Avg2_byAge_cytokinesNeri[row.names(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[rowSums(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc)!=0,]),]==0.01

## Heatmap_AverageExpression_0.1_2.pdf
pheatmap(t(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[rowSums(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc)!=0,]),angle_col = 90,scale = "none",color = colorRampPalette(c("blue", "white", "red"))(paletteLength-1), breaks=c(seq(-2,-0.1,length.out = paletteLength/2),seq(0.1,2,length.out = paletteLength/2)),cluster_rows = F)

## Heatmap_AverageExpression_0.05_1.5.pdf
pheatmap(t(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc[rowSums(myData_SIP.cytokines_Avg2_byAge_cytokinesNeri_filter_logfc)!=0,]),angle_col = 90,scale = "none",color = colorRampPalette(c("blue", "white", "red"))(paletteLength-1), breaks=c(seq(-1.5,-0.05,length.out = paletteLength/2),seq(0.05,1.5,length.out = paletteLength/2)),cluster_rows = F)


=================================================================
## tSNE Markers Plot
# Bcell_Cd79a_Biton.pdf
FeaturePlot(myData_SIP, features = c("Cd79a"),pt.size = 0.02,coord.fixed = 1,ncol = 2,cols = c("grey","red"))
# BcellEarly_Cd79b_Biton.pdf           CD4_Cytotoxic_Ccl5_Alon.pdf          ILC2_Gata3_RD.pdf                    ILC3_Tnfsf11_OmidProvided.pdf
# BcellEarly_Cd79b_Biton_2.pdf         CD8NK_Gzma_BitonCD8_IdoNK.pdf        ILC3_Il22_BZ.pdf                     Myeloid_Ifitm3_BitonMonocypte.pdf
# Bcell_Cd79a_Biton.pdf                CD8NK_Xcl1_BitonNK.pdf               ILC3_Il22_RD.pdf                     Plasma_Apoe_none.pdf
# CD4Naive_Lef1_Alon.pdf               ILC3_Klrk1_OmidProvided.pdf          Plasma_Igha_none.pdf
## Celltypes_Markers_References_All.pdf
FeaturePlot(myData_SIP, features = c("Gzma","Cd79a","Cd79b","Ccl5","Lef1","Xcl1","Gata3","Il22","Klrk1","Tnfsf11","Ifitm3","Apoe","Igha"),pt.size = 0.02,coord.fixed = 1,ncol = 4,cols = c("grey","red"))
# Celltypes_Markers_References_All_ScaleLegend.pdf
FeaturePlot(myData_SIP, features = c("Gzma","Cd79a","Cd79b","Ccl5","Lef1","Xcl1","Gata3","Il22","Klrk1","Tnfsf11","Ifitm3","Apoe","Igha"),pt.size = 0.02,coord.fixed = 1,ncol = 4,cols = c("grey","red"),max.cutoff = 4)