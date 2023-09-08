## Author: Jing
setwd("./")
immune.combined<-readRDS("./data/immune.combined.rds")
library(Seurat)
## SingleIfng_umap_treatment_celltype_lableF.pdf
DimPlot(immune.combined, reduction = "umap", label = TRUE,split.by ="treatment")
# rename celltype and show in umap
celltype_new<-c("TA","Enterocyte-Prog","Enterocyte","Stem","Enteroendocrine","Goblet-Paneth","Tuft")
immune.combined$celltype_new<-celltype_new[immune.combined$seurat_clusters]
Idents(immune.combined) <- "celltype_new"
DimPlot(immune.combined, reduction = "umap", label = TRUE,split.by ="treatment")+coord_fixed(ratio=1)
my.data<-FetchData(immune.combined,c("orig.ident","nCount_RNA","nFeature_RNA","seurat_clusters","celltype","treatment","Cd74_group","celltype.treatment","celltype.treatment.Cd74","celltype_new"))
umap<-Embeddings(immune.combined, reduction = "umap")[, 1:2]
my.data<-merge(my.data,umap,by="row.names")
rownames(my.data)<-my.data$Row.names
my.data<-my.data[,-1]
library(dplyr)
library(ggpubr)
library(stringr)
# add sampleID
my.data$sampleID<-str_split_fixed(row.names(my.data),"\\.",2)[,2]
# ratio2 is same with ratio so just keep one, both are based on samples
#data.ratio <- my.data% >% dplyr::group_by(sampleID,celltype_new) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count)) %>% dplyr::ungroup() %>% dplyr::group_by(sampleID) %>% dplyr::mutate(ratio2=count/sum(count))
#data.ratio <- my.data %>% dplyr::group_by(sampleID,celltype_new) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count)) %>% dplyr::ungroup() %>% dplyr::group_by(sampleID) %>% dplyr::mutate(ratio2=count/sum(count))
sum(data.ratio$count)
226/13054
data.ratio <- my.data %>% dplyr::group_by(sampleID,celltype_new) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
sum(data.ratio[data.ratio$sampleID=="1_shash3","count"])
226/932
264/932
data.ratio$hashtag<-str_split_fixed(data.ratio$sampleID,"_",2)[,2]
# add treatment
data.ratio$hashtag<-str_split_fixed(data.ratio$sampleID,"_",2)[,2]
data.ratio$treatment="Control"
data.ratio[(data.ratio$hashtag=="shash5"|data.ratio$hashtag=="shash9"),]$treatment="Ifng"
# have a check
unique(data.ratio[data.ratio$hashtag=="shash3","treatment"])
unique(data.ratio[data.ratio$hashtag=="shash5","treatment"])
unique(data.ratio[data.ratio$hashtag=="shash7","treatment"])
unique(data.ratio[data.ratio$hashtag=="shash9","treatment"])
## SingleIfng_CellCompostiton_BoxPlot.pdf
ggbarplot(data.ratio, x = "celltype_new", y = "ratio", add = "mean_se", color = "treatment", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = treatment))
# have a check
(data.ratio[data.ratio$treatment=="Control"&data.ratio$celltype_new=="Stem",])
(data.ratio[data.ratio$treatment=="Control"&data.ratio$celltype_new=="Stem","ratio"])
mean(c(0.149,0.145,0.112,0.129))

## test wilcox and t test with paired and nonpaired
data.ratio$celltype_new<-factor(data.ratio$celltype_new,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
ggbarplot(data.ratio, x = "celltype_new", y = "ratio", add = "mean_sd", color = "treatment", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = treatment))
ggbarplot(data.ratio, x = "celltype_new", y = "ratio", add = "mean_sd", color = "treatment", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = treatment),method="t.test",paired=T,label = "p.format")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))
ggbarplot(data.ratio, x = "celltype_new", y = "ratio", add = "mean_sd", color = "treatment", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = treatment),method="t.test",paired=F,label = "p.format")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))
# final
ggbarplot(data.ratio, x = "celltype_new", y = "ratio", add = "mean_sd", color = "treatment", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = treatment),method="t.test",paired=F,label = "p.format",label.y = 0.4)+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))+scale_color_manual(values = c("grey","grey0"))+ggtitle("t.test_pairedF")+xlab("")+ylab("Cell ratio(%)")
write.csv(data.ratio,"./data_output/data.ratio.csv")
## mhc2 expression check
mhc2<-c("H2-K1",
"H2-Ke6",
"H2-Oa",
"H2-DMa",
"H2-DMb2",
"H2-DMb1",
"H2-Ob",
"H2-Ab1",
"H2-Aa",
"H2-Eb1",
"H2-Eb2",
"H2-D1",
"H2-Q1",
"H2-Q2",
"H2-Q4",
"H2-Q6",
"H2-Q7",
"H2-Q10",
"H2-T24",
"H2-T23",
"H2-T22",
"H2-T3",
"H2-M10.2",
"H2-M10.1",
"H2-M10.3",
"H2-M10.4",
"H2-M11",
"H2-M9",
"H2-M1",
"H2-M10.5",
"H2-M10.6",
"H2-M5",
"H2-M3",
"H2-M2",
"Cd74")

# do heatmap choos assay='RNA' (defalut) not 'integrated', sacale all genes first
immune.combined2 <- ScaleData(object = immune.combined, features = rownames(immune.combined))
# check gene names
mhc2 %in% rownames(immune.combined)
mhc2 %in% rownames(immune.combined2)
length(rownames(immune.combined))
grep("^H2",rownames(immune.combined))
row.names(immune.combined)[grep("^H2",rownames(immune.combined))]
mhc2<-gsub("-","\\.",mhc2)
# for heatmap only immune.combined2 can work, while for vln, both immune.combined and immune.combined2 can work
#DoHeatmap(immune.combined, features = mhc2, size = 3,assay='RNA',group.by = "celltype_new")
DoHeatmap(immune.combined2, features = mhc2, size = 3,assay='RNA',group.by = "celltype_new")
#plots <- VlnPlot(immune.combined, features = mhc2, split.by = "treatment", group.by = "celltype_new", pt.size = 0, combine = FALSE)
plots <- VlnPlot(immune.combined2, features = mhc2, split.by = "treatment", group.by = "celltype_new", pt.size = 0, combine = FALSE)
patchwork::wrap_plots(plots = plots, ncol = 1)
VlnPlot(immune.combined2, features = "Cd74", split.by = "treatment", group.by = "celltype_new", pt.size = 0, combine = FALSE)
VlnPlot(immune.combined2, features = "Cd74", split.by = "treatment", group.by = "celltype_new", pt.size = 0,cols=c("grey", "grey0"))
VlnPlot(immune.combined2, features = "Cd74", split.by = "treatment", group.by = "celltype_new", pt.size = 0) + ggplot2::scale_fill_manual(values = c('grey', 'grey0')) # grey could be replaced by gray
immune.combined2$celltype_new<-factor(immune.combined2$celltype_new,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
for (i in mhc2$genename2){pdf(paste(i,".pdf",sep=""),width = 5,height = 5);try(print(VlnPlot(immune.combined2, features = i, split.by = "treatment", group.by = "celltype_new", pt.size = 0) + ggplot2::scale_fill_manual(values = c('gray', 'gray0'))));dev.off()}


## average expression
mhc2<-read.csv("./data/KEGG_mhc.genes",header = T)
markers10<-read.table("./data/10_top_markers_lnx.txt",sep="\t",header = T)
library(pheatmap)
immune.combined2$celltype_new.treatment<-paste(immune.combined2$celltype_new,immune.combined2$treatment,sep=".")
Idents(immune.combined2)<-immune.combined2$celltype_new.treatment
cluster.averages <- AverageExpression(immune.combined2,assays = 'RNA')
cluster.averages.RNA<-(cluster.averages$RNA)
mhc2$genename2<-gsub("-","\\.",mhc2$genename)
cluster.averages.RNA.mhc2<-(cluster.averages.RNA[as.character(mhc2$genename2),])

#ll<-c("Stem.sc_control","TA.sc_control","Enterocyte-Prog.sc_control","Enterocyte.sc_control","Enteroendocrine.sc_control","Goblet-Paneth.sc_control","Tuft.sc_control","Stem.sc_ifnr","TA.sc_ifnr","Enterocyte-Prog.sc_ifnr","Enterocyte.sc_ifnr","Enteroendocrine.sc_ifnr","Goblet-Paneth.sc_ifnr","Tuft.sc_ifnr")
ll2<-c("Stem.sc_control","Stem.sc_ifnr","TA.sc_control","TA.sc_ifnr","Enterocyte-Prog.sc_control","Enterocyte-Prog.sc_ifnr","Enterocyte.sc_control","Enterocyte.sc_ifnr","Enteroendocrine.sc_control","Enteroendocrine.sc_ifnr","Goblet-Paneth.sc_control","Goblet-Paneth.sc_ifnr","Tuft.sc_control","Tuft.sc_ifnr")


## SingleIfng_Mhc2_averageExpression_scaleRow_rankCol.pdf
cluster.averages.RNA.mhc2<-cluster.averages.RNA.mhc2[,ll2]
pheatmap(na.omit(cluster.averages.RNA.mhc2[apply(cluster.averages.RNA.mhc2,1,sd)!=0,]),cluster_rows = T,cluster_cols = F,fontsize_row= 6,scale="row")

cluster.averages.RNA.markers10<-(cluster.averages.RNA[as.character(markers10$gene),])
row.names(cluster.averages.RNA.markers10)<-paste(row.names(cluster.averages.RNA.markers10),markers10[match(row.names(cluster.averages.RNA.markers10),markers10$gene),"cluster"],sep=".")
cluster.averages.RNA.markers10<-cluster.averages.RNA.markers10[,ll2]
pheatmap(na.omit(cluster.averages.RNA.markers10[apply(cluster.averages.RNA.markers10,1,sd)!=0,]),cluster_rows = T,cluster_cols = F,fontsize_row= 3,scale="row")
## DEG
celltype_new<-unique(immune.combined2$celltype_new)
for(i in celltype_new){assign(paste("treatment.DEG.wilcox",i,sep = "_") , FindMarkers(immune.combined2, ident.1 = paste(i,"sc_ifnr",sep="."), ident.2 = paste(i,"sc_control",sep="."), verbose = FALSE,group.by='celltype_new.treatment'))}
for(i in celltype_new){write.csv(get(paste("treatment.DEG.wilcox",i,sep = "_")),paste(i,"treatment.DEG.wilcox.csv",sep = "_"))}

library(ggrepel)
library(RColorBrewer)

volcano<-function(data,name,cf=0.5) {
data$sig[(data$p_val_adj > 0.05|data$p_val_adj=="NA")|(data$avg_logFC < cf)& data$avg_logFC > -cf] <- "no"
data$sig[data$p_val_adj <= 0.05 & data$avg_logFC >= cf] <- "up"
data$sig[data$p_val_adj <= 0.05 & data$avg_logFC <= -cf] <- "down"
p<-ggplot(data[!is.na(data$sig),],aes(avg_logFC,-1*log10(p_val_adj),color = sig))+geom_point(size=2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(name)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
print(p)
}
volcano(treatment.DEG.wilcox_TA,"TA")
pdf("./figure/volcano.celltypes.pdf")
for(i in celltype_new){volcano(get(paste("treatment.DEG.wilcox",i,sep = "_")),i)}
dev.off()

volcano_label<-function(data,clname,cf=0.5) {
data$sig[(data$p_val_adj > 0.05|data$p_val_adj=="NA")|(data$avg_logFC < cf)& data$avg_logFC > -cf] <- "no"
data$sig[data$p_val_adj <= 0.05 & data$avg_logFC >= cf] <- "up"
data$sig[data$p_val_adj <= 0.05 & data$avg_logFC <= -cf] <- "down"
data$sign<-NA
data[data$sig!="no","sign"] <- row.names(data[data$sig!="no",])
p<-ggplot(data[!is.na(data$sig),],aes(avg_logFC,-1*log10(p_val_adj),color = sig))+geom_point(size=0.2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(clname)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)
print(p)
}
volcano_label(treatment.DEG.wilcox_TA,"TA")

pdf("./figure/volcano.celltypes.label.pdf")
for(i in celltype_new){volcano_label(get(paste("treatment.DEG.wilcox",i,sep = "_")),i)}
dev.off()


## venn on Ifng and Aging
### deg vene
# test
deg_aging<-read.csv("./data_output/deg_dresden/deg_wilcox_lfc0.1_DefaulstAssayRNA/Stem_age.DEG.wilcox.csv",header = T,stringsAsFactors = F)
deg_ifng<-read.csv("./data_output/deg_ifng/Stem_treatment_Ifng.DEG.wilcox.csv",header = T,stringsAsFactors = F)
deg_ifng_sig<-deg_ifng[deg_ifng$p_val_adj<0.05&abs(deg_ifng$avg_logFC)>0.5,]
deg_aging_sig<-deg_aging[deg_aging$p_val_adj<0.05,]
require(VennDiagram)
venn.diagram(list(deg_Ifng_sig_5=deg_ifng_sig$X,deg_aging_sig=deg_aging_sig$X),filename = "test.pdf",col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
deg_ifng_aging<-merge(deg_ifng_sig,deg_aging_sig,by.x="X",by.y="X",all=T,suffixes = c(".ifng",".aging"))
write.csv(deg_ifng_aging,"./data_output/deg_ifng_aging.csv")

venn_deg<-function(fn_ifng,fn_aging,fn_out){
deg_aging<-read.csv(paste("./data_output/deg_dresden/deg_wilcox_lfc0.1_DefaulstAssayRNA/",fn_aging,sep=""),header = T,stringsAsFactors = F)
deg_ifng<-read.csv(paste("./data_output/deg_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F)
deg_ifng_sig<-deg_ifng[deg_ifng$p_val_adj<0.05&abs(deg_ifng$avg_logFC)>0.5,]
deg_aging_sig<-deg_aging[deg_aging$p_val_adj<0.05,]
require(VennDiagram)
venn.diagram(list(deg_Ifng_sig_5=deg_ifng_sig$X,deg_aging_sig=deg_aging_sig$X),filename = paste(fn_out,"pdf",sep="."),col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
deg_ifng_aging<-merge(deg_ifng_sig,deg_aging_sig,by.x="X",by.y="X",all=T,suffixes = c(".ifng",".aging"))
write.csv(deg_ifng_aging,paste("./data_output/",fn_out,sep=""))
}
venn_deg("Stem_treatment_Ifng.DEG.wilcox.csv","Stem_age.DEG.wilcox.csv","venn_stem_deg")

##ipa deg
# test

ipa_aging<-read.table("./data_output/ipa_dresden/Stem_Aging.txt.2d",header = T,stringsAsFactors = F,sep = "\t")
ipa_ifng<-read.table("./data_output/ipa_ifng/Stem_Ifng.txt.2d",header = T,stringsAsFactors = F,sep = "\t")

venn_ipa<-function(fn_ifng,fn_aging,fn_out){
ipa_aging<-read.table(paste("./data_output/ipa_dresden/",fn_aging,sep=""),header = T,stringsAsFactors = F,sep="\t")
ipa_ifng<-read.table(paste("./data_output/ipa_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F,sep="\t")
ipa_ifng_sig<-ipa_ifng[ipa_ifng$X.log.p.value.>1.3,]
ipa_aging_sig<-ipa_aging[ipa_aging$X.log.p.value.>1.3,]
require(VennDiagram)
venn.diagram(list(ipa_Ifng_sig_5=na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),ipa_aging_sig=na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)),filename = paste(fn_out,"pdf",sep="."),col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
ipa_ifng_aging<-merge(ipa_ifng_sig,ipa_aging_sig,by.x="Ingenuity.Canonical.Pathways",by.y="Ingenuity.Canonical.Pathways",all=T,suffixes = c(".ifng",".aging"))
write.csv(ipa_ifng_aging,paste("./data_output/",fn_out,".csv",sep=""))
}
venn_ipa("Stem_Ifng.txt.2d","Stem_Aging.txt.2d","venn_stem_ipa")
###
venn_ipa("Stem_Ifng.txt.2d","Stem_Aging.txt.2d","venn_stem_ipa")
venn_ipa("TA_Ifng.txt.2d","TA_Aging.txt.2d","venn_TA_ipa")
venn_ipa("Enterocyte-Prog_Ifng.txt.2d","Early-Enterocytes-Prog_Aging.txt.2d","venn_Enterocyte-Prog_ipa")
venn_ipa("Enterocyte_Ifng.txt.2d","Enterocytes_Aging.txt.2d","venn_Enterocyte_ipa")
venn_ipa("Enteroendocrine_Ifng.txt.2d","Enteroendocrine_Aging.txt.2d","venn_Enteroendocrine_ipa")
venn_ipa("Goble-Paneth_Ifng.txt.2d","Goblet-Paneth_Aging.txt.unique.2d","venn_Goblet-Paneth_ipa")
venn_ipa("Tuft_Ifng.txt.2d","Tuft_Aging.txt.2d","venn_Tuft_ipa")
venn_deg("Stem_treatment_Ifng.DEG.wilcox.csv","Stem_age.DEG.wilcox.csv","venn_stem_deg")
venn_deg("TA_treatment_Ifng.DEG.wilcox.csv","TA_age.DEG.wilcox.csv","venn_TA_deg")
venn_deg("Enterocyte-Prog_treatment_Ifng.DEG.wilcox.csv","Early.Enterocytes.Progenitors_age.DEG.wilcox.csv","venn_Enterocyte-Prog_deg")
venn_deg("Enterocyte_treatment_Ifng.DEG.wilcox.csv","Enterocyte_age.DEG.wilcox.csv","venn_Enterocyte_deg")
venn_deg("Enteroendocrine_treatment_Ifng.DEG.wilcox.csv","Enteroendocrine_age.DEG.wilcox.csv","venn_Enteroendocrine_deg")
venn_deg("Goblet-Paneth_treatment_Ifng.DEG.wilcox.csv","Goblet-Paneth_age.DEG.wilcox.unique.csv","venn_Goblet-Paneth_deg")
venn_deg("Tuft_treatment_Ifng.DEG.wilcox.csv","Tuft_age.DEG.wilcox.csv","venn_Tuft_deg")

## venn pdf
venn_ipa_2<-function(fn_ifng,fn_aging,fn_out){
ipa_aging<-read.table(paste("./data_output/ipa_dresden/",fn_aging,sep=""),header = T,stringsAsFactors = F,sep="\t")
ipa_ifng<-read.table(paste("./data_output/ipa_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F,sep="\t")
ipa_ifng_sig<-ipa_ifng[ipa_ifng$X.log.p.value.>1.3,]
ipa_aging_sig<-ipa_aging[ipa_aging$X.log.p.value.>1.3,]
require(VennDiagram)
temp<-venn.diagram(list(ipa_Ifng_sig_5=na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),ipa_aging_sig=na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
pdf(file=paste(fn_out,"pdf",sep="."))
grid.draw(temp)
dev.off()
}

venn_deg_2<-function(fn_ifng,fn_aging,fn_out){
deg_aging<-read.csv(paste("./data_output/deg_dresden/deg_wilcox_lfc0.1_DefaulstAssayRNA/",fn_aging,sep=""),header = T,stringsAsFactors = F)
deg_ifng<-read.csv(paste("./data_output/deg_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F)
deg_ifng_sig<-deg_ifng[deg_ifng$p_val_adj<0.05&abs(deg_ifng$avg_logFC)>0.5,]
deg_aging_sig<-deg_aging[deg_aging$p_val_adj<0.05,]
require(VennDiagram)
temp<-venn.diagram(list(deg_Ifng_sig_5=deg_ifng_sig$X,deg_aging_sig=deg_aging_sig$X),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
pdf(file=paste(fn_out,"pdf",sep="."))
grid.draw(temp)
dev.off()
}

venn_ipa_2("Stem_Ifng.txt.2d","Stem_Aging.txt.2d","venn_stem_ipa")
venn_ipa_2("TA_Ifng.txt.2d","TA_Aging.txt.2d","venn_TA_ipa")
venn_ipa_2("Enterocyte-Prog_Ifng.txt.2d","Early-Enterocytes-Prog_Aging.txt.2d","venn_Enterocyte-Prog_ipa")
venn_ipa_2("Enterocyte_Ifng.txt.2d","Enterocytes_Aging.txt.2d","venn_Enterocyte_ipa")
venn_ipa_2("Enteroendocrine_Ifng.txt.2d","Enteroendocrine_Aging.txt.2d","venn_Enteroendocrine_ipa")
venn_ipa_2("Goble-Paneth_Ifng.txt.2d","Goblet-Paneth_Aging.txt.unique.2d","venn_Goblet-Paneth_ipa")
venn_ipa_2("Tuft_Ifng.txt.2d","Tuft_Aging.txt.2d","venn_Tuft_ipa")
venn_deg_2("Stem_treatment_Ifng.DEG.wilcox.csv","Stem_age.DEG.wilcox.csv","venn_stem_deg")
venn_deg_2("TA_treatment_Ifng.DEG.wilcox.csv","TA_age.DEG.wilcox.csv","venn_TA_deg")
venn_deg_2("Enterocyte-Prog_treatment_Ifng.DEG.wilcox.csv","Early.Enterocytes.Progenitors_age.DEG.wilcox.csv","venn_Enterocyte-Prog_deg")
venn_deg_2("Enterocyte_treatment_Ifng.DEG.wilcox.csv","Enterocyte_age.DEG.wilcox.csv","venn_Enterocyte_deg")
venn_deg_2("Enteroendocrine_treatment_Ifng.DEG.wilcox.csv","Enteroendocrine_age.DEG.wilcox.csv","venn_Enteroendocrine_deg")
venn_deg_2("Goblet-Paneth_treatment_Ifng.DEG.wilcox.csv","Goblet-Paneth_age.DEG.wilcox.unique.csv","venn_Goblet-Paneth_deg")
venn_deg_2("Tuft_treatment_Ifng.DEG.wilcox.csv","Tuft_age.DEG.wilcox.csv","venn_Tuft_deg")

## apoptosis
apoptosis<-read.table("./data/KEGG_Apoptosis_mouse_gene.txt",sep="\t",header = T)
cluster.averages.RNA.apoptosis<-(cluster.averages.RNA[as.character(apoptosis$genename),])
cluster.averages.RNA.apoptosis<-cluster.averages.RNA.apoptosis[,ll2]
pheatmap(na.omit(cluster.averages.RNA.apoptosis[apply(cluster.averages.RNA.apoptosis,1,sd)!=0,]),cluster_rows = T,cluster_cols = F,fontsize_row= 6,scale="row")

## venn adjusted
venn_ipa_2<-function(fn_ifng,fn_aging,fn_out){
    ipa_aging<-read.table(paste("./data_output/ipa_dresden/",fn_aging,sep=""),header = T,stringsAsFactors = F,sep="\t")
    ipa_ifng<-read.table(paste("./data_output/ipa_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F,sep="\t")
    ipa_ifng_sig<-ipa_ifng[ipa_ifng$X.log.p.value.>1.3,]
    ipa_aging_sig<-ipa_aging[ipa_aging$X.log.p.value.>1.3,]
    require(VennDiagram)
    temp<-venn.diagram(list(ipa_Ifng_sig_5=na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),ipa_aging_sig=na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test=T,total.population = (nrow(ipa_aging)+nrow(ipa_ifng)))
    pdf(file=paste(fn_out,"pdf",sep="."))
    grid.draw(temp)
    dev.off()
}

venn_deg_2<-function(fn_ifng,fn_aging,fn_out){
    deg_aging<-read.csv(paste("./data_output/deg_dresden/deg_wilcox_lfc0.1_DefaulstAssayRNA/",fn_aging,sep=""),header = T,stringsAsFactors = F)
    deg_ifng<-read.csv(paste("./data_output/deg_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F)
    deg_ifng_sig<-deg_ifng[deg_ifng$p_val_adj<0.05&abs(deg_ifng$avg_logFC)>0.5,]
    deg_aging_sig<-deg_aging[deg_aging$p_val_adj<0.05,]
    require(VennDiagram)
    temp<-venn.diagram(list(deg_Ifng_sig_5=deg_ifng_sig$X,deg_aging_sig=deg_aging_sig$X),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test=T,total.population = (nrow(deg_aging)+nrow(deg_ifnr)))
    pdf(file=paste(fn_out,"pdf",sep="."))
    grid.draw(temp)
    dev.off()
}
venn_ipa_2("Stem_Ifng.txt.2d","Stem_Aging.txt.2d","venn_stem_ipa")
venn_ipa_2("TA_Ifng.txt.2d","TA_Aging.txt.2d","venn_TA_ipa")
venn_ipa_2("Enterocyte-Prog_Ifng.txt.2d","Early-Enterocytes-Prog_Aging.txt.2d","venn_Enterocyte-Prog_ipa")
venn_ipa_2("Enterocyte_Ifng.txt.2d","Enterocytes_Aging.txt.2d","venn_Enterocyte_ipa")
venn_ipa_2("Enteroendocrine_Ifng.txt.2d","Enteroendocrine_Aging.txt.2d","venn_Enteroendocrine_ipa")
venn_ipa_2("Goble-Paneth_Ifng.txt.2d","Goblet-Paneth_Aging.txt.unique.2d","venn_Goblet-Paneth_ipa")
venn_ipa_2("Tuft_Ifng.txt.2d","Tuft_Aging.txt.2d","venn_Tuft_ipa")
venn_deg_2("Stem_treatment_Ifng.DEG.wilcox.csv","Stem_age.DEG.wilcox.csv","venn_stem_deg")
venn_deg_2("TA_treatment_Ifng.DEG.wilcox.csv","TA_age.DEG.wilcox.csv","venn_TA_deg")
venn_deg_2("Enterocyte-Prog_treatment_Ifng.DEG.wilcox.csv","Early.Enterocytes.Progenitors_age.DEG.wilcox.csv","venn_Enterocyte-Prog_deg")
venn_deg_2("Enterocyte_treatment_Ifng.DEG.wilcox.csv","Enterocyte_age.DEG.wilcox.csv","venn_Enterocyte_deg")
venn_deg_2("Enteroendocrine_treatment_Ifng.DEG.wilcox.csv","Enteroendocrine_age.DEG.wilcox.csv","venn_Enteroendocrine_deg")
venn_deg_2("Goblet-Paneth_treatment_Ifng.DEG.wilcox.csv","Goblet-Paneth_age.DEG.wilcox.unique.csv","venn_Goblet-Paneth_deg")
venn_deg_2("Tuft_treatment_Ifng.DEG.wilcox.csv","Tuft_age.DEG.wilcox.csv","venn_Tuft_deg")

venn_deg_3<-function(fn_ifng,fn_aging,fn_out){
    deg_aging<-read.csv(paste("./data_output/deg_dresden/deg_wilcox_lfc0.1_DefaulstAssayRNA/",fn_aging,sep=""),header = T,stringsAsFactors = F)
    deg_ifng<-read.csv(paste("./data_output/deg_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F)
    deg_ifng_sig<-deg_ifng[deg_ifng$p_val_adj<0.05&abs(deg_ifng$avg_logFC)>0.5,]
    deg_aging_sig<-deg_aging[deg_aging$p_val_adj<0.05,]
    require(VennDiagram)
    temp<-venn.diagram(list(deg_Ifng_sig_5=deg_ifng_sig$X,deg_aging_sig=deg_aging_sig$X),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test=T,total.population = (nrow(deg_aging)+nrow(deg_ifnr)),inverted=T)
    pdf(file=paste(fn_out,"pdf",sep="."))
    grid.draw(temp)
    dev.off()
}
venn_ipa_3<-function(fn_ifng,fn_aging,fn_out){
    ipa_aging<-read.table(paste("./data_output/ipa_dresden/",fn_aging,sep=""),header = T,stringsAsFactors = F,sep="\t")
    ipa_ifng<-read.table(paste("./data_output/ipa_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F,sep="\t")
    ipa_ifng_sig<-ipa_ifng[ipa_ifng$X.log.p.value.>1.3,]
    ipa_aging_sig<-ipa_aging[ipa_aging$X.log.p.value.>1.3,]
    require(VennDiagram)
    temp<-venn.diagram(list(ipa_Ifng_sig_5=na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),ipa_aging_sig=na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "Oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test=T,total.population = (nrow(ipa_aging)+nrow(ipa_ifng)),inverted=T)
    pdf(file=paste(fn_out,"pdf",sep="."))
    grid.draw(temp)
    dev.off()
}
venn_deg_3("Enterocyte-Prog_treatment_Ifng.DEG.wilcox.csv","Early.Enterocytes.Progenitors_age.DEG.wilcox.csv","venn_Enterocyte-Prog_deg3")
venn_ipa_3("Stem_Ifng.txt.2d","Stem_Aging.txt.2d","venn_stem_ipa3")
venn_deg_3("Stem_treatment_Ifng.DEG.wilcox.csv","Stem_age.DEG.wilcox.csv","venn_stem_deg3")
venn_deg_3("TA_treatment_Ifng.DEG.wilcox.csv","TA_age.DEG.wilcox.csv","venn_TA_deg3")



venn_ipa_p<-function(fn_ifng,fn_aging,fn_out){
    ipa_aging<-read.table(paste("./data_output/ipa_dresden/",fn_aging,sep=""),header = T,stringsAsFactors = F,sep="\t")
    ipa_ifng<-read.table(paste("./data_output/ipa_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F,sep="\t")
    ipa_ifng_sig<-ipa_ifng[ipa_ifng$X.log.p.value.>1.3,]
    ipa_aging_sig<-ipa_aging[ipa_aging$X.log.p.value.>1.3,]
    m=as.numeric(nrow(ipa_aging))
    n=as.numeric(nrow(ipa_ifng))
    q=as.numeric(length((intersect(na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)))))
    k=as.numeric(nrow(ipa_ifng_sig))
    print(paste(q,m,n,k,sep=","))
    print(paste(fn_out,"is:",sep=" "))
    print(phyper(q, m, n, k))
}
venn_ipa_p("Stem_Ifng.txt.2d","Stem_Aging.txt.2d","venn_stem_ipa")
venn_ipa_p("TA_Ifng.txt.2d","TA_Aging.txt.2d","venn_TA_ipa")
venn_ipa_p("Enterocyte-Prog_Ifng.txt.2d","Early-Enterocytes-Prog_Aging.txt.2d","venn_Enterocyte-Prog_ipa")
venn_ipa_p("Enterocyte_Ifng.txt.2d","Enterocytes_Aging.txt.2d","venn_Enterocyte_ipa")
venn_ipa_p("Enteroendocrine_Ifng.txt.2d","Enteroendocrine_Aging.txt.2d","venn_Enteroendocrine_ipa")
venn_ipa_p("Goble-Paneth_Ifng.txt.2d","Goblet-Paneth_Aging.txt.unique.2d","venn_Goblet-Paneth_ipa")
venn_ipa_p("Tuft_Ifng.txt.2d","Tuft_Aging.txt.2d","venn_Tuft_ipa")

venn_deg_p<-function(fn_ifng,fn_aging,fn_out){
    deg_aging<-read.csv(paste("./data_output/deg_dresden/deg_wilcox_lfc0.1_DefaulstAssayRNA/",fn_aging,sep=""),header = T,stringsAsFactors = F)
    deg_ifng<-read.csv(paste("./data_output/deg_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F)
    deg_ifng_sig<-deg_ifng[deg_ifng$p_val_adj<0.05&abs(deg_ifng$avg_logFC)>0.5,]
    deg_aging_sig<-deg_aging[deg_aging$p_val_adj<0.05,]
    m=as.numeric(nrow(deg_aging))
    n=as.numeric(nrow(deg_ifng))
    q=as.numeric(length((intersect(na.omit(deg_ifng_sig$X),na.omit(deg_aging_sig$X)))))
    k=as.numeric(nrow(deg_ifng_sig))
    print(paste(q,m,n,k,sep=","))
    print(paste(fn_out,"is:",sep=" "))
    print(phyper(q, m, n, k))
}
venn_deg_p("Stem_treatment_Ifng.DEG.wilcox.csv","Stem_age.DEG.wilcox.csv","venn_stem_deg")
venn_deg_p("TA_treatment_Ifng.DEG.wilcox.csv","TA_age.DEG.wilcox.csv","venn_TA_deg")
venn_deg_p("Enterocyte-Prog_treatment_Ifng.DEG.wilcox.csv","Early.Enterocytes.Progenitors_age.DEG.wilcox.csv","venn_Enterocyte-Prog_deg")
venn_deg_p("Enterocyte_treatment_Ifng.DEG.wilcox.csv","Enterocyte_age.DEG.wilcox.csv","venn_Enterocyte_deg")
venn_deg_p("Enteroendocrine_treatment_Ifng.DEG.wilcox.csv","Enteroendocrine_age.DEG.wilcox.csv","venn_Enteroendocrine_deg")
venn_deg_p("Goblet-Paneth_treatment_Ifng.DEG.wilcox.csv","Goblet-Paneth_age.DEG.wilcox.unique.csv","venn_Goblet-Paneth_deg")
venn_deg_p("Tuft_treatment_Ifng.DEG.wilcox.csv","Tuft_age.DEG.wilcox.csv","venn_Tuft_deg")

==============================
==============================
==============================
## redo
# ipa
venn_ipa_p<-function(fn_ifng,fn_aging,fn_out){
    ipa_aging<-read.delim(paste("./data_output/ipa_dresden/",fn_aging,sep=""),header = T,stringsAsFactors = F,sep="\t")
    ipa_ifng<-read.delim(paste("./data_output/ipa_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F,sep="\t")
    ipa_aging_sig<-ipa_aging[ipa_aging$X.log.p.value.>1.3,]
    ipa_ifng_sig<-ipa_ifng[ipa_ifng$X.log.p.value.>1.3,]
    m=as.numeric(nrow(ipa_aging))
    n=as.numeric(nrow(ipa_ifng))
    q=as.numeric(length((intersect(na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)))))
    k=as.numeric(nrow(ipa_ifng_sig))
    print(paste("q= overlap in sigs","m= all aging","n= all ifng","k= sig ifng",sep=","))
    print(paste(q,m,n,k,sep=","))
    print(paste(fn_out,"is:",sep=" "))
    print(phyper(q, m, n, k))
    pt <-phyper(q, m, n, k)
    require(VennDiagram)
    temp<-venn.diagram(list(ipa_Ifng_sig_5=na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),ipa_aging_sig=na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = paste("Oral pvalue=",pt,sep=" "),height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
    pdf(file=paste(fn_out,"pdf",sep="."))
    grid.draw(temp)
    dev.off()
    ipa_ifng_aging<-merge(ipa_ifng_sig,ipa_aging_sig,by.x="Ingenuity.Canonical.Pathways",by.y="Ingenuity.Canonical.Pathways",all=T,suffixes = c(".ifng",".aging"))
    write.csv(ipa_ifng_aging,paste("./data_output/",fn_out,".redo.csv",sep=""))   
}
venn_ipa_p("Stem_Ifng.txt.2d","Stem_Aging.txt.2d","venn_stem_ipa")
venn_ipa_p("TA_Ifng.txt.2d","TA_Aging.txt.2d","venn_TA_ipa")
venn_ipa_p("Enterocyte-Prog_Ifng.txt.2d","Early-Enterocytes-Prog_Aging.txt.2d","venn_Enterocyte-Prog_ipa")
venn_ipa_p("Enterocyte_Ifng.txt.2d","Enterocytes_Aging.txt.2d","venn_Enterocyte_ipa")
venn_ipa_p("Enteroendocrine_Ifng.txt.2d","Enteroendocrine_Aging.txt.2d","venn_Enteroendocrine_ipa")
venn_ipa_p("Goble-Paneth_Ifng.txt.2d","Goblet-Paneth_Aging.txt.unique.2d","venn_Goblet-Paneth_ipa")
venn_ipa_p("Tuft_Ifng.txt.2d","Tuft_Aging.txt.2d","venn_Tuft_ipa")


# deg
venn_deg_p<-function(fn_ifng,fn_aging,fn_out){
    deg_aging<-read.csv(paste("./data_output/deg_dresden/deg_wilcox_lfc0.1_DefaulstAssayRNA/",fn_aging,sep=""),header = T,stringsAsFactors = F)
    deg_ifng<-read.csv(paste("./data_output/deg_ifng/",fn_ifng,sep=""),header = T,stringsAsFactors = F)
    deg_aging_sig<-deg_aging[deg_aging$p_val_adj<0.05,]
    deg_ifng_sig<-deg_ifng[deg_ifng$p_val_adj<0.05&abs(deg_ifng$avg_logFC)>0.5,]  
    deg_aging_sig$d<-ifelse(deg_aging_sig$avg_logFC>0,"up","down")
    deg_ifng_sig$d<-ifelse(deg_ifng_sig$avg_logFC>0,"up","down")
    deg_aging_sig$X_d<-paste(deg_aging_sig$X,deg_aging_sig$d,sep="_")
    deg_ifng_sig$X_d<-paste(deg_ifng_sig$X,deg_ifng_sig$d,sep="_")
    m=as.numeric(nrow(deg_aging))
    n=as.numeric(nrow(deg_ifng))
    q=as.numeric(length((intersect(na.omit(deg_ifng_sig$X_d),na.omit(deg_aging_sig$X_d)))))
    k=as.numeric(nrow(deg_ifng_sig))
    print(paste(q,m,n,k,sep=","))
    print(paste(fn_out,"is:",sep=" "))
    print(phyper(q, m, n, k))
    pt<-phyper(q, m, n, k)
    require(VennDiagram)
    temp<-venn.diagram(list(deg_Ifng_sig_5=deg_ifng_sig$X_d,deg_aging_sig=deg_aging_sig$X_d),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = paste("Oral p=",pt,sep=" "),height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
    pdf(file=paste(fn_out,"redo.pdf",sep="."))
    grid.draw(temp)
    dev.off()
    deg_ifng_aging<-merge(deg_ifng_sig,deg_aging_sig,by.x="X_d",by.y="X_d",all=T,suffixes = c(".ifng",".aging"))
    write.csv(deg_ifng_aging,paste("./data_output/",fn_out,".redo",sep="")) 
}
venn_deg_p("Stem_treatment_Ifng.DEG.wilcox.csv","Stem_age.DEG.wilcox.csv","venn_stem_deg")
venn_deg_p("TA_treatment_Ifng.DEG.wilcox.csv","TA_age.DEG.wilcox.csv","venn_TA_deg")
venn_deg_p("Enterocyte-Prog_treatment_Ifng.DEG.wilcox.csv","Early.Enterocytes.Progenitors_age.DEG.wilcox.csv","venn_Enterocyte-Prog_deg")
venn_deg_p("Enterocyte_treatment_Ifng.DEG.wilcox.csv","Enterocyte_age.DEG.wilcox.csv","venn_Enterocyte_deg")
venn_deg_p("Enteroendocrine_treatment_Ifng.DEG.wilcox.csv","Enteroendocrine_age.DEG.wilcox.csv","venn_Enteroendocrine_deg")
venn_deg_p("Goblet-Paneth_treatment_Ifng.DEG.wilcox.csv","Goblet-Paneth_age.DEG.wilcox.unique.csv","venn_Goblet-Paneth_deg")
venn_deg_p("Tuft_treatment_Ifng.DEG.wilcox.csv","Tuft_age.DEG.wilcox.csv","venn_Tuft_deg")

## gsea

geneset<-read.csv("./data/geneset_Kim_Haber_etal.csv",head=T)
gsea<-function(count_geneset_cn3,geneset,celltype) {
    require(ggpubr)
    require(reshape2)
    cn<-c(paste(celltype,"sc_control",sep="."),paste(celltype,"sc_ifnr",sep="."))
    print(cn)
    count_geneset_cn3<-count_geneset_cn3[,cn]
    colnames(geneset)<-c("setName","gene","Direction")
    count_geneset_cn3<-count_geneset_cn3[geneset$gene,]
    csel<-count_geneset_cn3[apply(count_geneset_cn3 > 0, 1, sum)>ncol(count_geneset_cn3)*0.2,]
    csel<-t(scale(t(csel)))
    csel<-as.data.frame(csel)
    gsel<-geneset[geneset$gene %in% rownames(csel) & !is.na(geneset$gene),]
    gcsel<-merge(gsel,csel,by.x="gene",by.y="row.names",all.x = T)
    gcsel[,4:ncol(gcsel)]<- gcsel[,4:ncol(gcsel)]*gcsel$Direction
    gcsel_long<-melt(gcsel,
                     id.vars = "setName",
                     measure.vars = colnames(count_geneset_cn3),
                     variable.name='setNames',
                     value.name='exp_scaled')
    gcsel_long$group<-(str_split_fixed(gcsel_long$setNames,"_",2)[,2])
    gcsel_long$group<-factor(gcsel_long$group,levels = c("control","ifnr"))
    print(nrow(gcsel_long[gcsel_long$group=="control",]))
    print(nrow(gcsel_long[gcsel_long$group=="ifnr",]))
    pdf(paste(celltype,"pdf",sep="."),width=8, height=9)
    p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group", palette = "jco",facet.by = "setName", short.panel.labs = FALSE)
    q<-p + stat_compare_means(label = "p.format")+ggtitle(paste(celltype,"gene_number:",nrow(gcsel_long[gcsel_long$group=="control",]),sep=" "))
    print(q)
    dev.off()
}
cp<-unique(str_split_fixed(colnames(cluster.averages.RNA),"\\.",2)[,1])
for (i in cp){gsea(cluster.averages.RNA,geneset,i)}

geneset_markers10<-markers10
geneset_markers10$setName<-geneset_markers10$cluster
geneset_markers10$gene<-geneset_markers10$gene
geneset_markers10<-geneset_markers10[,c(3,1)]
geneset_markers10$Direction<-1
for (i in cp){gsea(cluster.averages.RNA,geneset_markers10,i)}

## heatmap of mahdi's markers10
markers10$cluster_gene<-paste(markers10$cluster,markers10$gene,sep="#")
immune.combined2.markers10<-subset(immune.combined2,features = as.character(markers10$gene))
DoHeatmap(immune.combined2.markers10, assay = "RNA", features = as.character(markers10$gene), size = 2)+theme(axis.text.y = element_text(size = 4))+ scale_y_discrete(limit=markers10$gene,labels=markers10$cluster_gene)




### heatmap
## Noticed do scale the data before heatmap, data <- ScaleData(object = data, features = rownames(data))
## https://github.com/satijalab/seurat/issues/1369

immune.combined2.markers10<-subset(immune.combined2,features = as.character(markers10$gene))
DoHeatmap(immune.combined2.markers10, features = as.character(markers10$gene), size = 3)
markers10$cluster_gene<-paste(markers10$cluster,markers10$gene,sep="_")
# adjust cell cluter rank
immune.combined2.markers10$celltype_new.treatment<-factor(immune.combined2.markers10$celltype_new.treatment,levels = c("Stem.sc_control","Stem.sc_ifnr","TA.sc_control","TA.sc_ifnr","Enterocyte-Prog.sc_control","Enterocyte-Prog.sc_ifnr","Enterocyte.sc_control","Enterocyte.sc_ifnr","Enteroendocrine.sc_control","Enteroendocrine.sc_ifnr","Goblet-Paneth.sc_control","Goblet-Paneth.sc_ifnr","Tuft.sc_control","Tuft.sc_ifnr"))
# adjust gene rank
ct<-c( "E_Stem","D_TA","F_Secretory_precursos","C_Early_Enterocytes_Prog","B_Late_Enterocytes_Prog","A_Enterocyte","I_Enteroendocrine","G_Goblet","H_Paneth","J_Tuft")
cts<-data.frame(ct=ct,orders=1:10)
markers10s<-markers10
for(i in 1:nrow(markers10s)){markers10s$num[i]<-cts[cts$ct==markers10s$cluster[i],"orders"]}
markers10s<-markers10s[order(markers10s$num),]
markers10s<-markers10s[order(markers10s$num,decreasing = T),]

p<-DoHeatmap(immune.combined2.markers10, assay = "RNA", features = as.character(markers10s[markers10s$gene!="Gm45123","gene"]), size = 2, group.by='celltype_new.treatment')+theme(axis.text.y = element_text(size = 4))+ scale_y_discrete(limit=markers10s[markers10s$gene!="Gm45123","gene"],labels=markers10s[markers10s$gene!="Gm45123","cluster_gene"])+ theme(
legend.title = element_text(size = 4),
legend.text = element_text(size = 4)
)
q<-p+ guides(shape = guide_legend(override.aes = list(size = 1)))
print(q)
## change heatmap color
DoHeatmap(immune.combined2.markers10, assay = "RNA", features = as.character(markers10s[markers10s$gene!="Gm45123","gene"]), size = 2, group.by='celltype_new.treatment')+theme(axis.text.y = element_text(size = 4))+ scale_y_discrete(limit=markers10s[markers10s$gene!="Gm45123","gene"],labels=markers10s[markers10s$gene!="Gm45123","cluster_gene"])+ theme(
legend.title = element_text(size = 4),
legend.text = element_text(size = 4)
) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))

## redo GSEA with all cell expression
### GSEA all cells
immune.combined2.matrix <- as.matrix(GetAssayData(immune.combined2, slot = "counts"))
colnames(immune.combined2.matrix)<-as.character(as.data.frame(immune.combined2$celltype_new.treatment)[,1])
immune.combined2.matrix<-as.data.frame(immune.combined2.matrix)

## with p1p2 markers
gsea2<-function(count_geneset_cn3,geneset,celltype) {
require(ggpubr)
require(reshape2)
cn<-c(paste(celltype,"sc_control",sep="."),paste(celltype,"sc_ifnr",sep="."))
print(cn)
count_geneset_cn3<-count_geneset_cn3[,colnames(count_geneset_cn3) %in% cn]
colnames(geneset)<-c("setName","gene","Direction")
count_geneset_cn3<-count_geneset_cn3[geneset$gene,]
csel<-count_geneset_cn3[apply(count_geneset_cn3 > 0, 1, sum)>ncol(count_geneset_cn3)*0.2,]
csel<-t(scale(t(csel)))
csel<-as.data.frame(csel)
gsel<-geneset[geneset$gene %in% rownames(csel) & !is.na(geneset$gene),]
gcsel<-merge(gsel,csel,by.x="gene",by.y="row.names",all.x = T)
gcsel[,4:ncol(gcsel)]<- gcsel[,4:ncol(gcsel)]*gcsel$Direction
gcsel_long<-melt(gcsel,
id.vars = "setName",
measure.vars = colnames(count_geneset_cn3),
variable.name='setNames',
value.name='exp_scaled')
gcsel_long$setNames<-as.character(gcsel_long$setNames)
gcsel_long$group<-(str_split_fixed(gcsel_long$setNames,"_",2)[,2])
gcsel_long$group<-(str_split_fixed(gcsel_long$group,"\\.",2)[,1])
gcsel_long$group<-factor(gcsel_long$group,levels = c("control","ifnr"))
print(nrow(gcsel_long[gcsel_long$group=="control",]))
print(nrow(gcsel_long[gcsel_long$group=="ifnr",]))
print(ncol(gcsel_long))
pdf(paste(celltype,"pdf",sep="."),width=3, height=28)
ylim1 = boxplot.stats(gcsel_long$exp_scaled)$stats[c(1, 5)]
p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group", outlier.shape=NA, palette = "jco",facet.by = "setName",ncol=1, short.panel.labs = FALSE)+scale_color_manual(values = c("gray","gray0"))
q<-p + stat_compare_means(label = "p.format",method="t.test",label.y=2)+ggtitle(paste(celltype,"gene_number:",nrow(gcsel),sep=" "))+  coord_cartesian(ylim = ylim1*1.05)
print(q)
dev.off()
}
for (i in cp){gsea2(immune.combined2.matrix,geneset,i)}
## with Dresden top10 markers
gsea2<-function(count_geneset_cn3,geneset,celltype) {
require(ggpubr)
require(reshape2)
cn<-c(paste(celltype,"sc_control",sep="."),paste(celltype,"sc_ifnr",sep="."))
print(cn)
count_geneset_cn3<-count_geneset_cn3[,colnames(count_geneset_cn3) %in% cn]
colnames(geneset)<-c("setName","gene","Direction")
count_geneset_cn3<-count_geneset_cn3[geneset$gene,]
csel<-count_geneset_cn3[apply(count_geneset_cn3 > 0, 1, sum)>ncol(count_geneset_cn3)*0.2,]
csel<-t(scale(t(csel)))
csel<-as.data.frame(csel)
gsel<-geneset[geneset$gene %in% rownames(csel) & !is.na(geneset$gene),]
gcsel<-merge(gsel,csel,by.x="gene",by.y="row.names",all.x = T)
gcsel[,4:ncol(gcsel)]<- gcsel[,4:ncol(gcsel)]*gcsel$Direction
gcsel_long<-melt(gcsel,
id.vars = "setName",
measure.vars = colnames(count_geneset_cn3),
variable.name='setNames',
value.name='exp_scaled')
gcsel_long$setNames<-as.character(gcsel_long$setNames)
gcsel_long$group<-(str_split_fixed(gcsel_long$setNames,"_",2)[,2])
gcsel_long$group<-(str_split_fixed(gcsel_long$group,"\\.",2)[,1])
gcsel_long$group<-factor(gcsel_long$group,levels = c("control","ifnr"))
print(nrow(gcsel_long[gcsel_long$group=="control",]))
print(nrow(gcsel_long[gcsel_long$group=="ifnr",]))
print(ncol(gcsel_long))
pdf(paste(celltype,"pdf",sep="."),width=3, height=28)
ylim1 = boxplot.stats(gcsel_long$exp_scaled)$stats[c(1, 5)]
p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group",outlier.shape=NA, palette = "jco",facet.by = "setName",ncol=1, short.panel.labs = FALSE)+scale_color_manual(values = c("gray","gray0"))
q<-p + stat_compare_means(label = "p.format",method="t.test",label.y=1)+ggtitle(paste(celltype,"gene_number:",nrow(gcsel),sep=" "))+  coord_cartesian(ylim = ylim1*1.05)
print(q)
dev.off()
}
for (i in cp){gsea2(immune.combined2.matrix,geneset_markers10,i)}


#### boxplot of cell-type markers expression in each cell type
markers_exp_boxplot<-function(count_geneset_cn3,geneset,celltype) {
require(ggpubr)
require(reshape2)
cn<-c(paste(celltype,"sc_control",sep="."),paste(celltype,"sc_ifnr",sep="."))
print(cn)
count_geneset_cn3<-count_geneset_cn3[,colnames(count_geneset_cn3) %in% cn]
colnames(geneset)<-c("setName","gene","Direction")
count_geneset_cn3<-count_geneset_cn3[geneset$gene,]
csel<-count_geneset_cn3[apply(count_geneset_cn3 > 0, 1, sum)>ncol(count_geneset_cn3)*0.2,]
csel<-t(scale(t(csel)))
csel<-as.data.frame(csel)
gsel<-geneset[geneset$gene %in% rownames(csel) & !is.na(geneset$gene),]
gcsel<-merge(gsel,csel,by.x="gene",by.y="row.names",all.x = T)
gcsel[,4:ncol(gcsel)]<- gcsel[,4:ncol(gcsel)]*gcsel$Direction
gcsel_long<-melt(gcsel,
id.vars = "setName",
measure.vars = colnames(count_geneset_cn3),
variable.name='setNames',
value.name='exp_scaled')
gcsel_long$setNames<-as.character(gcsel_long$setNames)
gcsel_long$group<-(str_split_fixed(gcsel_long$setNames,"_",2)[,2])
gcsel_long$group<-(str_split_fixed(gcsel_long$group,"\\.",2)[,1])
gcsel_long$group<-factor(gcsel_long$group,levels = c("control","ifnr"))
print(nrow(gcsel_long[gcsel_long$group=="control",]))
print(nrow(gcsel_long[gcsel_long$group=="ifnr",]))
print(ncol(gcsel_long))
pdf(paste(celltype,"pdf",sep="."),width=3, height=28)
ylim1 = boxplot.stats(gcsel_long$exp_scaled)$stats[c(1, 5)]
p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group",outlier.shape=NA, palette = "jco",facet.by = "setName",ncol=1, short.panel.labs = FALSE)+scale_color_manual(values = c("gray","gray0"))
q<-p + stat_compare_means(label = "p.format",method="t.test",label.y=1)+ggtitle(paste(celltype,"gene_number:",nrow(gcsel),sep=" "))+  coord_cartesian(ylim = ylim1*1.05)
print(q)
dev.off()
}


## boxplot of markers10
####
geneset_markers10$cp<-str_split_fixed(geneset_markers10$setName,"_",2)[,2]
geneset_markers10[geneset_markers10$cp=="Early_Enterocytes_Prog","cp"]<-"Enterocyte-Prog"
geneset_markers10[geneset_markers10$cp=="Goblet","cp"]<-"Goblet-Paneth"
geneset_markers10[geneset_markers10$cp=="Paneth","cp"]<-"Goblet-Paneth"
immune.combined2.matrix.markers10<-immune.combined2.matrix[geneset_markers10$gene,]
gsel<-geneset_markers10[geneset_markers10$gene %in% rownames(immune.combined2.matrix.markers10) & !is.na(geneset_markers10$gene),]
gcsel<-merge(gsel,immune.combined2.matrix.markers10,by.x="gene",by.y="row.names",all.x = T)
gcsel$cp_gene<-paste(gcsel$cp,gcsel$gene,sep="#")

gcsel_long<-melt(gcsel,
id.vars = "cp_gene",
measure.vars = colnames(gcsel)[5:ncol(gcsel)],
variable.name='setNames',
value.name='exp_scaled')
gcsel_long$celltype<-str_split_fixed(gcsel_long$setNames,"\\.",2)[,1]
gcsel_long$cp<-str_split_fixed(gcsel_long$cp_gene,"#",2)[,1]
gcsel_long$treatment<-str_split_fixed(gcsel_long$setNames,"sc_",2)[,2]
gcsel_long$treatment<-str_split_fixed(gcsel_long$treatment,"\\.",2)[,1]
gcsel_long2<-gcsel_long[gcsel_long$cp==gcsel_long$celltype,]
gcsel_long2$exp_scaled_log2<-log2(as.numeric(gcsel_long2$exp_scaled)+0.1)
p <- ggboxplot(gcsel_long2, x = "cp_gene", y = "exp_scaled_log2",
color = "treatment", palette = "jco",
facet.by = "celltype", scales = "free",short.panel.labs = FALSE,ncol=1)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p + stat_compare_means(aes(group = treatment),label = "p.signif")+scale_color_manual(values = c("gray","gray0"))+ggtitle("ns: p > 0.05; *: p <= 0.05; **: p <= 0.01; ***: p <= 0.001; ****: p <= 0.0001")+theme(plot.title = element_text(size = 8))
p + stat_compare_means(aes(group = treatment),label = "p.signif",method = "t.test")+scale_color_manual(values = c("gray","gray0"))+ggtitle("ns: p > 0.05; *: p <= 0.05; **: p <= 0.01; ***: p <= 0.001; ****: p <= 0.0001")+theme(plot.title = element_text(size = 8))

### check Kit+ cell ratio

VlnPlot(immune.combined2,c("Kit"),group.by="celltype_new")
immune.combined2$celltype_new.treatment<-factor(immune.combined2$celltype_new.treatment,levels = c("Stem.sc_control","Stem.sc_ifnr","TA.sc_control","TA.sc_ifnr","Enterocyte-Prog.sc_control","Enterocyte-Prog.sc_ifnr","Enterocyte.sc_control","Enterocyte.sc_ifnr","Enteroendocrine.sc_control","Enteroendocrine.sc_ifnr","Goblet-Paneth.sc_control","Goblet-Paneth.sc_ifnr","Tuft.sc_control","Tuft.sc_ifnr"))
VlnPlot(immune.combined2,c("Kit"),group.by="celltype_new.treatment")
Kit_pos<-subset(immune.combined2, subset = Kit > 0)

md_pos<-data.frame("Var1"=Kit_pos$celltype_new.treatment)
## count celltype_treatment number and calculate the ratio of which in total
data.ratio.pos <- dplyr::group_by(md_pos,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratio.pos$treatment<-str_split_fixed(data.ratio.pos$Var1,"_",2)[,2]
data.ratio.pos$celltype<-str_split_fixed(data.ratio.pos$Var1,"\\.",2)[,1]
## recalculate the ratio in treatment 
data.ratio2 <- data.ratio.pos %>% dplyr::group_by(treatment) %>% dplyr::mutate(ratio2=count/sum(count))
library(RColorBrewer)
## Kit.pos.pdf
data.ratio2$celltype<-factor(data.ratio2$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
ggplot(data.ratio2, aes(x = treatment, y = ratio2, fill = celltype, label = round(ratio2,3))) +
geom_bar(stat = "identity") +
geom_text(size = 3, position = position_stack(vjust = 0.5)) + scale_fill_brewer(palette="Set3")+ylab("Kit+ cell ratio")

## Kit+ ratio in each cell type
md_all<-data.frame("Var2"=immune.combined2$celltype_new.treatment)
data.ratio.all <- dplyr::group_by(md_all,Var2) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratios<-merge(data.ratio.pos,data.ratio.all,by.x="Var1",by.y = "Var2",all = T, suffixes = c(".kitPos",".all"))
data.ratios$Freq<-data.ratios$count.kitPos/data.ratios$count.all
data.ratios$celltype<-str_split_fixed(data.ratios$Var1,"\\.",2)[,1]
data.ratios$treatment<-str_split_fixed(data.ratios$Var1,"_",2)[,2]
data.ratios$celltype<-factor(data.ratios$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
ggplot(data=data.ratios, aes(x=celltype, y=Freq, fill=treatment)) +
geom_bar(stat="identity", position=position_dodge())+
geom_text(aes(label=round(Freq,3)), vjust=-0.3, color="black",
position = position_dodge(0.9), size=3)+
scale_fill_manual(values = c("gray","gray0"))

## draw cell-type composition in all cells
data.ratio.all$celltype<-str_split_fixed(data.ratio.all$Var2,"\\.",2)[,1]
data.ratio.all$treatment<-str_split_fixed(data.ratio.all$Var2,"_",2)[,2]
data.ratio.all$celltype<-factor(data.ratio.all$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
ggplot(data=data.ratio.all, aes(x=celltype, y=ratio, fill=treatment)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(ratio,3)), vjust=-0.3, color="black",
              position = position_dodge(0.9), size=3)+
    scale_fill_manual(values = c("gray","gray0"))

## version 2
# redo Dresden top10 markers
load("./data/Point5.RData")
library(Seurat)
markers10s_genes<-markers10s[(markers10s$gene!="Gm45123")&(markers10s$cluster!="F_Secretory_precursos"),]
DoHeatmap(immune.combined2.markers10, assay = "RNA", features = as.character(markers10s_genes$gene), size = 2, group.by='celltype_new.treatment')+theme(axis.text.y = element_text(size = 4))+ scale_y_discrete(limit=markers10s_genes$gene,labels=markers10s_genes$cluster_gene)+ theme(
legend.title = element_text(size = 8),
legend.text = element_text(size = 8)
) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
# do top10 markers
top10<-read.csv("./data/top10-conserved_markers.csv",header = T,stringsAsFactors = F)
library(stringr)
top10$cluster<-str_split_fixed(top10$cluster,"-",2)[,1]
top10$celltype_new<-celltype_new[(as.numeric(top10$cluster)+1)]
immune.combined<-readRDS("./data/immune.combined.rds")
immune.combined2 <- ScaleData(object = immune.combined, features = rownames(immune.combined))
rm(immune.combined)
top10$celltype_gene<-paste(top10$celltype_new,top10$gene,sep="_")
immune.combined2$celltype_new<-celltype_new[immune.combined2$seurat_clusters]
immune.combined2$celltype_new.treatment<-paste(immune.combined2$celltype_new,immune.combined2$treatment,sep=".")
immune.combined2.top10$celltype_new.treatment<-factor(immune.combined2.top10$celltype_new.treatment,levels = c("Stem.sc_control","Stem.sc_ifnr","TA.sc_control","TA.sc_ifnr","Enterocyte-Prog.sc_control","Enterocyte-Prog.sc_ifnr","Enterocyte.sc_control","Enterocyte.sc_ifnr","Enteroendocrine.sc_control","Enteroendocrine.sc_ifnr","Goblet-Paneth.sc_control","Goblet-Paneth.sc_ifnr","Tuft.sc_control","Tuft.sc_ifnr"))
immune.combined2.top10<-subset(immune.combined2,features = as.character(top10$gene))
immune.combined2.top10$celltype_new.treatment<-factor(immune.combined2.top10$celltype_new.treatment,levels = c("Stem.sc_control","Stem.sc_ifnr","TA.sc_control","TA.sc_ifnr","Enterocyte-Prog.sc_control","Enterocyte-Prog.sc_ifnr","Enterocyte.sc_control","Enterocyte.sc_ifnr","Enteroendocrine.sc_control","Enteroendocrine.sc_ifnr","Goblet-Paneth.sc_control","Goblet-Paneth.sc_ifnr","Tuft.sc_control","Tuft.sc_ifnr"))
top10s<-top10
top10s$num<-"num"
ct2<-c( "Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft")
cts2<-data.frame(ct2=ct2,orders=1:7)
for(i in 1:nrow(top10s)){top10s$num[i]<-cts2[cts2$ct2==top10s$celltype_new[i],"orders"]}
top10s<-top10s[order(top10s$num,decreasing = T),

#not needed but i did: immune.combined2.top10 <- ScaleData(object = immune.combined2.top10, features = as.character(top10s$gene))
DoHeatmap(immune.combined2.top10, assay = "RNA", features = rev(as.character(top10s$gene)), size = 2, group.by='celltype_new.treatment')+theme(axis.text.y = element_text(size = 4))+scale_y_discrete(labels=top10s$celltype_gene)+ theme(
legend.title = element_text(size = 8),
legend.text = element_text(size = 8)
)+ scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))

DoHeatmap(immune.combined2.top10, assay = "RNA", features = rev(as.character(top10s$gene)), size = 2, group.by='celltype_new.treatment')+theme(axis.text.y = element_text(size = 4))+ theme(
legend.title = element_text(size = 8),
legend.text = element_text(size = 8)
)+ scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))

## average expression version2
cluster.averages.RNA.markers10_v2<-cluster.averages.RNA.markers10[!grepl("Secretory",row.names(cluster.averages.RNA.markers10)),]
pheatmap(na.omit(cluster.averages.RNA.markers10_v2[apply(cluster.averages.RNA.markers10_v2,1,sd)!=0,]),cluster_rows = T,cluster_cols = F,fontsize_row= 3,scale="row")

### odds ratio
table(data.ratio$Var1)
library(ggpubr)
##
data.ratio <- my.data %>% dplyr::group_by(sampleID,celltype_new) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratio$hashtag<-str_split_fixed(data.ratio$sampleID,"_",2)[,2]
data.ratio$treatment="Control"
data.ratio[(data.ratio$hashtag=="shash5"|data.ratio$hashtag=="shash9"),]$treatment="Ifng"
ggbarplot(data.ratio, x = "celltype_new", y = "ratio", add = "mean_sd", color = "treatment", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = treatment),method="t.test",paired=F,label = "p.format",label.y = 0.4)+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))+scale_color_manual(values = c("grey","grey0"))+ggtitle("t.test_pairedF")+xlab("")+ylab("Cell ratio(%)")
data.ratio$celltype_new<-factor(data.ratio$celltype_new,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
ggbarplot(data.ratio, x = "celltype_new", y = "ratio", add = "mean_sd", color = "treatment", palette = "jco", position = position_dodge(0.8))+ stat_compare_means(aes(group = treatment),method="t.test",paired=F,label = "p.format",label.y = 0.4)+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))+scale_color_manual(values = c("grey","grey0"))+ggtitle("t.test_pairedF")+xlab("")+ylab("Cell ratio(%)")
###go
data.ratio2<-data.ratio %>% dplyr::group_by(celltype_new,treatment)  %>% dplyr::mutate(celltypes_sum = sum(count))
data.ratio2<-unique(data.ratio2[,c(2,6,7)])
data.ratio2<-data.ratio2 %>% dplyr::group_by(treatment)  %>% dplyr::mutate(celltypes_sum_ratio = celltypes_sum/sum(celltypes_sum))
library(reshape2)
data.ratio2_wide<-dcast(data.ratio2,celltype_new~data.ratio2$treatment,value.var = 'celltypes_sum_ratio')
data.ratio2_wide$odds_ratio<-(data.ratio2_wide$Ifng/(1-data.ratio2_wide$Ifng))/(data.ratio2_wide$Control/(1-data.ratio2_wide$Control))
data.ratio2<-data.ratio2 %>% dplyr::group_by(treatment)  %>% dplyr::mutate(treatment_sum = sum(celltypes_sum))
data.ratio2_wide2<-dcast(data.ratio2,celltype_new~data.ratio2$treatment,value.var = c("celltypes_sum"))
data.ratio2_wide2_2<-dcast(data.ratio2,celltype_new~data.ratio2$treatment,value.var = c("treatment_sum"))
data.ratio2_wide2<-merge(data.ratio2_wide2,data.ratio2_wide2_2,by="celltype_new",suffixes=c(".specific",".total"))
data.ratio2_wide2$k<-(data.ratio2_wide2$Control.specific+data.ratio2_wide2$Ifng.specific)
colnames(data.ratio2_wide2)<-c("celltype_new","x_Control.specific","q_Ifng.specific","n_Control.total","m_Ifng.total","k")
fisher_conf<-function(q,m,n,x){tmp<-fisher.test(matrix(c(q,(m-q),x,(n-x)),nrow=2));return(tmp$conf.int[c(1,2)])}
for (i in 1:nrow(data.ratio2_wide2)){data.ratio2_wide2$conf_1[i]<-fisher_conf(data.ratio2_wide2[i,"q_Ifng.specific"],data.ratio2_wide2[i,"m_Ifng.total"],data.ratio2_wide2[i,"n_Control.total"],data.ratio2_wide2[i,"x_Control.specific"])[1];data.ratio2_wide2$conf_2[i]<-fisher_conf(data.ratio2_wide2[i,"q_Ifng.specific"],data.ratio2_wide2[i,"m_Ifng.total"],data.ratio2_wide2[i,"n_Control.total"],data.ratio2_wide2[i,"x_Control.specific"])[2]}
data.ratio2_wide2_m<-merge(data.ratio2_wide,data.ratio2_wide2,by = "celltype_new",all = T)
library(ggplot2)
## SingleIfng_OddsRatio_v2.pdf
data.ratio2_wide2$hypergeometric_p_upper<-phyper(data.ratio2_wide2$q_Ifng.specific,data.ratio2_wide2$m_Ifng.total,data.ratio2_wide2$n_Control.total,data.ratio2_wide2$k,lower.tail = F)
data.ratio2_wide2_m$hypergeometric_p_upper<-data.ratio2_wide2$hypergeometric_p_upper
data.ratio2_wide2$hypergeometric_p_lower<-phyper(data.ratio2_wide2$q_Ifng.specific,data.ratio2_wide2$m_Ifng.total,data.ratio2_wide2$n_Control.total,data.ratio2_wide2$k)
data.ratio2_wide2_m$hypergeometric_p_lower<-data.ratio2_wide2$hypergeometric_p_lower
data.ratio2_wide2_m$hypergeometric_p_tails<-pmin(data.ratio2_wide2_m$hypergeometric_p_lower,data.ratio2_wide2_m$hypergeometric_p_upper)
ggplot(data.ratio2_wide2_m, aes(x=celltype_new, y=odds_ratio)) +
geom_bar(stat="identity", color="black") +
geom_errorbar(aes(ymin=conf_1, ymax=conf_2), width=.2,
position=position_dodge(.9)) + coord_flip()+scale_x_discrete(limits = rev(levels(factor(data.ratio2_wide2_m$celltype_new))))+theme_classic(base_size = 20) + geom_hline(yintercept=1, linetype="dashed", color = "black")+geom_text(aes(label=signif(hypergeometric_p_tails,digits = 3)),size = 3, position=position_dodge(width=1),hjust = -1,vjust=2)

## TA secretory bias
immune.combined2$celltype_new.treatment<-factor(immune.combined2$celltype_new.treatment,levels = c("Stem.sc_control","Stem.sc_ifnr","TA.sc_control","TA.sc_ifnr","Enterocyte-Prog.sc_control","Enterocyte-Prog.sc_ifnr","Enterocyte.sc_control","Enterocyte.sc_ifnr","Enteroendocrine.sc_control","Enteroendocrine.sc_ifnr","Goblet-Paneth.sc_control","Goblet-Paneth.sc_ifnr","Tuft.sc_control","Tuft.sc_ifnr"))
Idents(immune.combined2)<-immune.combined2$celltype_new.treatment
sec<-c("Muc2","Kit","Lyz1","Chga","Dclk1")
VlnPlot(immune.combined2,features=sec,idents = c("Stem.sc_control","Stem.sc_ifnr","TA.sc_control","TA.sc_ifnr"),ncol=2)
VlnPlot(immune.combined2,features=sec,idents = c("TA.sc_control","TA.sc_ifnr"),ncol=2)
gene_pos_composition<-function(gene){
    expr <- FetchData(object = immune.combined2, vars = gene)
    gene_pos<-immune.combined2[, which(x = expr > 0)]
    md_pos<-data.frame("Var1"=gene_pos$celltype_new.treatment)
    data.ratio.pos <- dplyr::group_by(md_pos,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
    library(stringr)
    data.ratio.pos$treatment<-str_split_fixed(data.ratio.pos$Var1,"_",2)[,2]
    data.ratio.pos$celltype<-str_split_fixed(data.ratio.pos$Var1,"\\.",2)[,1]
    data.ratio2 <- data.ratio.pos %>% dplyr::group_by(treatment) %>% dplyr::mutate(ratio2=count/sum(count))
    library(RColorBrewer)
    data.ratio2$celltype<-factor(data.ratio2$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
    ggplot(data.ratio2, aes(x = treatment, y = ratio2, fill = celltype, label = round(ratio2,3))) +
        geom_bar(stat = "identity") +
        geom_text(size = 3, position = position_stack(vjust = 0.5)) + scale_fill_brewer(palette="Set3")+ylab("gene+ cell ratio")+ylab(paste(gene,"-positive-cell composition",sep=""))
}
gene_pos_composition("Kit")
gene_pos_ratio<-function(gene){
    require(stringr)
    require(ggplot2)
    expr <- FetchData(object = immune.combined2, vars = gene)
    gene_pos<-immune.combined2[, which(x = expr > 0)]
    md_pos<-data.frame("Var1"=gene_pos$celltype_new.treatment)
    data.ratio.pos <- dplyr::group_by(md_pos,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
    data.ratio.pos$treatment<-str_split_fixed(data.ratio.pos$Var1,"_",2)[,2]
    data.ratio.pos$celltype<-str_split_fixed(data.ratio.pos$Var1,"\\.",2)[,1]
    md_all<-data.frame("Var2"=immune.combined2$celltype_new.treatment)
    data.ratio.all <- dplyr::group_by(md_all,Var2) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
    data.ratios<-merge(data.ratio.pos,data.ratio.all,by.x="Var1",by.y = "Var2",all = T, suffixes = c(".genePos",".all"))
    data.ratios$Freq<-data.ratios$count.genePos/data.ratios$count.all
    data.ratios$celltype<-str_split_fixed(data.ratios$Var1,"\\.",2)[,1]
    data.ratios$treatment<-str_split_fixed(data.ratios$Var1,"_",2)[,2]
    data.ratios$celltype<-factor(data.ratios$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
    ggplot(data=data.ratios, aes(x=celltype, y=Freq, fill=treatment)) +
        geom_bar(stat="identity", position=position_dodge())+
        geom_text(aes(label=round(Freq,3)), vjust=-0.3, color="black",
                  position = position_dodge(0.9), size=3)+
        scale_fill_manual(values = c("gray","gray0"))+ylab(paste(gene,"-positive-cell ratio",sep=""))+xlab("")
}
gene_pos_ratio("Kit")
for(i in sec){pdf(paste(i,"pos.pdf",sep="."));print(gene_pos_composition(i));print(gene_pos_ratio(i));dev.off()}
## for TA
gene_pos_ratio_TA<-function(gene){
require(stringr)
require(ggplot2)
expr <- FetchData(object = immune.combined2, vars = gene)
gene_pos<-immune.combined2[, which(x = expr > 0)]
md_pos<-data.frame("Var1"=gene_pos$celltype_new.treatment)
data.ratio.pos <- dplyr::group_by(md_pos,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratio.pos$treatment<-str_split_fixed(data.ratio.pos$Var1,"_",2)[,2]
data.ratio.pos$celltype<-str_split_fixed(data.ratio.pos$Var1,"\\.",2)[,1]
md_all<-data.frame("Var2"=immune.combined2$celltype_new.treatment)
data.ratio.all <- dplyr::group_by(md_all,Var2) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratios<-merge(data.ratio.pos,data.ratio.all,by.x="Var1",by.y = "Var2",all = T, suffixes = c(".genePos",".all"))
data.ratios$Freq<-data.ratios$count.genePos/data.ratios$count.all
data.ratios$celltype<-str_split_fixed(data.ratios$Var1,"\\.",2)[,1]
data.ratios$treatment<-str_split_fixed(data.ratios$Var1,"_",2)[,2]
data.ratios$celltype<-factor(data.ratios$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
ggplot(data=data.ratios[data.ratios$celltype=="TA",], aes(x=celltype, y=Freq, fill=treatment)) +
geom_bar(stat="identity", position=position_dodge())+
geom_text(aes(label=round(Freq,3)), vjust=-0.3, color="black",
position = position_dodge(0.9), size=3)+
scale_fill_manual(values = c("gray","gray0"))+ylab(paste(gene,"-positive-cell ratio",sep=""))+xlab("")
}
for(i in sec){pdf(paste(i,"pos.TA.pdf",sep="."));print(gene_pos_ratio_TA(i));dev.off()}
## in sampleID with SD
immune.combined2$sampleID<-str_split_fixed(row.names(as.data.frame(immune.combined2$orig.ident)),"\\.",2)[,2]
immune.combined2$celltype_new.treatment.sid<-paste(immune.combined2$celltype_new.treatment,immune.combined2$sampleID,sep="_")
# all
gene_pos_ratio_sd<-function(gene){
require(stringr)
require(ggplot2)
expr <- FetchData(object = immune.combined2, vars = gene)
gene_pos<-immune.combined2[, which(x = expr > 0)]
md_pos<-data.frame("Var1"=gene_pos$celltype_new.treatment.sid)
data.ratio.pos <- dplyr::group_by(md_pos,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratio.pos$treatment<-str_split_fixed(data.ratio.pos$Var1,"_",3)[,2]
data.ratio.pos$celltype<-str_split_fixed(data.ratio.pos$Var1,"\\.",2)[,1]
md_all<-data.frame("Var2"=immune.combined2$celltype_new.treatment.sid)
data.ratio.all <- dplyr::group_by(md_all,Var2) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratios<-merge(data.ratio.pos,data.ratio.all,by.x="Var1",by.y = "Var2",all = T, suffixes = c(".genePos",".all"))
data.ratios$Freq<-data.ratios$count.genePos/data.ratios$count.all
data.ratios$celltype<-str_split_fixed(data.ratios$Var1,"\\.",2)[,1]
data.ratios$treatment<-str_split_fixed(data.ratios$Var1,"_",3)[,2]
data.ratios$celltype<-factor(data.ratios$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
plot1<-ggbarplot(data.ratios, x = "treatment", y = "Freq", add = "mean_sd",
color = "treatment", palette = "jco", facet.by = "celltype", short.panel.labs = FALSE,scales="free_y",
position = position_dodge(0.8))+stat_compare_means(aes(group = treatment), method = "t.test",label = "p.format", paired=T,label.x.npc=c("middle"))+ggtitle("t.test.paired")+scale_color_manual(values = c("gray","gray0"))
plot2<-ggbarplot(data.ratios, x = "treatment", y = "Freq", add = "mean_sd",
color = "treatment", palette = "jco", facet.by = "celltype", short.panel.labs = FALSE,
position = position_dodge(0.8))+stat_compare_means(aes(group = treatment), method = "t.test",label = "p.format", paired=T,label.x.npc=c("middle"))+ggtitle("t.test.paired")+scale_color_manual(values = c("gray","gray0"))
print(plot1)
print(plot2)
}
for(i in sec){pdf(paste(i,"pos.sd.pdf",sep="."));gene_pos_ratio_sd(i);dev.off()}
# TA
gene_pos_ratio_sd_TA<-function(gene){
require(stringr)
require(ggplot2)
expr <- FetchData(object = immune.combined2, vars = gene)
gene_pos<-immune.combined2[, which(x = expr > 0)]
md_pos<-data.frame("Var1"=gene_pos$celltype_new.treatment.sid)
data.ratio.pos <- dplyr::group_by(md_pos,Var1) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratio.pos$treatment<-str_split_fixed(data.ratio.pos$Var1,"_",3)[,2]
data.ratio.pos$celltype<-str_split_fixed(data.ratio.pos$Var1,"\\.",2)[,1]
md_all<-data.frame("Var2"=immune.combined2$celltype_new.treatment.sid)
data.ratio.all <- dplyr::group_by(md_all,Var2) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(ratio=count/sum(count))
data.ratios<-merge(data.ratio.pos,data.ratio.all,by.x="Var1",by.y = "Var2",all = T, suffixes = c(".genePos",".all"))
data.ratios$Freq<-data.ratios$count.genePos/data.ratios$count.all
data.ratios$celltype<-str_split_fixed(data.ratios$Var1,"\\.",2)[,1]
data.ratios$treatment<-str_split_fixed(data.ratios$Var1,"_",3)[,2]
data.ratios$celltype<-factor(data.ratios$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
plot2<-ggbarplot(data.ratios[data.ratios$celltype=="TA",], x = "treatment", y = "Freq", add = "mean_sd",
color = "treatment", palette = "jco", facet.by = "celltype", short.panel.labs = FALSE,
position = position_dodge(0.8))+stat_compare_means(aes(group = treatment), method = "t.test",label = "p.format", paired=T,label.x.npc=c("middle"))+ggtitle(paste(gene,"_t.test.paired"))+scale_color_manual(values = c("gray","gray0")) + theme(plot.title = element_text(hjust = 0.1,size=10),legend.title = element_text(size = 0))+ylab(paste(gene,"_positive_cell ratio"))
print(plot2)
}
pdf("TA.sec.pdf",width = 2,height = 3)
for(i in sec){gene_pos_ratio_sd_TA(i)}
dev.off()

## geneset Kim_P1/Kim_P2 on heatmap
# modify setName for genes with multi-sources
library(data.table)
geneset_dt<-dcast(setDT(geneset), gene~rowid(gene, prefix="sth"), value.var="setName")
geneset_dt$setName<-paste(geneset_dt$sth1,geneset_dt$sth2,geneset_dt$sth3,sep="#")
geneset_dt$setName<-gsub("#NA","",geneset_dt$setName)
geneset_dt<-geneset_dt[-1,]

cluster.averages.RNA.geneset<-na.omit(cluster.averages.RNA[as.character(geneset_dt$gene),])
# row.names(cluster.averages.RNA.geneset)<-paste(row.names(cluster.averages.RNA.geneset),geneset_dt[match(row.names(cluster.averages.RNA.geneset),geneset_dt$gene),"setName"],sep=".")
tmp<-geneset_dt[match(row.names(cluster.averages.RNA.geneset),geneset_dt$gene),"setName"]
row.names(cluster.averages.RNA.geneset)<-paste(row.names(cluster.averages.RNA.geneset),tmp$setName,sep=".")
cluster.averages.RNA.geneset<-cluster.averages.RNA.geneset[,ll2]
cluster.averages.RNA.geneset.Kim<-cluster.averages.RNA.geneset[grepl("Kim",row.names(cluster.averages.RNA.geneset)),]
cluster.averages.RNA.geneset.Kim.P<-cluster.averages.RNA.geneset.Kim[grepl("Kim_P",row.names(cluster.averages.RNA.geneset.Kim)),]
cluster.averages.RNA.geneset.Kim.others<-cluster.averages.RNA.geneset.Kim[!grepl("Kim_P",row.names(cluster.averages.RNA.geneset.Kim)),]
pheatmap(na.omit(cluster.averages.RNA.geneset.Kim.P[apply(cluster.averages.RNA.geneset.Kim.P,1,sd)!=0,]),cluster_rows = T,cluster_cols = F,fontsize_row= 3,scale="row")
pheatmap(na.omit(cluster.averages.RNA.geneset.Kim.others[apply(cluster.averages.RNA.geneset.Kim.others,1,sd)!=0,]),cluster_rows = T,cluster_cols = F,fontsize_row= 3,scale="row")
pheatmap(na.omit(cluster.averages.RNA.geneset[apply(cluster.averages.RNA.geneset,1,sd)!=0,]),cluster_rows = T,cluster_cols = F,fontsize_row= 3,scale="row")

## Vln for markers of Enteroendocrine/Goble-Paneth/Tuft in this data set

immune.combined2$celltype_new<-factor(immune.combined2$celltype_new,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))

pdf("SingleIfng_Tuft_Vln.pdf")
for(i in top10s[top10s$celltype_new=="Tuft","gene"]){plot<-VlnPlot(immune.combined2, features = i, split.by = "treatment", group.by = "celltype_new", pt.size = 0) + ggplot2::scale_fill_manual(values = c('grey', 'grey0'))+ggtitle(paste("Tuft_",i));print(plot)}
dev.off()
pdf("SingleIfng_Enteroendocrine_Vln.pdf")
for(i in top10s[top10s$celltype_new=="Enteroendocrine","gene"]){plot<-VlnPlot(immune.combined2, features = i, split.by = "treatment", group.by = "celltype_new", pt.size = 0) + ggplot2::scale_fill_manual(values = c('grey', 'grey0'))+ggtitle(paste("Enteroendocrine_",i));print(plot)}
dev.off()
pdf("SingleIfng_Goblet-Paneth_Vln.pdf")
for(i in top10s[top10s$celltype_new=="Goblet-Paneth","gene"]){plot<-VlnPlot(immune.combined2, features = i, split.by = "treatment", group.by = "celltype_new", pt.size = 0) + ggplot2::scale_fill_manual(values = c('grey', 'grey0'))+ggtitle(paste("Goblet-Paneth_",i));print(plot)}
dev.off()

## top10 markers in this data set on heatmap
cluster.averages.RNA.top10s<-na.omit(cluster.averages.RNA[as.character(top10s$gene),])
tmps<-top10s[match(row.names(cluster.averages.RNA.top10s),top10s$gene),"celltype_gene"]
tmps[49]<-"Enterocyte-Prog_Sis"
row.names(cluster.averages.RNA.top10s)<-tmps
cluster.averages.RNA.top10s<-cluster.averages.RNA.top10s[,ll2]
pheatmap(na.omit(cluster.averages.RNA.top10s[apply(cluster.averages.RNA.top10s,1,sd)!=0,]),cluster_rows = T,cluster_cols = F,fontsize_row= 3,scale="row")
aaa<-cluster.averages.RNA.top10s[rev(row.names(cluster.averages.RNA.top10s)),]
pheatmap(na.omit(aaa[apply(aaa,1,sd)!=0,]),cluster_rows = F,cluster_cols = F,fontsize_row= 3,scale="row")

##### check apoptotic_genes
load("./data/Point5.RData")
library(Seurat)
immune.combined2<-readRDS("./data/immune.combined2.rds")
# Check_apoptotic_genes.pdf
VlnPlot(immune.combined2, features = c("Bcl2","Bcl2l1","Bak1"), split.by = "treatment", group.by = "celltype_new", pt.size = 0)
# Check_PositiveRatio_apoptotic_genes_bcl2.pdf Check_PositiveRatio_apoptotic_genes_bcl2l1.pdf Check_PositiveRatio_apoptotic_genes_bak1.pdf
library(ggpubr)
gene_pos_ratio_sd("Bcl2")
gene_pos_ratio_sd("Bcl2l1")
gene_pos_ratio_sd("Bak1")
gene_pos_ratio("Bak1")
gene_pos_ratio("Bcl2")
# pseudo-bulk_apoptotic_genes
immune.combined2.RawCounts<-as.data.frame(immune.combined2[["RNA"]]@counts)
colnames(immune.combined2.RawCounts)<-immune.combined2$celltype_new.treatment.sid
library(stringr)
rm(immune.combined2)
# sum by sample_celltypes
immune.combined2.RawCounts<-t(rowsum(t(immune.combined2.RawCounts), colnames(immune.combined2.RawCounts)))
# sum by sample
immune.combined2.RawCounts2<-t(rowsum(t(immune.combined2.RawCounts), str_split_fixed(colnames(immune.combined2.RawCounts),"sc_",2)[,2]))
# normalize raw counts by total counts of sample, and scale it with median total counts os 8 samples
immune.combined2.RawCounts3<-sapply(colnames(immune.combined2.RawCounts),function(x){(immune.combined2.RawCounts[,x]/colSums(immune.combined2.RawCounts2)[str_split_fixed(x,"sc_",2)[,2]])*median(colSums(immune.combined2.RawCounts2))})
immune.combined2.RawCounts4<-t(rowsum(t(immune.combined2.RawCounts3), str_split_fixed(colnames(immune.combined2.RawCounts3),"sc_",2)[,2]))
immune.combined2.RawCounts4<-as.data.frame(immune.combined2.RawCounts4)
immune.combined2.RawCounts3<-as.data.frame(immune.combined2.RawCounts3)

# function for treatment with sd
cytokine_bar<-function(immune.combined2.RawCounts4,gene){
test<-as.data.frame(t(immune.combined2.RawCounts4[gene,]))
test$treatment<-str_split_fixed(row.names(test),"_",2)[,1]
test$treatment<-factor(test$treatment,levels = c("control","ifnr"))
p<-ggbarplot(test, x = "treatment", y = gene, add = "mean_sd",color = "treatment", palette = c("grey","grey0"))+ stat_compare_means(method = "t.test")+ylab(paste(gene," normalized counts (by sample),sd"))
print(p)
}
# function for celltype-treatment with sd
cytokine_bar_2<-function(immune.combined2.RawCounts3,gene){
test<-as.data.frame(t(immune.combined2.RawCounts3[gene,]))
test$treatment<-str_split_fixed(row.names(test),"_",4)[,2]
test$celltype<-str_split_fixed(row.names(test),"\\.",2)[,1]
test$celltype<-factor(test$celltype,levels = c("Stem","TA","Enterocyte-Prog","Enterocyte","Enteroendocrine","Goblet-Paneth","Tuft"))
test$treatment<-factor(test$treatment,levels = c("control","ifnr"))
p<-ggbarplot(test, x = "celltype", y = gene, add = "mean_sd",color = "treatment", palette = c("grey","grey0"),position = position_dodge(0.8))+ stat_compare_means(aes(group = treatment),method = "t.test",size=2)+ylab(paste(gene," normalized counts (by sample),sd"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
}

pdf("pseudo-bulk_apoptotic_genes.pdf")
cytokine_bar_2(immune.combined2.RawCounts3,"Bcl2")
cytokine_bar_2(immune.combined2.RawCounts3,"Bcl2l1")
cytokine_bar_2(immune.combined2.RawCounts3,"Bak1")
dev.off()
pdf("pseudo-bulk_apoptotic_genes_2.pdf")
cytokine_bar(immune.combined2.RawCounts4,"Bcl2")
cytokine_bar(immune.combined2.RawCounts4,"Bcl2l1")
cytokine_bar(immune.combined2.RawCounts4,"Bak1")
dev.off()

# Check_Mki67_PositiveRatio.pdf
gene_pos_ratio_sd("Mki67")
# Check_Mki67.pdf
VlnPlot(immune.combined2, features = c("Mki67"), split.by = "treatment", group.by = "celltype_new", pt.size = 0)
# pseudo-bulk_Mki67.pdf
pdf("pseudo-bulk_Mki67.pdf")
cytokine_bar_2(immune.combined2.RawCounts3,"Mki67")
cytokine_bar(immune.combined2.RawCounts4,"Mki67")
dev.off()

> saveRDS(immune.combined2.RawCounts3,"./data/immune.combined2.RawCounts3.RDS")
> saveRDS(immune.combined2.RawCounts4,"./data/immune.combined2.RawCounts4.RDS")

> FeaturePlot(immune.combined, features = c("Lgr5"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Olfm4"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Mki67"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Zg16"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Gm7849"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Lyz1"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Defa3"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Ang4"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Lyz1"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Agr2"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Fcgbp"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Gale"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Clca1"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Bak1"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Bcl2"),split.by ="treatment")
> FeaturePlot(immune.combined, features = c("Bcl2l1"),split.by ="treatment")


======================================================================
======================================================================
## Modify1
load("./data/Point5.RData")
library(pheatmap)
## heatmap of mhc2 genes from mahdi's dataset and order: SingleIfng_mhc2_mahdi.pdf
mhc2_mahdi<-c("H2.M3","H2.T23","Tap2","Tapbp","Cd74","H2.Aa","H2.Eb1","H2.Ab1","Ciita","H2.DMb1","H2.DMa","H2.Q7","Psmb8","B2m","Tap1","Calr","Pdia3","Psmb6")
cluster.averages.RNA.mhc2_mahdi<-(cluster.averages.RNA[as.character(mhc2_mahdi),])
cluster.averages.RNA.mhc2_mahdi<-cluster.averages.RNA.mhc2_mahdi[,ll2]
pheatmap(na.omit(cluster.averages.RNA.mhc2_mahdi[apply(cluster.averages.RNA.mhc2_mahdi,1,sd)!=0,]),cluster_rows = F,cluster_cols = F,fontsize_row= 6,scale="row",border_color = NA)

## Ifngr1 Ifngr2: RnaAnna_Ifngr1.pdf RnaAnna_Ifngr2.pdf
## wechat screenshot could pick up corlor RGB number
## use the same color and order: 
library(Seurat)
immune.combined2<-readRDS("./data/immune.combined2.rds")
mycolor<-c("#8a7db9","#8a7db9","#f16962","#f16962","#bd9129","#bd9129","#1db086","#1db086","#1aace6","#1aace6","#44a93e","#44a93e","#cd609e","#cd609e")
VlnPlot(immune.combined2,pt.size = 0, features = c("Ifngr1"))+ scale_fill_manual(values = mycolor)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
VlnPlot(immune.combined2,pt.size = 0, features = c("Ifngr2"))+ scale_fill_manual(values = mycolor)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

## umap plot for gene expression

# SingleIfng_UMAP_Markers.pdf
FeaturePlot(immune.combined, features = c("Lgr5", "Olfm4","Apoa1", "Muc2", "Defa26","Aplp1"),pt.size = 0.02,coord.fixed = 1,ncol = 2,cols = c("grey","red"))
# SingleIfng_UMAP_Markers_scaleMax3.pdf
FeaturePlot(immune.combined, features = c("Lgr5", "Olfm4","Apoa1", "Muc2", "Defa26","Aplp1"),pt.size = 0.02,coord.fixed = 1,ncol = 2,cols = c("grey","red"),max.cutoff = 3)

# SingleIfng_UMAP_MHC2_scaleMax3.pdf
FeaturePlot(immune.combined, features = c("H2.Aa","H2.Eb1"),pt.size = 0.02,split.by = "treatment",coord.fixed = 1,ncol = 2,cols = c("grey","red"),max.cutoff = 3)


===================================================
# Attached
# For further check
===================================================

## results were in:
## /Users/jlu/Desktop/Pro_Mine/Paper_Omid/Point5/figure/gsea_p1p2_allcells_ttest/
load("./data/Point5.RData")
library(stringr)
# merging on data.table type were not working very well, by.y="row.names" could not work, so tranfrome gsel to data.frame

gsea2(immune.combined2.matrix,geneset,"TA")
gsea2<-function(count_geneset_cn3,geneset,celltype) {
require(ggpubr)
require(reshape2)
cn<-c(paste(celltype,"sc_control",sep="."),paste(celltype,"sc_ifnr",sep="."))
print(cn)
count_geneset_cn3<-count_geneset_cn3[,colnames(count_geneset_cn3) %in% cn]
colnames(geneset)<-c("setName","gene","Direction")
count_geneset_cn3<-count_geneset_cn3[geneset$gene,]
csel<-count_geneset_cn3[apply(count_geneset_cn3 > 0, 1, sum)>ncol(count_geneset_cn3)*0.2,]
csel<-t(scale(t(csel)))
csel<-as.data.frame(csel)
gsel<-geneset[geneset$gene %in% rownames(csel) & !is.na(geneset$gene),]
gsel<-as.data.frame(gsel)
gcsel<-merge(gsel,csel,by.x="gene",by.y="row.names",all.x = T)
gcsel[,4:ncol(gcsel)]<- gcsel[,4:ncol(gcsel)]*gcsel$Direction
gcsel_long<-melt(gcsel,
id.vars = "setName",
measure.vars = colnames(count_geneset_cn3),
variable.name='setNames',
value.name='exp_scaled')
gcsel_long$setNames<-as.character(gcsel_long$setNames)
gcsel_long$group<-(str_split_fixed(gcsel_long$setNames,"_",2)[,2])
gcsel_long$group<-(str_split_fixed(gcsel_long$group,"\\.",2)[,1])
gcsel_long$group<-factor(gcsel_long$group,levels = c("control","ifnr"))
print(nrow(gcsel_long[gcsel_long$group=="control",]))
print(nrow(gcsel_long[gcsel_long$group=="ifnr",]))
print(ncol(gcsel_long))
pdf(paste(celltype,"pdf",sep="."),width=3, height=28)
ylim1 = boxplot.stats(gcsel_long$exp_scaled)$stats[c(1, 5)]
p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group", outlier.shape=NA, palette = "jco",facet.by = "setName",ncol=1, short.panel.labs = FALSE)+scale_color_manual(values = c("gray","gray0"))
q<-p + stat_compare_means(label = "p.format",method="t.test",label.y=2)+ggtitle(paste(celltype,"gene_number:",nrow(gcsel),sep=" "))+  coord_cartesian(ylim = ylim1*1.05)
print(q)
dev.off()
}
colnames(immune.combined2.matrix)[colnames(immune.combined2.matrix)=="Stem.sc_ifnr"]<-"StemTA.sc_ifnr"
colnames(immune.combined2.matrix)[colnames(immune.combined2.matrix)=="Stem.sc_control"]<-"StemTA.sc_control"
colnames(immune.combined2.matrix)[colnames(immune.combined2.matrix)=="TA.sc_control"]<-"StemTA.sc_control"
colnames(immune.combined2.matrix)[colnames(immune.combined2.matrix)=="TA.sc_ifnr"]<-"StemTA.sc_ifnr"
## StemTA.pdf
gsea2(immune.combined2.matrix,geneset,"StemTA")
library(Seurat)
immune.combined2<-readRDS("./data/immune.combined2.rds")
Lgr5_pos<-subset(immune.combined2, subset = Lgr5 > 0)
Olfm4_pos<-subset(immune.combined2, subset = Olfm4 > 0)

gsea3<-function(count_geneset_cn3,geneset,celltype) {
require(ggpubr)
require(reshape2)
cn<-c(paste(celltype,"sc_control",sep=""),paste(celltype,"sc_ifnr",sep=""))
print(cn)
count_geneset_cn3<-count_geneset_cn3[,colnames(count_geneset_cn3) %in% cn]
colnames(geneset)<-c("setName","gene","Direction")
count_geneset_cn3<-count_geneset_cn3[geneset$gene,]
csel<-count_geneset_cn3[apply(count_geneset_cn3 > 0, 1, sum)>ncol(count_geneset_cn3)*0.2,]
csel<-t(scale(t(csel)))
csel<-as.data.frame(csel)
gsel<-geneset[geneset$gene %in% rownames(csel) & !is.na(geneset$gene),]
gsel<-as.data.frame(gsel)
gcsel<-merge(gsel,csel,by.x="gene",by.y="row.names",all.x = T)
gcsel[,4:ncol(gcsel)]<- gcsel[,4:ncol(gcsel)]*gcsel$Direction
gcsel_long<-melt(gcsel,
id.vars = "setName",
measure.vars = colnames(count_geneset_cn3),
variable.name='setNames',
value.name='exp_scaled')
gcsel_long$setNames<-as.character(gcsel_long$setNames)
gcsel_long$group<-(str_split_fixed(gcsel_long$setNames,"_",2)[,2])
gcsel_long$group<-(str_split_fixed(gcsel_long$group,"\\.",2)[,1])
gcsel_long$group<-factor(gcsel_long$group,levels = c("control","ifnr"))
print(nrow(gcsel_long[gcsel_long$group=="control",]))
print(nrow(gcsel_long[gcsel_long$group=="ifnr",]))
print(ncol(gcsel_long))
pdf(paste("gene_positive","pdf",sep="."),width=3, height=28)
ylim1 = boxplot.stats(gcsel_long$exp_scaled)$stats[c(1, 5)]
p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group", outlier.shape=NA, palette = "jco",facet.by = "setName",ncol=1, short.panel.labs = FALSE)+scale_color_manual(values = c("gray","gray0"))
q<-p + stat_compare_means(label = "p.format",method="t.test",label.y=2)+ggtitle(paste(celltype,"gene_number:",nrow(gcsel),sep=" "))+  coord_cartesian(ylim = ylim1*1.05)
print(q)
dev.off()
}
# Lgr5_positive.pdf
immune.combined2.lgr5 <- as.matrix(GetAssayData(Lgr5_pos, slot = "counts"))
colnames(immune.combined2.lgr5)<-as.character(as.data.frame(Lgr5_pos$treatment)[,1])
immune.combined2.lgr5<-as.data.frame(immune.combined2.lgr5)
gsea3(immune.combined2.lgr5,geneset,"")
# Vln_Lgr5.pdf
VlnPlot(immune.combined2, features = "Lgr5", split.by = "treatment", group.by = "celltype_new", pt.size = 0, combine = FALSE)
# Vln_Olfm4.pdf
VlnPlot(immune.combined2, features = "Olfm4", split.by = "treatment", group.by = "celltype_new", pt.size = 0, combine = FALSE)
Olfm4_pos_2<-subset(immune.combined2, subset = Olfm4 > 2)
# Lgr5_pos_0.2<-subset(immune.combined2, subset = Lgr5 > 0.2) # could not work and cutoff=0.2 has no big difference with cutoff=0
rm(Lgr5_pos_0.2)
rm(immune.combined2)
# Olfm4_positive.pdf
immune.combined2.olfm4 <- as.matrix(GetAssayData(Olfm4_pos, slot = "counts"))
colnames(immune.combined2.olfm4)<-as.character(as.data.frame(Olfm4_pos$treatment)[,1])
immune.combined2.olfm4<-as.data.frame(immune.combined2.olfm4)
gsea3(immune.combined2.olfm4,geneset,"")
# Olfm4_cutoff2.pdf
immune.combined2.olfm4_2 <- as.matrix(GetAssayData(Olfm4_pos_2, slot = "counts"))
colnames(immune.combined2.olfm4_2)<-as.character(as.data.frame(Olfm4_pos_2$treatment)[,1])
immune.combined2.olfm4_2<-as.data.frame(immune.combined2.olfm4_2)
gsea3(immune.combined2.olfm4_2,geneset,"")
