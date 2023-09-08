## Author: Jing

setwd("./")
## pca plot on rnaseq of ifng treatment and control
count<-read.csv("./data/Counts_norm_DESeq2_AllGenes.csv",header = T,row.names = 1,check.names = F)
countPCA <- count[apply(count > 10, 1, sum)>2.99 & apply(log2(count+0.1), 1, IQR) >= 1.5, ]
t.log<-log(countPCA)
df<-t.log[ , ! apply( t.log , 2 , function(x) all(is.na(x)) ) ]
df <- df[!is.infinite(rowSums(df)),]
pca<-prcomp(t(df),center = TRUE,scale. = TRUE)
ss<-read.csv("./data/SampleSheet2.csv")
p1<-cbind(as.data.frame(pca$x[,1:2]),Group=ss$Group)
p2<-p1
p1<-p1[p1$Group!="Nrn1",]
ggplot(p1[p1$Group!="Nrn1",],aes(x=p1[,1],y=p1[,2],col=Group))+geom_point()+xlab("PC1")+ylab("PC2")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))+coord_fixed(ratio=1)+stat_ellipse(level = 0.8)
## add variance explanation percentage 
aa<-summary(pca)
ggplot(p1[p1$Group!="Nrn1",],aes(x=p1[,1],y=p1[,2],col=Group))+geom_point()+xlab(paste("PC1",aa$importance[2,"PC1"],sep = "_"))+ylab(paste("PC2",aa$importance[2,"PC2"],sep="_"))+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))+coord_fixed(ratio=1)

## volcano plot abs(log2fc)>0.5
deg<-read.csv("./data/Deseq2_AllGenes_ConVSIfng_omid.csv",header = T,row.names = 1,check.names = F)
mhc2<-read.csv("./data/KEGG_mhc.genes",header = T)
deg2<-deg
deg<-deg[!is.na(deg),]
library(RColorBrewer)
logFC<-deg2$log2FoldChange
padj<-deg2$padj
data <- data.frame(logFC=logFC,padj=padj)
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < 0.5)& data$logFC > -0.5] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 0.5] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -0.5] <- "down"
x_lim <- max(logFC,-logFC)
p<-ggplot(data[!is.na(data$sig),],aes(logFC,-1*log10(padj),color = sig))+geom_point(size=0.2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-0.5,0.5),linetype=4)
## label top20 genes
top20_ifng<-head(row.names(data[order(data$logFC),]),n=20)
top20_ifng<-c(top20_ifng,head(row.names(data[order(-data$logFC),]),n=20))
data$sign<-NA
row.names(data)<-row.names(deg2)
data[row.names(data) %in% top20_ifng,"sign"]=row.names(data[row.names(data) %in% top20_ifng,])
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < 5)& data$logFC > -5] <- "no"
data$sig[data$padj <= 0.05 & data$logFC >= 5] <- "up"
data$sig[data$padj <= 0.05 & data$logFC <= -5] <- "down"
p<-ggplot(data[!is.na(data$sig),],aes(logFC,-1*log10(padj),color = sig))+geom_point(size=0.2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-5,5),linetype=4)
p+geom_text_repel(aes(label = sign), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
####

library(ggrepel)
row.names(data)<-row.names(deg2)
data$sign<-NA
data[row.names(data) %in% mhc2$genename,"sign"]=row.names(data[row.names(data) %in% mhc2$genename,])
p<-ggplot(data[!is.na(data$sig),],aes(logFC,-1*log10(padj),color = sig))+geom_point(size=0.2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-0.5,0.5),linetype=4)
p+geom_text_repel(aes(label = sign), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)

##heatmap
diff<-deg[(!is.na(deg$log2FoldChange))&(!is.na(deg$padj)),]
diff_up<-diff[(diff$log2FoldChange>5)&(diff$padj<0.05),]
dim(diff_up)
diff<-deg[(!is.na(deg$log2FoldChange))&(!is.na(deg$padj)),]
diff_down<-diff[(diff$log2FoldChange<(-5))&(diff$padj<0.05),]
dim(diff_down)
# remove Nrn1
count2<-count[,c(1:8)]
diff_up_order<-diff_up[order(-diff_up$log2FoldChange),]
diff_down_order<-diff_down[order(diff_down$log2FoldChange),]
library(pheatmap)

# mhc2 heatmap
# mhc2 gene name should be removed from factor to match the target dataframe
pheatmap(log2(count2[mhc2$genename,]+0.1),cluster_rows = T)
head(row.names(count2))
head(row.names(count2[mhc2$genename,]))
head(mhc2$genename)
head(row.names(count2[as.character(mhc2$genename),]))
# count should be log2 and remove na value
count3<-na.omit(log2(count2[as.character(mhc2$genename),]+0.1))
# before heatmap, hcluster require the sd is not Zero, so remove rows with sd=0, and scale row to show heatmap better
pheatmap(count3[apply(count3,1,sd)!=0,],cluster_rows = T,cluster_cols = F,fontsize_row= 5,scale="row")

heatmap_log2exp<-function(genes){
count3<-na.omit(log2(count2[as.character(genes),]+0.1))
pheatmap(count3[apply(count3,1,sd)!=0,],cluster_rows = T,cluster_cols = F,fontsize_row= 2,scale="row")
}
heatmap_log2exp(row.names(diff_up_order))
heatmap_log2exp(row.names(diff_down_order))


## redo mhc2 expression in anna data
Counts_Norm<-read.csv("./data/Counts_norm_DESeq2_AllGenes.csv",head=T,row.names = 1)
Counts_Norm_MHCII<-Counts_Norm[row.names(Counts_Norm)%in%mhc2$genename,]
cn<-c("young_F_crypts_oral_rep1",
"young_F_crypts_oral_rep2",
"young_F_crypts_oral_rep3",
"young_F_crypts_oral_rep4",
"young_M_crypts_oral_rep1",
"young_M_crypts_oral_rep2",
"young_M_crypts_oral_rep3",
"young_M_crypts_oral_rep4",
"geriatric_F_crypts_oral_rep1",
"geriatric_F_crypts_oral_rep2",
"geriatric_F_crypts_oral_rep3",
"geriatric_F_crypts_oral_rep4",
"old_M_crypts_oral_rep1",
"old_M_crypts_oral_rep2",
"old_M_crypts_oral_rep3",
"old_M_crypts_oral_rep4",
"young_F_organoids_oral_rep1",
"young_F_organoids_oral_rep2",
"young_F_organoids_oral_rep3",
"young_F_organoids_oral_rep4",
"young_M_organoids_oral_rep1",
"young_M_organoids_oral_rep2",
"young_M_organoids_oral_rep3",
"young_M_organoids_oral_rep4",
"geriatric_F_organoids_oral_rep1",
"geriatric_F_organoids_oral_rep2",
"geriatric_F_organoids_oral_rep3",
"geriatric_F_organoids_oral_rep4",
"old_M_organoids_oral_rep1",
"old_M_organoids_oral_rep2",
"old_M_organoids_oral_rep3",
"old_M_organoids_oral_rep4")
Counts_Norm_MHCII_cn<-Counts_Norm_MHCII[,cn]
library(pheatmap)
library(stringr)
colnames(Counts_Norm_MHCII_cn)<-gsub("geriatric","old",colnames(Counts_Norm_MHCII_cn))
anno_col<-data.frame(str_split_fixed(colnames(Counts_Norm_MHCII_cn),"_",5)[,2],str_split_fixed(colnames(Counts_Norm_MHCII_cn),"_",5)[,1],str_split_fixed(colnames(Counts_Norm_MHCII_cn),"_",5)[,3],row.names=colnames(Counts_Norm_MHCII_cn))
colnames(anno_col) <- c("Gender","Age","Invivo_Invitro")
xxx<-log2(Counts_Norm_MHCII_cn+0.1)
pheatmap(xxx[apply(xxx,1,sd)!=0,],cluster_cols = F,fontsize_row=4,legend_labels=c("log"),annotation_col = anno_col,border_color = F,fontsize = 6,main="MHC2 genes expression in bulk RNAseq, log2(exp+1)",scale = "row")

>-log10(0.05)
1.3

## change cutoff for logfoldchange
data$sig2[(data$padj > 0.05|data$padj=="NA")|(data$logFC < 5)& data$logFC > -5] <- "no"
data$sig2[data$padj <= 0.05 & data$logFC >= 5] <- "up"
data$sig2[data$padj <= 0.05 & data$logFC <= -5] <- "down"
ggplot(data[!is.na(data$sig2),],aes(logFC,-1*log10(padj),color = sig2))+geom_point(size=0.2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-5,5),linetype=4)+geom_text_repel(aes(label = sign), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)
##ipa
#ipa<-read.table("./data/ipa_deg_ifng_logfc5_sig.txt",header = T,sep="\t",row.names = NULL)
ipa<-read.table("./data/ipa_deg_ifng_logfc5_sig.txt",header = T,sep="\t")
ipa_order<-ipa[order(-ipa$X.log.p.value.),]
ipa_order$Ingenuity.Canonical.Pathways<-factor(ipa_order$Ingenuity.Canonical.Pathways,levels = rev(levels(ipa_order$Ingenuity.Canonical.Pathways)))

ggplot(ipa_order[1:20,],aes(x=Ingenuity.Canonical.Pathways,y=X.log.p.value.))+
geom_bar(stat="identity")+
coord_flip()+
theme(text = element_text(size=15))

## celltype markers expression change after ifng treatment
markers10<-read.table("./data/10_top_markers_lnx.txt",sep="\t",header = T)
diff_markers10<-diff[as.character(markers10$gene),]
diff_markers10$x<-paste(row.names(diff_markers10),markers10[match(row.names(diff_markers10),markers10$gene),"cluster"],sep="#")
diff_markers10$X<-paste(markers10[match(row.names(diff_markers10),markers10$gene),"cluster"],row.names(diff_markers10),sep="#")

ggplot(na.omit(diff_markers10),aes(x=X,y=log2FoldChange))+
geom_bar(stat="identity")+
theme(text = element_text(size=15))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=8))

## markers10 heatmap Rna Ifng
count3<-na.omit(log2(count2[as.character(markers10$gene),]+0.1))
row.names(count3)<-paste(row.names(count3),markers10[match(row.names(count3),markers10$gene),"cluster"],sep="#")
pheatmap(count3[apply(count3,1,sd)!=0,],cluster_rows = T,cluster_cols = F,fontsize_row= 5,scale="row") 

##########13.05.2020

# anna rnaseq volcano plot label mhc2 genes
deg_aging<-read.csv("./data/Deseq2_AllGenes_young_M_organoids_oralVSold_M_organoids_oral.csv.title.csv",head=F)
colnames(deg_aging)<-colnames(deg_aging_oral)[1:8]
library(ggrepel)
volcano_label<-function(data,clname,cf=0.5) {
library(ggrepel)
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$log2FoldChange < cf)& data$log2FoldChange > (-cf)] <- "no"
data$sig[data$padj <= 0.05 & data$log2FoldChange >= cf] <- "up"
data$sig[data$padj <= 0.05 & data$log2FoldChange <= (-cf)] <- "down"
data$sign<-NA
data[data$sig!="no","sign"] <- row.names(data[data$sig!="no",])
#data[data$X %in% mhc2$genename,"sign"] <- data[data$X %in% mhc2$genename,"X"]
p<-ggplot(data[!is.na(data$sig),],aes(log2FoldChange,-1*log10(padj),color = sig))+geom_point(size=0.2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(clname)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)
print(p)
}
mhc2$genename<-as.character(mhc2$genename)
volcano_label(deg_aging,"DegAging")

# anna rnaseq volcano plot label deg genes

deg_aging_nona<-(na.omit(deg_aging))
volcano_label<-function(data,clname,cf=0.5) {
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$log2FoldChange < cf)& data$log2FoldChange > (-cf)] <- "no"
data$sig[data$padj <= 0.05 & data$log2FoldChange >= cf] <- "up"
data$sig[data$padj <= 0.05 & data$log2FoldChange <= (-cf)] <- "down"
data$sign<-NA
data[data$sig!="no","sign"] <- data[data$sig!="no","X"]
#data[data$X %in% mhc2$genename,"sign"] <- data[data$X %in% mhc2$genename,"X"]
p<-ggplot(data[!is.na(data$sig),],aes(log2FoldChange,-1*log10(padj),color = sig))+geom_point(size=0.2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(clname)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)
print(p)
}
volcano_label(deg_aging_nona,"DegAging")

## volcano without label
volcano<-function(data,clname,cf=0.5) {
data$sig[(data$padj > 0.05|data$padj=="NA")|(data$log2FoldChange < cf)& data$log2FoldChange > (-cf)] <- "no"
data$sig[data$padj <= 0.05 & data$log2FoldChange >= cf] <- "up"
data$sig[data$padj <= 0.05 & data$log2FoldChange <= (-cf)] <- "down"
data$sign<-NA
data[data$sig!="no","sign"] <- data[data$sig!="no","X"]
#data[data$X %in% mhc2$genename,"sign"] <- data[data$X %in% mhc2$genename,"X"]
p<-ggplot(data[!is.na(data$sig),],aes(log2FoldChange,-1*log10(padj),color = sig))+geom_point(size=0.2)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(clname)+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
#p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)
print(p)
}
volcano(deg_aging_nona,"DegAging")

## anna aging sig degs
cf=0.5
deg_aging_nona$sig[(deg_aging_nona$padj > 0.05|deg_aging_nona$padj=="NA")|(deg_aging_nona$log2FoldChange < cf)& deg_aging_nona$log2FoldChange > (-cf)] <- "no"
deg_aging_nona$sig[deg_aging_nona$padj <= 0.05 & deg_aging_nona$log2FoldChange >= cf] <- "up"
deg_aging_nona$sig[deg_aging_nona$padj <= 0.05 & deg_aging_nona$log2FoldChange <= (-cf)] <- "down"
rm(cf)
deg_aging_nona_sig<-deg_aging_nona[deg_aging_nona$sig!="no",]
# omid ifng sig degs
diff$sig<-"no"
diff[(diff$log2FoldChange>5)&(diff$padj<0.05),"sig"]<-"up"
diff[(diff$log2FoldChange<(-5))&(diff$padj<0.05),"sig"]<-"down"
diff_sig<-diff[diff$sig!="no",]
# venn plot
diff_sig$XXX<-paste(row.names(diff_sig),diff_sig$sig,sep="_")
deg_aging_nona_sig$YYY<-paste(deg_aging_nona_sig$X,deg_aging_nona_sig$sig,sep="_")
require(VennDiagram)
venn.diagram(list(deg_Ifng_sig_5=diff_sig$XXX,deg_aging_sig_0.5=deg_aging_nona_sig$YYY),filename = "test.png",col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
intersect(diff_sig$XXX,deg_aging_nona_sig$YYY)

deg_ifng_aging<-merge(diff_sig,deg_aging_nona_sig,by.x="row.names",by.y="X",all=T,suffixes = c(".ifng",".aging"))
write.csv(deg_ifng_aging,"./data/deg_ifng_aging.csv")

# anna aging ipa
# global setting for stringsAsFactors: options(stringsAsFactors = FALSE)
ipa_aging<-read.table("./data/ipa_deg_aging_male_organoids.pathways.txt",header = T,sep="\t",stringsAsFactors = FALSE,dec=",")
ipa_aging_sig<-ipa_aging[ipa_aging$X.log.p.value.>1.3,]
ipa_aging_sig_order<-ipa_aging_sig[order(-ipa_aging_sig$X.log.p.value.),]
ipa_ifng_aging<-merge(ipa_order,ipa_aging_sig_order,by.x="Ingenuity.Canonical.Pathways",by.y="Ingenuity.Canonical.Pathways",all=T,suffixes = c(".ifng",".aging"))
write.csv(ipa_ifng_aging,"./data/ipa_ifng_aging.csv")
venn.diagram(list(ipa_Ifng_sig_5=ipa_order$Ingenuity.Canonical.Pathways,ipa_aging_sig_0.5=ipa_aging_sig_order$Ingenuity.Canonical.Pathways),filename = "test.png",col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))

## crypts part
deg_aging_invivo<-read.csv("./data/Deseq2_AllGenes_young_M_crypts_oralVSold_M_crypts_oral.csv.title.csv",head=F)
colnames(deg_aging_invivo)<-colnames(deg_aging_oral)[1:8]
deg_aging_invivo_nona<-(na.omit(deg_aging_invivo))
cf=0.5
deg_aging_invivo_nona$sig[(deg_aging_invivo_nona$padj > 0.05|deg_aging_invivo_nona$padj=="NA")|(deg_aging_invivo_nona$log2FoldChange < cf)& deg_aging_invivo_nona$log2FoldChange > (-cf)] <- "no"
deg_aging_invivo_nona$sig[deg_aging_invivo_nona$padj <= 0.05 & deg_aging_invivo_nona$log2FoldChange >= cf] <- "up"
deg_aging_invivo_nona$sig[deg_aging_invivo_nona$padj <= 0.05 & deg_aging_invivo_nona$log2FoldChange <= (-cf)] <- "down"
rm(cf)
deg_aging_invivo_nona_sig<-deg_aging_invivo_nona[deg_aging_invivo_nona$sig!="no",]
deg_aging_invivo_nona_sig$YYY<-paste(deg_aging_invivo_nona_sig$X,deg_aging_invivo_nona_sig$sig,sep="_")
venn.diagram(list(deg_Ifng_sig_5=diff_sig$XXX,deg_aging_invivo_sig_0.5=deg_aging_invivo_nona_sig$YYY),filename = "venn_ifng_agingInvivo.png",col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral_crypts",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
deg_ifng_aging_invivo<-merge(diff_sig,deg_aging_invivo_nona_sig,by.x="row.names",by.y="X",all=T,suffixes = c(".ifng",".agingInvivo"))
write.csv(deg_ifng_aging_invivo,"./data_out/deg_ifng_aging_invivo.csv")

ipa_aging_invivo<-read.table("./data_out/ipa_deg_aging_male_crypts.pathways.txt",header = T,sep="\t",stringsAsFactors = FALSE,dec=",")
ipa_aging_invivo_sig<-ipa_aging_invivo[ipa_aging_invivo$X.log.p.value.>1.3,]
ipa_aging_invivo_sig_order<-ipa_aging_invivo_sig[order(-ipa_aging_invivo_sig$X.log.p.value.),]
ipa_ifng_aging_invivo<-merge(ipa_order,ipa_aging_invivo_sig_order,by.x="Ingenuity.Canonical.Pathways",by.y="Ingenuity.Canonical.Pathways",all=T,suffixes = c(".ifng",".aging"))
write.csv(ipa_ifng_aging_invivo,"./data/ipa_ifng_aging_invivo.csv")
venn.diagram(list(ipa_Ifng_sig_5=ipa_order$Ingenuity.Canonical.Pathways,ipa_aging_invivo_sig_0.5=na.omit(ipa_aging_invivo_sig_order$Ingenuity.Canonical.Pathways)),filename = "ipa_aging_invivo.png",col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))

## venn pdf
temp<-venn.diagram(list(ipa_Ifng_sig_5=ipa_order$Ingenuity.Canonical.Pathways,ipa_aging_invivo_sig_0.5=na.omit(ipa_aging_invivo_sig_order$Ingenuity.Canonical.Pathways)),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
pdf(file="venn.ipa_ifng_aging_invivo.pdf")
grid.draw(temp)
dev.off()

temp<-venn.diagram(list(deg_Ifng_sig_5=diff_sig$XXX,deg_aging_invivo_sig_0.5=deg_aging_invivo_nona_sig$YYY),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral_crypts",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
pdf(file="venn.deg_ifng_aging_invivo.pdf")
grid.draw(temp)
dev.off()

## venn pdf adjusted
temp<-venn.diagram(list(deg_Ifng_sig_5=diff_sig$XXX,deg_aging_invivo_sig_0.5=deg_aging_invivo_nona_sig$YYY),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral_crypts",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test = T,total.population = (nrow(diff)+nrow(deg_aging_invivo_nona)))
pdf(file="venn.deg_ifng_aging_invivo_2.pdf")
grid.draw(temp)
dev.off()
temp<-venn.diagram(list(deg_Ifng_sig_5=diff_sig$XXX,deg_aging_invivo_sig_0.5=deg_aging_invivo_nona_sig$YYY),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral_crypts",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test = T,total.population = (nrow(diff)+nrow(deg_aging_invivo_nona)),inverted=T)
pdf(file="venn.deg_ifng_aging_invivo_3.pdf")
grid.draw(temp)
dev.off()
temp<-venn.diagram(list(ipa_Ifng_sig_5=ipa_order$Ingenuity.Canonical.Pathways,ipa_aging_invivo_sig_0.5=na.omit(ipa_aging_invivo_sig_order$Ingenuity.Canonical.Pathways)),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test = T,total.population = (nrow(ipa_aging_invivo)+nrow(ipa)))
pdf(file="venn.ipa_ifng_aging_invivo_2.pdf")
grid.draw(temp)
dev.off()
## p value venn
venn_ipa_p<-function(ipa_aging,ipa_ifng,ipa_aging_sig,ipa_ifng_sig,fn_out){
    m=as.numeric(nrow(ipa_aging))
    n=as.numeric(nrow(ipa_ifng))
    q=as.numeric(length((intersect(na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)))))
    k=as.numeric(nrow(ipa_ifng_sig))
    print(paste(q,m,n,k,sep=","))
    print(paste(fn_out,"is:",sep=" "))
    print(phyper(q, m, n, k))
}
ipa_ifng<-read.table("./data/ipa_deg_ifng_logfc0.5_sig.txt",header = T,sep="\t")
ipa_ifng_sig<-read.table("./data/ipa_deg_ifng_logfc5_sig.txt",header = T,sep="\t")
venn_ipa_p(ipa_aging_invivo,ipa_ifng,ipa_aging_invivo_sig,ipa_ifng_sig,"rnaseq_ipa_p")
venn_deg_p<-function(deg_aging,deg_ifng,deg_aging_sig,deg_ifng_sig,fn_out){
    m=as.numeric(nrow(deg_aging))
    n=as.numeric(nrow(deg_ifng))
    q=as.numeric(length(intersect(na.omit(row.names(deg_ifng_sig)),na.omit(deg_aging_sig$X))))
    k=as.numeric(nrow(deg_ifng_sig))
    print(paste(q,m,n,k,sep=","))
    print(paste(fn_out,"is:",sep=" "))
    print(phyper(q, m, n, k))
}
venn_deg_p(deg_aging_invivo_nona,diff,deg_aging_invivo_nona,diff_sig,"rnaseq_deg_p")


===================================
===================================
===================================
## redo venn and p value
# ipa
ipa_ifng<-read.delim("./data/ipa_deg_ifng_logfc5.txt",header = T,sep="\t",stringsAsFactors = F)
ipa_ifng_sig<-read.delim("./data/ipa_deg_ifng_logfc5_sig.txt",header = T,sep="\t")
ipa_aging_invivo<-read.delim("./data/ipa_deg_aging_male_crypts.pathways.txt",header = T,sep="\t",stringsAsFactors = FALSE,dec=",")
ipa_aging_invivo_sig<-ipa_aging_invivo[ipa_aging_invivo$X.log.p.value.>1.3,]
venn_ipa_p<-function(ipa_ifng,ipa_aging,ipa_ifng_sig,ipa_aging_sig,fn_out){
    n=as.numeric(nrow(ipa_ifng))
    m=as.numeric(nrow(ipa_aging))
    q=as.numeric(length((intersect(na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),na.omit(ipa_aging_sig$Ingenuity.Canonical.Pathways)))))
    k=as.numeric(nrow(ipa_ifng_sig))
    print(paste("q= overlap in sigs","m= all aging","n= all ifng","k= sig ifng",sep=","))
    print(paste(q,m,n,k,sep=","))
    print(paste(fn_out,"is:",sep=" "))
    print(phyper(q, m, n, k))
    return(phyper(q, m, n, k))
}
pt<-venn_ipa_p(ipa_ifng,ipa_aging_invivo,ipa_ifng_sig,ipa_aging_invivo_sig,"rnaseq_ipa_p")
#temp<-venn.diagram(list(ipa_Ifng_sig=na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),ipa_aging_invivo_sig=na.omit(ipa_aging_invivo_sig$Ingenuity.Canonical.Pathways)),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = paste("M_oral_crypts, p=",pt,sep=" "),height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test = T,total.population = (nrow(ipa_ifng)+nrow(ipa_aging_invivo)))
temp<-venn.diagram(list(ipa_Ifng_sig=na.omit(ipa_ifng_sig$Ingenuity.Canonical.Pathways),ipa_aging_invivo_sig=na.omit(ipa_aging_invivo_sig$Ingenuity.Canonical.Pathways)),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = paste("M_oral_crypts, p=",pt,sep=" "),height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
pdf(file="venn.ipa_ifng_aging_invivo_redo.pdf")
grid.draw(temp)
dev.off()
ipa_ifng_aging<-merge(ipa_ifng_sig,ipa_aging_invivo_sig,by.x="Ingenuity.Canonical.Pathways",by.y="Ingenuity.Canonical.Pathways",all=T,suffixes = c(".ifng",".aging"))
write.csv(ipa_ifng_aging,"./data/ipa_ifng_aging_redo.csv")


# deg

### ignored if done before
deg_aging_invivo<-read.csv("./data/Deseq2_AllGenes_young_M_crypts_oralVSold_M_crypts_oral.csv.title.csv",head=F)
colnames(deg_aging_invivo)<-colnames(deg_aging_oral)[1:8]
deg_aging_invivo_nona<-(na.omit(deg_aging_invivo))
cf=0.5
deg_aging_invivo_nona$sig[(deg_aging_invivo_nona$padj > 0.05|deg_aging_invivo_nona$padj=="NA")|(deg_aging_invivo_nona$log2FoldChange < cf)& deg_aging_invivo_nona$log2FoldChange > (-cf)] <- "no"
deg_aging_invivo_nona$sig[deg_aging_invivo_nona$padj <= 0.05 & deg_aging_invivo_nona$log2FoldChange >= cf] <- "up"
deg_aging_invivo_nona$sig[deg_aging_invivo_nona$padj <= 0.05 & deg_aging_invivo_nona$log2FoldChange <= (-cf)] <- "down"
deg_aging_invivo_nona_sig<-deg_aging_invivo_nona[deg_aging_invivo_nona$sig!="no",]
deg_aging_invivo_nona_sig$YYY<-paste(deg_aging_invivo_nona_sig$X,deg_aging_invivo_nona_sig$sig,sep="_")

deg<-read.csv("./data/Deseq2_AllGenes_ConVSIfng_omid.csv",header = T,row.names = 1,check.names = F)
deg<-deg[!is.na(deg),]
diff<-deg[(!is.na(deg$log2FoldChange))&(!is.na(deg$padj)),]
diff$sig<-"no"
diff[(diff$log2FoldChange>5)&(diff$padj<0.05),"sig"]<-"up"
diff[(diff$log2FoldChange<(-5))&(diff$padj<0.05),"sig"]<-"down"
diff_sig<-diff[diff$sig!="no",]
diff_sig$XXX<-paste(row.names(diff_sig),diff_sig$sig,sep="_")
### ignored

venn_deg_p<-function(deg_ifng,deg_aging,deg_ifng_sig,deg_aging_sig,fn_out){
    n=as.numeric(nrow(deg_ifng))
    m=as.numeric(nrow(deg_aging))
    q=as.numeric(length(intersect(na.omit(deg_ifng_sig$XXX),na.omit(deg_aging_sig$YYY))))
    k=as.numeric(nrow(deg_ifng_sig))
    print(paste("q= overlap in sigs","m= all aging","n= all ifng","k= sig ifng",sep=","))
    print(paste(q,m,n,k,sep=","))
    print(paste(fn_out,"is:",sep=" "))
    print(phyper(q, m, n, k))
    return(phyper(q, m, n, k))
}
pt<-venn_deg_p(diff,deg_aging_invivo_nona,diff_sig,deg_aging_invivo_nona_sig,"rnaseq_deg_p")
#temp<-venn.diagram(list(deg_Ifng_sig_5=diff_sig$XXX,deg_aging_invivo_sig_0.5=deg_aging_invivo_nona_sig$YYY),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = "M_oral_crypts",height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20),hyper.test = T,total.population = (nrow(diff)+nrow(deg_aging_invivo_nona)))
temp<-venn.diagram(list(deg_Ifng_sig_5=diff_sig$XXX,deg_aging_invivo_sig_0.5=deg_aging_invivo_nona_sig$YYY),filename = NULL,col=c("#440154ff", 'blue'),fill = c(alpha("#440154ff",0.3), alpha('blue',0.3)),main = paste("M_oral_crypts p=",pt,sep=" "),height = 1000, width = 1000,main.cex = 0.8,cat.cex=0.4,cat.pos=c(-20,20))
pdf(file="venn.deg_ifng_aging_invivo_redo.pdf")
grid.draw(temp)
dev.off()
deg_ifng_aging_invivo_redo<-merge(diff_sig,deg_aging_invivo_nona_sig,by.x="XXX",by.y="YYY",all=T,suffixes = c(".ifng",".agingInvivo"))
write.csv(deg_ifng_aging_invivo_redo,"./data_output/deg_ifng_aging_invivo_redo.csv")

## check markers10 in anna data
# markers10
Counts_Norm_Markers10<-Counts_Norm[row.names(Counts_Norm)%in%markers10$gene,]
Counts_Norm_Markers10_cn<-Counts_Norm_Markers10[,cn]
colnames(Counts_Norm_Markers10_cn)<-gsub("geriatric","old",colnames(Counts_Norm_Markers10_cn))
anno_col<-data.frame(str_split_fixed(colnames(Counts_Norm_Markers10_cn),"_",5)[,2],str_split_fixed(colnames(Counts_Norm_Markers10_cn),"_",5)[,1],str_split_fixed(colnames(Counts_Norm_Markers10_cn),"_",5)[,3],row.names=colnames(Counts_Norm_Markers10_cn))
colnames(anno_col) <- c("Gender","Age","Invivo_Invitro")
xxx<-log2(Counts_Norm_Markers10_cn+0.1)
pheatmap(xxx[apply(xxx,1,sd)!=0,],cluster_cols = F,fontsize_row=4,legend_labels=c("log"),annotation_col = anno_col,border_color = F,fontsize = 6,main="Markers10 genes expression in bulk RNAseq, log2(exp+1)",scale = "row")
# 
cn2<-c("young_F_crypts_oral_rep1",
        "young_F_crypts_oral_rep2",
        "young_F_crypts_oral_rep3",
        "young_F_crypts_oral_rep4",
        "young_M_crypts_oral_rep1",
        "young_M_crypts_oral_rep2",
        "young_M_crypts_oral_rep3",
        "young_M_crypts_oral_rep4",
        "geriatric_F_crypts_oral_rep1",
        "geriatric_F_crypts_oral_rep2",
        "geriatric_F_crypts_oral_rep3",
        "geriatric_F_crypts_oral_rep4",
        "old_M_crypts_oral_rep1",
        "old_M_crypts_oral_rep2",
        "old_M_crypts_oral_rep3",
        "old_M_crypts_oral_rep4")
Counts_Norm_Markers10_cn2<-Counts_Norm_Markers10[,cn2]
colnames(Counts_Norm_Markers10_cn2)<-gsub("geriatric","old",colnames(Counts_Norm_Markers10_cn2))
anno_col<-data.frame(str_split_fixed(colnames(Counts_Norm_Markers10_cn2),"_",5)[,2],str_split_fixed(colnames(Counts_Norm_Markers10_cn2),"_",5)[,1],str_split_fixed(colnames(Counts_Norm_Markers10_cn2),"_",5)[,3],row.names=colnames(Counts_Norm_Markers10_cn2))
colnames(anno_col) <- c("Gender","Age","Invivo_Invitro")
xxxx<-log2(Counts_Norm_Markers10_cn2+0.1)
pheatmap(xxxx[apply(xxxx,1,sd)!=0,],cluster_cols = F,fontsize_row=4,legend_labels=c("log"),annotation_col = anno_col,border_color = F,fontsize = 6,main="Markers10 genes expression in bulk RNAseq, log2(exp+1)",scale = "row")
rownames(xxxx)<-(markers10[match(rownames(xxxx),markers10$gene),"gene_cluster"])
pheatmap(xxxx[apply(xxxx,1,sd)!=0,],cluster_cols = F,fontsize_row=4,legend_labels=c("log"),annotation_col = anno_col,border_color = F,fontsize = 6,main="Markers10 genes expression in bulk RNAseq, log2(exp+1)",scale = "row")
## boxplot markers10
library(reshape2)
xxxxx$gene<-row.names(xxxxx)
xxxxx_long<-melt(xxxxx,
                 id.vars = "gene",
                 measure.vars = c( "young_M_crypts_oral_rep1",
                                   "young_M_crypts_oral_rep2",
                                   "young_M_crypts_oral_rep3",
                                   "young_M_crypts_oral_rep4",
                                   "old_M_crypts_oral_rep1",
                                   "old_M_crypts_oral_rep2",
                                   "old_M_crypts_oral_rep3",
                                   "old_M_crypts_oral_rep4"),
                 variable.name='sample',
                 value.name='exp_log2m')
xxxxx_long$group<-(str_split_fixed(xxxxx_long$sample,"_rep",2)[,1])
#ggplot(xxxxx_long, aes(x=gene, y=exp_log2m, fill=group)) + geom_boxplot()+theme(axis.text=element_text(size=4,face = "bold"),axis.title.x=element_text(size=4))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
xxxxx_long$group<-factor(xxxxx_long$group,levels = c("young_M_crypts_oral","old_M_crypts_oral"))
library(ggpubr)
ggboxplot(xxxxx_long, x="gene", y="exp_log2m", color="group",palette = "jco",outlier.size=0.2)+theme(axis.text=element_text(size=4,face = "bold"),axis.title.x=element_text(size=4))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ stat_compare_means(aes(group = group),label = "p.signif",size=2)
xxxxx_long$celltype<-str_split_fixed(xxxxx_long$gene,":",2)[,2]
cts<-data.frame(ct=c("E_Stem","D_TA","C_Early_Enterocytes_Prog","B_Late_Enterocytes_Prog","A_Enterocyte","F_Secretory_precursos","I_Enteroendocrine","G_Goblet","J_Tuft"."H_Paneth"),num=0:9)
for (i in 1:nrow(xxxxx_long)){xxxxx_long$num[i]<-cts[cts$ct==xxxxx_long$celltype[i],"num"]}
xxxxx_long2<-xxxxx_long[order(xxxxx_long$num),]
ggboxplot(xxxxx_long2, x="gene", y="exp_log2m", color="group",palette = "jco",outlier.size=0.2)+theme(axis.text=element_text(size=4,face = "bold"),axis.title.x=element_text(size=4))+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+ stat_compare_means(aes(group = group),label = "p.signif",size=2)+scale_color_manual(values = c("gray","gray0"))

## Anna RNAseq enrichment analysis
## RnaAnna_old_GSEA_v2.pdf
options(stringsAsFactors = F)
geneset<-read.csv("./data/geneset_Kim_Haber_etal.csv",head=T)
Counts_Norm_Geneset<-Counts_Norm[row.names(Counts_Norm)%in%geneset$gene,]
Counts_Norm_Geneset_cn3<-Counts_Norm_Geneset[,cn3]

csel<-Counts_Norm_Geneset_cn3[apply(Counts_Norm_Geneset_cn3 > 0, 1, sum)>ncol(Counts_Norm_Geneset_cn3)*0.2,]
csel<-t(scale(t(csel)))
csel<-csel*gsel$Direction
csel<-as.data.frame(csel)

gsel<-geneset[geneset$gene %in% rownames(csel) & !is.na(geneset$gene),]
gcsel<-merge(gsel,csel,by.x="gene",by.y="row.names",all.x = T)
gcsel_long<-melt(gcsel,
                 id.vars = "setName",
                 measure.vars = c( "young_M_crypts_oral_rep1",
                                   "young_M_crypts_oral_rep2",
                                   "young_M_crypts_oral_rep3",
                                   "young_M_crypts_oral_rep4",
                                   "old_M_crypts_oral_rep1",
                                   "old_M_crypts_oral_rep2",
                                   "old_M_crypts_oral_rep3",
                                   "old_M_crypts_oral_rep4"),
                 variable.name='setNames',
                 value.name='exp_scaled')
gcsel_long$group<-(str_split_fixed(gcsel_long$setNames,"_",2)[,1])
gcsel_long$group<-factor(gcsel_long$group,levels = c("young","old"))
p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group", palette = "jco",
facet.by = "setName", short.panel.labs = FALSE)
p + stat_compare_means(label = "p.format")

## ifng RNAseq gsea

count_geneset<-count[row.names(count)%in%geneset$gene,]
count_geneset_cn3<-count_geneset[,1:8]
gsea<-function(count_geneset_cn3,geneset) {
require(ggpubr)
require(reshape2)
colnames(geneset)<-c("setName","gene","Direction")
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
gcsel_long$group<-factor(gcsel_long$group,levels = c("Con","Ifng"))
p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group", palette = "jco",facet.by = "setName", nrow=length(gcsel_long$group),short.panel.labs = FALSE)
p + stat_compare_means(label = "p.format",method="t.test")+ggtitle("t.test_pairedF") + scale_color_manual(values=c("gray","gray0"))
}
gsea(count_geneset_cn3,geneset)


### anna M oral mhc2
mo<-c(
     "young_M_crypts_oral_rep1",
     "young_M_crypts_oral_rep2",
     "young_M_crypts_oral_rep3",
     "young_M_crypts_oral_rep4",
     "old_M_crypts_oral_rep1",
     "old_M_crypts_oral_rep2",
     "old_M_crypts_oral_rep3",
     "old_M_crypts_oral_rep4",
     "young_M_organoids_oral_rep1",
     "young_M_organoids_oral_rep2",
     "young_M_organoids_oral_rep3",
     "young_M_organoids_oral_rep4",
     "old_M_organoids_oral_rep1",
     "old_M_organoids_oral_rep2",
     "old_M_organoids_oral_rep3",
     "old_M_organoids_oral_rep4")
Counts_Norm_MHCII_mo<-Counts_Norm_MHCII[,mo]
Counts_Norm_MHCII_mo_log<-log2(Counts_Norm_MHCII_mo+0.1)
anno_col<-data.frame(factor(str_split_fixed(colnames(Counts_Norm_MHCII_mo),"_",5)[,1],levels = c("young","old")),str_split_fixed(colnames(Counts_Norm_MHCII_mo),"_",5)[,3],row.names=colnames(Counts_Norm_MHCII_mo))
colnames(anno_col) <- c("Age","Invivo_Invitro")
pheatmap(Counts_Norm_MHCII_mo_log[apply( Counts_Norm_MHCII_mo_log,1,sd)!=0,],cluster_cols = F,fontsize_row=4,legend_labels=c("log"),annotation_col = anno_col,border_color = F,fontsize = 6,main="MHC2 genes expression in bulk RNAseq, log2(exp+0.1)",scale = "row")


## anna M oral GSEA
Counts_Norm_Geneset_mo<-Counts_Norm_Geneset[,mo]
csel<-Counts_Norm_Geneset_mo[apply(Counts_Norm_Geneset_mo > 0, 1, sum)>ncol(Counts_Norm_Geneset_mo)*0.2,]
csel<-t(scale(t(csel)))
csel<-as.data.frame(csel)

gsel<-geneset[geneset$gene %in% rownames(csel) & !is.na(geneset$gene),]
csel<-csel*gsel$Direction
gcsel<-merge(gsel,csel,by.x="gene",by.y="row.names",all.x = T)

gcsel_long<-melt(gcsel,
                 id.vars = "setName",
                 measure.vars = mo,
                 variable.name='setNames',
                 value.name='exp_scaled')
gcsel_long$group<-(str_split_fixed(gcsel_long$setNames,"_",2)[,1])
gcsel_long$group<-factor(gcsel_long$group,levels = c("young_M_crypts_oral","old_M_crypts_oral","young_M_organoids_oral","old_M_organoids_oral"))
p <- ggboxplot(gcsel_long, x = "group", y = "exp_scaled",color = "group", palette = "jco",
facet.by = "setName", short.panel.labs = FALSE,nrow=length(gcsel_long$group))
p + stat_compare_means(label = "p.format")

my_comparisons <- list( c("young_M_crypts_oral","old_M_crypts_oral"), c("young_M_organoids_oral","old_M_organoids_oral"), c("young_M_crypts_oral","old_M_crypts_oral","young_M_organoids_oral","old_M_organoids_oral") )
p + stat_compare_means(label = "p.format",comparisons = my_comparisons)
p + stat_compare_means(label = "p.format",comparisons = my_comparisons)+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+ scale_color_manual(values=c("lightpink", "hotpink", "lightblue1","deepskyblue1")) + theme(legend.position="right")


## version 2
## GSEA only keep old M 
gcsel_long_old<-gcsel_long[gcsel_long$group %in% c("old_M_crypts_oral","old_M_organoids_oral"),]
gcsel_long_old$group<-droplevels(gcsel_long_old$group)
library(ggpubr)
p1<-ggboxplot(gcsel_long_old, x = "group", y = "exp_scaled",color = "group", palette = "jco",facet.by = "setName", short.panel.labs = FALSE,nrow=length(gcsel_long$group))
p1 + stat_compare_means(label = "p.format")+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+ scale_color_manual(values=c("gray", "gray0")) + theme(legend.position="right")
## antigen presentation pathway's sigDE genes in old M bulk RNAseq
mhc2_sigDEG<-read.csv("./data/KEGG_mhc.genes.sigDEG",header = T)
Counts_Norm_MHCII_mo_log_old_sigDEG<-Counts_Norm_MHCII_mo_log[row.names(Counts_Norm_MHCII_mo_log) %in% as.character(mhc2_sigDEG$gene),grepl("old",colnames(Counts_Norm_MHCII_mo_log))]
library(stringr)
library(pheatmap)
anno_col_old<-data.frame(factor(str_split_fixed(colnames(Counts_Norm_MHCII_mo_log_old_sigDEG),"_",5)[,3],levels = c("crypts","organoids")),row.names=colnames(Counts_Norm_MHCII_mo_log_old_sigDEG))
colnames(anno_col_old)<-"Invivo_Invitro"
annotation_colors = list(
Invivo_Invitro = c(crypts="gray", organoids="gray0")
)
pheatmap(Counts_Norm_MHCII_mo_log_old_sigDEG[apply(Counts_Norm_MHCII_mo_log_old_sigDEG,1,sd)!=0,],cluster_cols = F,fontsize_row=6,legend_labels=c("log"),annotation_col = anno_col_old,border_color = F,annotation_colors =annotation_colors,fontsize = 6,main="Antigen-presentation pathway's sigDE genes expression in bulk RNAseq, log2(exp+0.1)",scale = "row")
## antigen presentation pathway's Dresden genes in old M bulk RNAseq
mhc2_mahdi<-read.csv("./data/mhc2_mahdi.genes",header = T,stringsAsFactors = F)
Counts_Norm<-read.csv("./data/Counts_norm_DESeq2_AllGenes.csv",head=T,row.names = 1)

cn_old<-c(
"old_M_crypts_oral_rep1",
"old_M_crypts_oral_rep2",
"old_M_crypts_oral_rep3",
"old_M_crypts_oral_rep4",
"old_M_organoids_oral_rep1",
"old_M_organoids_oral_rep2",
"old_M_organoids_oral_rep3",
"old_M_organoids_oral_rep4")
Counts_Norm_MHCII_mahdi<-Counts_Norm[row.names(Counts_Norm)%in%mhc2_mahdi$gene,]
Counts_Norm_MHCII_mahdi_cno<-Counts_Norm_MHCII_mahdi[,cn_old]
Counts_Norm_MHCII_mahdi_cno_log<-log2(Counts_Norm_MHCII_mahdi_cno+0.1)
anno_col_mahdi_old<-data.frame(factor(str_split_fixed(colnames(Counts_Norm_MHCII_mahdi_cno_log),"_",5)[,3],levels = c("crypts","organoids")),row.names=colnames(Counts_Norm_MHCII_mahdi_cno_log))
colnames(anno_col_mahdi_old)<-"Invivo_Invitro"
pheatmap(Counts_Norm_MHCII_mahdi_cno_log[apply(Counts_Norm_MHCII_mahdi_cno_log,1,sd)!=0,],cluster_cols = F,fontsize_row=6,legend_labels=c("log"),annotation_col = anno_col_mahdi_old,border_color = F,annotation_colors =annotation_colors,fontsize = 6,main="Antigen-presentation pathway's Dresden genes expression in bulk RNAseq, log2(exp+0.1)",scale = "row")
## remove not-sigDE genes in old M: Psmb6
tmp<-Counts_Norm_MHCII_mahdi_cno_log[row.names(Counts_Norm_MHCII_mahdi_cno_log)!="Psmb6",]
pheatmap(tmp[apply(tmp,1,sd)!=0,],cluster_cols = F,fontsize_row=6,legend_labels=c("log"),annotation_col = anno_col_mahdi_old,border_color = F,annotation_colors =annotation_colors,fontsize = 6,main="Antigen-presentation pathway's Dresden genes expression in bulk RNAseq, log2(exp+0.1)",scale = "row")
save.image("./data/Point3.RData")

=========================================
=========================================
#### Modify_1 ###########################
load("./data/Point3.RData")
library(pheatmap)
## orginal one: Counts_Norm_MHCII_mahdi_cno_log_sigDEG_v2.pdf
pheatmap(tmp[apply(tmp,1,sd)!=0,],cluster_cols = F,fontsize_row=6,legend_labels=c("log"),annotation_col = anno_col_mahdi_old,border_color = F,annotation_colors =annotation_colors,fontsize = 6,main="Antigen-presentation pathway's Dresden genes expression in bulk RNAseq, log2(exp+0.1)",scale = "row")
## display all mahdi's mhc2 genes by his order: RnaAnna_mhc2_mahdi.pdf
mhc2_mahdi<-c("H2-M3","H2-T23","Tap2","Tapbp","Cd74","H2-Aa","H2-Eb1","H2-Ab1","Ciita","H2-DMb1","H2-DMa","H2-Q7","Psmb8","B2m","Tap1","Calr","Pdia3","Psmb6")
tmp<-Counts_Norm_MHCII_mahdi_cno_log[mhc2_mahdi,]
pheatmap(tmp[apply(tmp,1,sd)!=0,],cluster_cols = F,cluster_rows = F, fontsize_row=6,legend_labels=c("log"),annotation_col = anno_col_mahdi_old,border_color = F,annotation_colors =annotation_colors,fontsize = 6,main="Antigen-presentation pathway's Dresden genes expression in bulk RNAseq, log2(exp+0.1)",scale = "row")

==========================================
==========================================
#### Francesco's comments
load("./data/Point3.RData")
library(ggpubr)
Counts_Norm_Markers10_mo<-Counts_Norm_Markers10[,mo]
csel_Markers10<-Counts_Norm_Markers10_mo[apply(Counts_Norm_Markers10_mo > 0, 1, sum)>ncol(Counts_Norm_Markers10_mo)*0.2,]
csel_Markers10<-t(scale(t(csel_Markers10)))
csel_Markers10<-as.data.frame(csel_Markers10)
gsel_Markers10<-markers10[markers10$gene %in% rownames(csel_Markers10) & !is.na(markers10$gene),]
csel_Markers10<-csel_Markers10*1
gcsel_Markers10<-merge(gsel_Markers10,csel_Markers10,by.x="gene",by.y="row.names",all.x = T)
gcsel_Markers10_long<-reshape2::melt(gcsel_Markers10,
id.vars = "cluster",
measure.vars = mo,
variable.name='cluster',
value.name='exp_scaled')

library(stringr)
colnames(gcsel_Markers10_long)[1]<-"setName"
gcsel_Markers10_long$group<-(str_split_fixed(gcsel_Markers10_long$cluster,"_rep",2)[,1])
gcsel_Markers10_long$group<-factor(gcsel_Markers10_long$group,levels = c("young_M_crypts_oral","old_M_crypts_oral","young_M_organoids_oral","old_M_organoids_oral"))
gcsel_Markers10_long_old<-gcsel_Markers10_long[gcsel_Markers10_long$group %in% c("old_M_crypts_oral","old_M_organoids_oral"),]
gcsel_Markers10_long_old$group<-droplevels(gcsel_Markers10_long_old$group)

## RnaAnna_old_M_GSEA_Markers10.pdf
p1_Markers10<-ggboxplot(gcsel_Markers10_long_old, x = "group", y = "exp_scaled",color = "group", palette = "jco",facet.by = "setName", short.panel.labs = FALSE,nrow=length(gcsel_Markers10_long$group))
p1_Markers10 + stat_compare_means(label = "p.format")+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+ scale_color_manual(values=c("gray", "gray0")) + theme(legend.position="right")
## RnaAnna_old_M_GSEA_Markers10s.pdf
ggboxplot(gcsel_Markers10_long, x = "group", y = "exp_scaled",color = "group", palette = "jco",facet.by = "setName", short.panel.labs = FALSE,nrow=length(gcsel_Markers10_long$group)) + stat_compare_means(label = "p.format",comparisons = my_comparisons)+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+ scale_color_manual(values=c("lightpink", "hotpink", "lightblue1","deepskyblue1")) + theme(legend.position="right")
saveRDS(gcsel_Markers10_long,"./data/gcsel_Markers10_long.rds")

### Check with all the cell type markers from Dresden's data
MarkersMore<-read.csv("./data/CellTypeMarkers_GseaInput.csv",header = T,stringsAsFactors = F)
colnames(MarkersMore)[2]<-"gene"
load("./data/Point3.RData")
Counts_Norm<-read.csv("./data/Counts_norm_DESeq2_AllGenes.csv",head=T,row.names = 1)
Counts_Norm_MarkersMore<-Counts_Norm[row.names(Counts_Norm)%in%MarkersMore$gene,]
Counts_Norm_MarkersMore_mo<-Counts_Norm_MarkersMore[,mo]
csel_MarkersMore<-Counts_Norm_MarkersMore_mo[apply(Counts_Norm_MarkersMore_mo > 0, 1, sum)>ncol(Counts_Norm_MarkersMore_mo)*0.2,]
csel_MarkersMore<-t(scale(t(csel_MarkersMore)))
csel_MarkersMore<-as.data.frame(csel_MarkersMore)
gsel_MarkersMore<-MarkersMore[MarkersMore$gene %in% rownames(csel_MarkersMore) & !is.na(MarkersMore$gene),]
csel_MarkersMore<-csel_MarkersMore*1
gcsel_MarkersMore<-merge(gsel_MarkersMore,csel_MarkersMore,by.x="gene",by.y="row.names",all.x = T)

gcsel_MarkersMore_long<-reshape2::melt(gcsel_MarkersMore,
id.vars = "Cluster",
measure.vars = mo,
variable.name='Cluster',
value.name='exp_scaled')
library(stringr)
colnames(gcsel_MarkersMore_long)[1]<-"setName"
gcsel_MarkersMore_long$group<-(str_split_fixed(gcsel_MarkersMore_long$Cluster,"_rep",2)[,1])
gcsel_MarkersMore_long$group<-factor(gcsel_MarkersMore_long$group,levels = c("young_M_crypts_oral","old_M_crypts_oral","young_M_organoids_oral","old_M_organoids_oral"))
gcsel_MarkersMore_long_old<-gcsel_MarkersMore_long[gcsel_MarkersMore_long$group %in% c("old_M_crypts_oral","old_M_organoids_oral"),]
gcsel_MarkersMore_long_old$group<-droplevels(gcsel_MarkersMore_long_old$group)
library(ggpubr)
## RnaAnna_old_M_GSEA_MarkersMore.pdf RnaAnna_old_M_GSEA_MarkersMores.pdf
ggboxplot(gcsel_MarkersMore_long_old, x = "group", y = "exp_scaled",color = "group", palette = "jco",facet.by = "setName", short.panel.labs = FALSE,nrow=length(gcsel_MarkersMore_long$group))+ stat_compare_means(label = "p.format")+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+ scale_color_manual(values=c("gray", "gray0")) + theme(legend.position="right")
ggboxplot(gcsel_MarkersMore_long, x = "group", y = "exp_scaled",color = "group", palette = "jco",facet.by = "setName", short.panel.labs = FALSE,nrow=length(gcsel_MarkersMore_long$group)) + stat_compare_means(label = "p.format",comparisons = my_comparisons)+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+ scale_color_manual(values=c("lightpink", "hotpink", "lightblue1","deepskyblue1")) + theme(legend.position="right")

### crypts and organoids
csel_MarkersMore<-Counts_Norm_MarkersMore[apply(Counts_Norm_MarkersMore > 0, 1, sum)>ncol(Counts_Norm_MarkersMore)*0.2,]
csel_MarkersMore<-t(scale(t(csel_MarkersMore)))
csel_MarkersMore<-as.data.frame(csel_MarkersMore)
gsel_MarkersMore<-MarkersMore[MarkersMore$gene %in% rownames(csel_MarkersMore) & !is.na(MarkersMore$gene),]
csel_MarkersMore<-csel_MarkersMore*1
gcsel_MarkersMore<-merge(gsel_MarkersMore,csel_MarkersMore,by.x="gene",by.y="row.names",all.x = T)

gcsel_MarkersMore_long<-reshape2::melt(gcsel_MarkersMore,
id.vars = "Cluster",
measure.vars = colnames(gcsel_MarkersMore)[c(3:98)],  ## change from some specific samples to all samples
variable.name='Cluster',
value.name='exp_scaled')

colnames(gcsel_MarkersMore_long)[1]<-"setName"
gcsel_MarkersMore_long$group<-(str_split_fixed(gcsel_MarkersMore_long$Cluster,"_rep",2)[,1])
gcsel_MarkersMore_long$group<-paste(sapply(str_split(gcsel_MarkersMore_long$group,"_"),'[',1),sapply(str_split(gcsel_MarkersMore_long$group,"_"),'[',3),sep="_")
gcsel_MarkersMore_long$group<-gsub("geriatric","old",gcsel_MarkersMore_long$group)
gcsel_MarkersMore_long$group<-factor(gcsel_MarkersMore_long$group,levels = c("young_crypts","old_crypts","young_organoids","old_organoids"))
### RnaAnna_All_GSEA_MarkersMores.pdf  RnaAnna_All_GSEA_MarkersMores_2.pdf
my_comparisons_all<-list(c("young_crypts","old_crypts"),c("young_organoids","old_organoids"),c("young_crypts","young_organoids"),c("old_crypts","old_organoids"))
ggboxplot(gcsel_MarkersMore_long, x = "group", y = "exp_scaled",color = "group", palette = "jco",facet.by = "setName", short.panel.labs = FALSE,nrow=length(gcsel_MarkersMore_long$group)) + stat_compare_means(label = "p.format",comparisons = my_comparisons_all)+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+ scale_color_manual(values=c("lightpink", "hotpink", "lightblue1","deepskyblue1")) + theme(legend.position="right")
