## Author: Jing


setwd("./")
stat1<-read.csv("./data/stat1.csv",header = T)
stat1$types<-str_split_fixed(stat1$type,"_wilcox",2)[,1]
stat1$geneID<-factor(stat1$geneID,levels = c("Stat1","Cd74"))
ggplot(data=stat1, aes(x=types, y=-avg_logFC,fill=geneID)) + geom_bar(stat="identity", position=position_dodge()) +coord_flip() + scale_fill_grey()+guides(fill = guide_legend(reverse = TRUE))

## volcano
volcano<-function(data,clname,cf=0.5) {
    require(ggrepel)
    data$log2FoldChange<-(-data$avg_logFC)
    data$sig[(data$p_val_adj > 0.05|data$p_val_adj=="NA")|(data$log2FoldChange < cf)& data$log2FoldChange > (-cf)] <- "no"
    data$sig[data$p_val_adj <= 0.05 & data$log2FoldChange >= cf] <- "up"
    data$sig[data$p_val_adj <= 0.05 & data$log2FoldChange <= (-cf)] <- "down"
    data$sign<-NA
    data[data$geneID %in% c("Cd74","Stat1"),"sign"] <- as.character(data[data$geneID %in% c("Cd74","Stat1"),"geneID"])
    p<-ggplot(data[!is.na(data$sig),],aes(log2FoldChange,-1*log10(p_val_adj),color = sig))+geom_point(size=1)+labs(x="-log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(colnames(data)[1])+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
    p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.5, "lines"), show.legend = F, size = 3)
    print(p)
}
o_ss_cd74<-read.csv("./data/Stem.Cd74n_Control_vs_Stem_Cd74p_IFNg_wilcox_lfc0.1.csv",header = T)
volcano(o_ss_cd74)
o_st_cd74<-read.csv("./data/Stem.Cd74n_Control_vs_TA_Cd74p_IFNg_wilcox_lfc0.1.csv",header = T)
volcano(o_st_cd74)
d_ss_cd74<-read.csv("./data/Stem.Cd74n_Young_vs_Cd74p_Old_wilcox_lfc0.1.csv",header = T )
volcano(d_ss_cd74)

volcano<-function(data,clname,cf=0.5) {
    require(ggrepel)
    data$log2FoldChange<-(-data$avg_logFC)
    data$sig[(data$p_val_adj > 0.05|data$p_val_adj=="NA")|(data$log2FoldChange < cf)& data$log2FoldChange > (-cf)] <- "no"
    data$sig[data$p_val_adj <= 0.05 & data$log2FoldChange >= cf] <- "up"
    data$sig[data$p_val_adj <= 0.05 & data$log2FoldChange <= (-cf)] <- "down"
    data$sign<-NA
    data[data$geneID %in% c("Cd74","Stat1"),"sign"] <- as.character(data[data$geneID %in% c("Cd74","Stat1"),"geneID"])
    p<-ggplot(data[!is.na(data$sig),],aes(log2FoldChange,-1*log10(p_val_adj),color = sig))+geom_point(size=0.5)+labs(x="log2(FoldChange)",y="-log10(FDR)")+ scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-cf,cf),linetype=4)+ggtitle(colnames(data)[1])+theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
    p<-p+geom_text_repel(aes(label = sign), box.padding = unit(0.5, "lines"), show.legend = F, size = 3)
    print(p)
}
volcano(d_ss_cd74,cf=0.5)
volcano(o_ss_cd74,cf=0.5)
volcano(o_st_cd74,cf=0.5)

######## Stat1 Vln in organoids

immune.combined<-readRDS("./data/immune.combined.rds")
# Vln_Stat1_StemCD74.pdf
VlnPlot(object = subset(x = immune.combined, subset = celltype.treatment.Cd74 == c("E_Stem_cluster3_sc_control_negative_Cd74","E_Stem_cluster3_sc_ifnr_positive_Cd74")), features = c('Cd74','Stat1'), group.by = 'celltype.treatment.Cd74')
# HalfVln_Stat1_AllCellTypes.pdf
Idents(immune.combined)<-factor(Idents(immune.combined),levels = c("E_Stem_cluster3_sc_control","C_Early_Enterocytes_Prog_D_TA_cluster0_sc_control","A_Enterocyte_B_Late_Enterocytes_Prog_cluster1_sc_control","A_Enterocyte_cluster2_sc_control","I_Enteroendocrine_cluster4_sc_control","G_Goblet_H_Paneth_cluster5_sc_control","J_Tuft_cluster6_sc_control","E_Stem_cluster3_sc_ifnr","C_Early_Enterocytes_Prog_D_TA_cluster0_sc_ifnr","A_Enterocyte_B_Late_Enterocytes_Prog_cluster1_sc_ifnr","A_Enterocyte_cluster2_sc_ifnr","I_Enteroendocrine_cluster4_sc_ifnr","G_Goblet_H_Paneth_cluster5_sc_ifnr","J_Tuft_cluster6_sc_ifnr"))
immune.combined$celltype<-factor(immune.combined$celltype,levels = c("E_Stem_cluster3","C_Early_Enterocytes_Prog_D_TA_cluster0","A_Enterocyte_B_Late_Enterocytes_Prog_cluster1","A_Enterocyte_cluster2","I_Enteroendocrine_cluster4","G_Goblet_H_Paneth_cluster5","J_Tuft_cluster6"))
immune.combined$treatment<-factor(immune.combined$treatment,levels = c("sc_control","sc_ifnr"))
VlnPlot(object = immune.combined, features = c('Stat1'),group.by = "celltype",split.by = "treatment",cols=c("forestgreen", "orange"),pt.size = 0.2)+aes(color=immune.combined$treatment)+ scale_color_manual(values=c("forestgreen", "orange"))
# Vln_Stat1_AllCellTypes.pdf Vln_Stat1_AllCellTypes_2.pdf
Idents(immune.combined)<-factor(Idents(immune.combined),levels = c("E_Stem_cluster3_sc_control","E_Stem_cluster3_sc_ifnr","C_Early_Enterocytes_Prog_D_TA_cluster0_sc_control","C_Early_Enterocytes_Prog_D_TA_cluster0_sc_ifnr","A_Enterocyte_B_Late_Enterocytes_Prog_cluster1_sc_control","A_Enterocyte_B_Late_Enterocytes_Prog_cluster1_sc_ifnr","A_Enterocyte_cluster2_sc_control","A_Enterocyte_cluster2_sc_ifnr","I_Enteroendocrine_cluster4_sc_control","I_Enteroendocrine_cluster4_sc_ifnr","G_Goblet_H_Paneth_cluster5_sc_control","G_Goblet_H_Paneth_cluster5_sc_ifnr","J_Tuft_cluster6_sc_control","J_Tuft_cluster6_sc_ifnr"))
ColorScheme<-c("#9489be","#9489be","#f07570","#f07570","#c19938","#c19938","#33b792","#33b792","#23b5e7","#23b5e7","#56b24e","#56b24e","#d36ea9","#d36ea9")
VlnPlot(object = immune.combined, features = c('Stat1'),pt.size = 0,cols = ColorScheme)+NoLegend()
VlnPlot(object = immune.combined, features = c('Stat1'),pt.size = 0.2,cols = ColorScheme)+NoLegend()


########## Stat1 vln in crypts
######### aged in vivo
age.combined<-readRDS("./data/age.combined.rds")
DefaultAssay(age.combined) <- "RNA"
## Vln_Stat1_StemCD74_Age.pdf
age.combined.subset<-subset(x = age.combined, subset = Celltypes_Age_Cd74 == c("Stem_Young_Cd74.negative","Stem_Old_Cd74.positive"))
age.combined.subset$Celltypes_Age_Cd74<-factor(age.combined.subset$Celltypes_Age_Cd74,levels = c("Stem_Young_Cd74.negative","Stem_Old_Cd74.positive"))
VlnPlot(object = age.combined.subset, features = c('Cd74','Stat1'), group.by = 'Celltypes_Age_Cd74')
## Vln_Stat1_AllCellTypes_Age.pdf Vln_Stat1_AllCellTypes_Age_2.pdf
age.combined$celltypes_age<-factor(age.combined$celltypes_age,levels = c("Stem_Young","Stem_Old","TA_Young","TA_Old","Secretory.Progenitors_Young","Secretory.Progenitors_Old","Early.Enterocytes.Progenitors_Young","Early.Enterocytes.Progenitors_Old","Late.Enterocytes.Progenitors_Young","Late.Enterocytes.Progenitors_Old","Enterocyte_Young","Enterocyte_Old","Enteroendocrine_Young","Enteroendocrine_Old","Goblet_Young","Goblet_Old","Paneth_Young","Paneth_Old","Tuft_Young","Tuft_Old"))
VlnPlot(object = age.combined, features = c('Stat1'),pt.size = 0.2,cols = rep("black",length(unique(age.combined$celltypes_age))),group.by = "celltypes_age")+NoLegend()
VlnPlot(object = age.combined, features = c('Stat1'),pt.size = 0,cols = rep("black",length(unique(age.combined$celltypes_age))),group.by = "celltypes_age")+NoLegend()
