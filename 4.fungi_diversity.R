library(stringr)
library(ggplot2)
library(reshape)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(gridExtra)
library(vegan)
library(dplyr)
library(car)
library(devtools)

setwd("/home/user/")

# load plotting styles
source("Plotting_styles.R")
# load custom functions
source_url("https://github.com/mchialva/myNGS_tools/blob/master/TaxR.R?raw=TRUE")

#### load Sample metadata and normalized reads counts ####
load("my_colData")
load("Fungi_norm")

Fungi_norm[Fungi_norm==0]<-NA

#### Cut LCA table at phylum level ####
Fungi_phylum<-LCA_level(Fungi_norm, "phylum", 1:18)

#### Relative abundance ####
Fungi_phylum_reduced<-reduce_taxa(Fungi_phylum,3)
Fungi_phylum_abundance_by<-merge(Fungi_phylum_reduced, my_colData, by="row.names", all.x=T)
Fungi_phylum_abundance_byCond<-data.frame(lapply(Fungi_phylum_abundance_by[,1:length(colnames(Fungi_phylum_reduced))+1], function(x) tapply(x, Fungi_phylum_abundance_by$condition, mean, na.rm=TRUE)))
Fungi_phylum_abundance_byCond<-data.frame(melt.data.frame(Fungi_phylum_abundance_byCond/rowSums(Fungi_phylum_abundance_byCond)), value.names=rownames(Fungi_phylum_abundance_byCond))
Fungi_phylum_abundance_byCond$variable = with(Fungi_phylum_abundance_byCond, factor(variable, levels = rev(levels(variable))))
Fungi_phylum_abundance_byCond$value.names = factor(Fungi_phylum_abundance_byCond$value.names, levels=factor(c("CAL", "BAL", "CRO", "BRO","CCONT", "BCONT")))

f_phylum_plot<-ggplot(Fungi_phylum_abundance_byCond,aes(value.names,value,fill=variable))+
  geom_bar(stat="identity")+
  labs(title="Fungi")+
  scale_fill_manual(values=c("#B4B4B4B3","#4C005C","#1464C8B3","#7fc97f","#ef3b2c","#fdb462","#a6cee3", "#984ea3","#662506", "#a6cee3"),labels=rev(str_replace(unique(Fungi_phylum_abundance_byCond$variable), "root", "unclassified")))+
  theme_bw()+ylab("Reads Mean Relative Abundance")+xlab("Sample")+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(nrow=2,byrow=FALSE))+
  theme(panel.spacing = unit(0, "lines"))

save(f_phylum_plot, file = "f_phylum_plot.RData")

#### Cut LCA table at family level ####
Fungi_family<-LCA_level(Fungi_norm, "family", 1:18)

#### PERMANOVA ####
dis <- vegdist(t(Fungi_family), method="bray")

beta_disp<-with(my_colData, betadisper(dis,condition)) 
beta_disp
plot(beta_disp)### Plot the groups and distances to centroids on the first two PCoA axes
boxplot(beta_disp)
permutest(beta_disp)
TukeyHSD(beta_disp)
set.seed(1)
adonis<-adonis(t(Fungi_class) ~ genotype*soil, data=my_colData, permutations=9999, method="bray")
adonis$aov.tab$expl_variance<-(adonis$aov.tab$SumsOfSqs/adonis$aov.tab[dim(adonis$aov.tab)[1],2])*100
adonis

#### DEseq2 differential abundance ####
sig=0.05
load("data_host_filtered")
Fungi_raw<-filter(data_host_filtered, kingdom=="Fungi")

#### Phylum level ####
Fungi_phylum_raw<-LCA_level(Fungi_raw, "phylum", 1:18)
Fungi_phylum_raw<-subset(Fungi_phylum_raw, rownames(Fungi_phylum_raw)!="root")

### AL vs RO
Fungi_phylum_dds_AL_RO<-DESeqDataSetFromMatrix(countData= Fungi_phylum_raw, colData= my_colData[1:3], design=~soil)
Fungi_phylum_dds_AL_RO<-DESeq(Fungi_phylum_dds_AL_RO, fitType="parametric", betaPrior=T)
DESeq2_f_phylum_AL_RO<-results(Fungi_phylum_dds_AL_RO, contrast=c("soil", "AL","RO"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_f_phylum_AL_RO_filtered = subset(DESeq2_f_phylum_AL_RO, padj < sig)

### C vs B
Fungi_phylum_dds_C_B<-DESeqDataSetFromMatrix(countData= Fungi_phylum_raw, colData= my_colData[1:3], design=~genotype)
Fungi_phylum_dds_C_B<-DESeq(Fungi_phylum_dds_C_B, fitType="parametric", betaPrior=T)
DESeq2_f_phylum_C_B<-results(Fungi_phylum_dds_C_B, contrast=c("genotype", "C","B"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_f_phylum_C_B_filtered = subset(DESeq2_f_phylum_C_B, padj < sig)

#### Native vs Control 
Fungi_phylum_dds_N_C<-DESeqDataSetFromMatrix(countData= Fungi_phylum_raw, colData= my_colData[1:4], design=~cond)
Fungi_phylum_dds_N_C<-DESeq(Fungi_phylum_dds_N_C, fitType="parametric", betaPrior=T)
DESeq2_f_phylum_N_C<-results(Fungi_phylum_dds_N_C, contrast=c("cond", "N","CONT"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_f_phylum_N_C_filtered = subset(DESeq2_f_phylum_N_C, padj < sig & abs(log2FoldChange)>fc)

#### Family level ####
Fungi_family_raw<-LCA_level(Fungi_raw, "family", 1:18)
Fungi_family_raw<-subset(Fungi_family_raw, rownames(Fungi_family_raw)!="root")
              
### AL vs RO
Fungi_family_dds_AL_RO<-DESeqDataSetFromMatrix(countData= Fungi_family_raw, colData= my_colData[1:3], design=~soil)
Fungi_family_dds_AL_RO<-DESeq(Fungi_family_dds_AL_RO, fitType="parametric", betaPrior=T)

DESeq2_f_family_AL_RO<-results(Fungi_family_dds_AL_RO, contrast=c("soil", "AL","RO"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_f_family_AL_RO_filtered = subset(DESeq2_f_family_AL_RO, padj < sig)

DESeq2_f_family_AL_RO$Significant <- ifelse(rownames(DESeq2_f_family_AL_RO) %in% rownames(DESeq2_f_family_AL_RO_filtered) , "Yes", "No")

S_C_DET_fung <- ggplot(data = as.data.frame(DESeq2_f_family_AL_RO), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 0.4) +
  ggtitle("AL vs RO")+
  scale_x_log10() +
  geom_hline(yintercept = 0, color="firebrick", linetype=2)+
  scale_color_manual(values=c("gray20", pale_red), guide=FALSE) +
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(padj <0.05, rownames(DESeq2_f_family_AL_RO), ""), size=-log(padj,10)-0),show.legend = FALSE, segment.size = 0, force=0.2, max.iter = 100000)+
  scale_size(range=c(0,2))

S_C_DET_fung<-S_C_DET_fung + geom_point(data = as.data.frame(DESeq2_f_phylum_AL_RO_filtered), aes(x = baseMean, y = log2FoldChange), shape=16, color=light_blue)+
                            geom_text_repel(data = as.data.frame(DESeq2_f_phylum_AL_RO_filtered), aes(x = baseMean, y = log2FoldChange, label = rownames(DESeq2_f_phylum_AL_RO_filtered), size=-log(padj,10)), colour = light_blue, show.legend = FALSE, segment.size = 0, force=0.2, max.iter = 100000)
  
### C vs B
Fungi_family_dds_C_B<-DESeqDataSetFromMatrix(countData= Fungi_family_raw, colData= my_colData[1:3], design=~genotype)
Fungi_family_dds_C_B<-DESeq(Fungi_family_dds_C_B, fitType="parametric", betaPrior=T)

DESeq2_f_family_C_B<-results(Fungi_family_dds_C_B, contrast=c("genotype", "C","B"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_f_family_C_B_filtered = subset(DESeq2_f_family_C_B, padj < sig)

DESeq2_f_family_C_B$Significant <- ifelse(rownames(DESeq2_f_family_C_B) %in% rownames(DESeq2_f_family_C_B_filtered) , "Yes", "No")

S_R_DET_fung <- ggplot(data = as.data.frame(DESeq2_f_family_C_B), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 0.4) +
  ggtitle("'Cuore di Bue' vs 'Battito'")+
  scale_x_log10() +
  geom_hline(yintercept = 0, color="firebrick", linetype=2)+
  scale_color_manual(values=c("gray20", pale_red), guide=FALSE) +
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(padj <0.05, rownames(DESeq2_f_family_C_B), ""), size=-log(padj,10)-0), show.legend = FALSE, segment.size = 0, force=0.5, max.iter = 500000)+
  scale_size(range=c(1,2.3))

S_R_DET_fung<-S_R_DET_fung + geom_point(data = as.data.frame(DESeq2_f_phylum_C_B_filtered), aes(x = baseMean, y = log2FoldChange), shape=16, color=light_blue)+
  geom_text_repel(data = as.data.frame(DESeq2_f_phylum_C_B_filtered), aes(x = baseMean, y = log2FoldChange, label = rownames(DESeq2_f_phylum_C_B_filtered), size=-log(padj,10)), colour = light_blue, show.legend = FALSE, segment.size = 0, force=1.2, max.iter = 500000)

#### native vs Control 
fc=abs(2.5)
sig=0.05
Fungi_family_dds_N_C<-DESeqDataSetFromMatrix(countData= Fungi_family_raw, colData= my_colData[1:4], design=~cond)
Fungi_family_dds_N_C<-DESeq(Fungi_family_dds_N_C, fitType="parametric", betaPrior=T)

DESeq2_f_family_N_C<-results(Fungi_family_dds_N_C, contrast=c("cond", "N","CONT"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_f_family_N_C_filtered = subset(DESeq2_f_family_N_C, padj < 0.05 & abs(log2FoldChange)>abs(2.5))

plot(DESeq2_f_family_N_C$baseMean, DESeq2_f_family_N_C$log2FoldChange, pch=19, log="x", main="AL vs RO\nFungi (family)", xlab="Mean expression", ylab=expression(paste(log[2], "(fold-change)")),ylim=c(min(DESeq2_f_family_N_C$log2FoldChange)-0.5,max(DESeq2_f_family_N_C$log2FoldChange)+0.5),col=ifelse(DESeq2_f_family_N_C$padj<sig, "firebrick", "gray20"), cex=0.8)
text(DESeq2_f_family_N_C$baseMean+5, DESeq2_f_family_N_C$log2FoldChange+0.2, ifelse(DESeq2_f_family_N_C$padj<sig, rownames(DESeq2_f_family_N_C), ""), cex=0.8, col="firebrick")

DESeq2_f_family_N_C$Significant <- ifelse(rownames(DESeq2_f_family_N_C) %in% rownames(DESeq2_f_family_N_C_filtered) , "Yes", "No")

N_C_DET_fung <- ggplot(data = as.data.frame(DESeq2_f_family_N_C), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 0.4) +
  ggtitle("Native vs Control")+
  scale_x_log10() +
  geom_hline(yintercept = 0, color="firebrick", linetype=2)+
  scale_color_manual(values=c("gray20", pale_red), guide=FALSE) +
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(Significant=="Yes", rownames(DESeq2_f_family_N_C), ""), size=-log(padj,10)-0),show.legend = FALSE, segment.size = 0, force=0.5, max.iter = 500000)+
  scale_size(range=c(1,2.3))

N_C_DET_fung<-N_C_DET_fung + geom_point(data = as.data.frame(DESeq2_f_phylum_N_C_filtered), aes(x = baseMean, y = log2FoldChange), shape=16, color=light_blue)+
  geom_text_repel(data = as.data.frame(DESeq2_f_phylum_N_C_filtered), aes(x = baseMean, y = log2FoldChange, label = rownames(DESeq2_f_phylum_N_C_filtered), size=-log(padj,10)), colour = light_blue, show.legend = FALSE, segment.size = 0, force=1.2, max.iter = 500000, nudge_x = 0.9, )

#### Export Plots ####
save(S_C_DET_fung, S_R_DET_fung, N_C_DET_fung, DESeq2_f_family_AL_RO, DESeq2_f_family_N_C, DESeq2_f_family_C_B, DESeq2_f_phylum_AL_RO_filtered, DESeq2_f_phylum_C_B_filtered, DESeq2_f_phylum_N_C_filtered, file = "DEseq2_fung.RData")
