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
source("~/Plotting_styles.R")
# load custom functions
source_url("https://github.com/mchialva/myNGS_tools/blob/master/TaxR.R?raw=TRUE")

#### load Sample metadata and normalized reads counts ####
load("~my_colData")
load("~Bacteria_norm")

Bacteria_norm[Bacteria_norm==0]<-NA

#### Cut LCA table at phylum level ####
Bacteria_phylum<-LCA_level(Bacteria_norm, "phylum", 1:18)

#### Relative abundance ####
Bacteria_phylum_reduced<-reduce_taxa(Bacteria_phylum,3)
Bacteria_phylum_abundance_by<-merge(Bacteria_phylum_reduced, my_colData, by="row.names", all.x=T)
Bacteria_phylum_abundance_byCond<-data.frame(lapply(Bacteria_phylum_abundance_by[,1:length(colnames(Bacteria_phylum_reduced))+1], function(x) tapply(x, Bacteria_phylum_abundance_by$condition, mean, na.rm=TRUE)))
Bacteria_phylum_abundance_byCond<-data.frame(melt.data.frame(Bacteria_phylum_abundance_byCond/rowSums(Bacteria_phylum_abundance_byCond)), value.names=rownames(Bacteria_phylum_abundance_byCond))
Bacteria_phylum_abundance_byCond$variable = with(Bacteria_phylum_abundance_byCond, factor(variable, levels = rev(levels(variable))))
Bacteria_phylum_abundance_byCond$value.names = factor(Bacteria_phylum_abundance_byCond$value.names, levels=factor(c("CAL", "BAL", "CRO", "BRO","CCONT", "BCONT")))

b_phylum_plot<-ggplot(Bacteria_phylum_abundance_byCond,aes(value.names,value,fill=variable))+
  geom_bar(stat="identity")+
  labs(title="Bacteria")+
  scale_fill_manual(values=c("#B4B4B4B3","#ef3b2c","#fdb462","#a6cee3","darkorange", "forestgreen", "royalblue3", "firebrick", "salmon", "#FFCC99", "#0075DC","gray80","gray60","#0075DC","#4C005C"), labels=rev(str_replace(unique(Bacteria_phylum_abundance_byCond$variable), "root", "unclassified")))+
  theme_bw()+ylab("Reads Mean Relative Abundance")+xlab("Sample")+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(nrow=2,byrow=FALSE))

save(b_phylum_plot, file = "b_phylum_plot.RData")

#### Cut LCA table at family level ####
Bacteria_family<-LCA_level(Bacteria_norm, "family", 1:18)

#### PERMANOVA ####
dis <- vegdist(t(Bacteria_family), method="bray")

beta_disp<-with(my_colData, betadisper(dis,condition)) 
beta_disp
### Plot the groups and distances to centroids on the first two PCoA axes
plot(beta_disp)
boxplot(beta_disp)
permutest(beta_disp)
TukeyHSD(beta_disp)
set.seed(1)
adonis<-adonis(t(Bacteria_family) ~ genotype*soil, data=my_colData, permutations=9999, method="bray")
adonis$aov.tab$expl_variance<-(adonis$aov.tab$SumsOfSqs/adonis$aov.tab[dim(adonis$aov.tab)[1],2])*100
adonis

#### DEseq2 differential abundance ####
sig=0.05
load("data_host_filtered")
Bacteria_raw<-filter(data_host_filtered, superkingdom=="Bacteria")

#### Phylum level ####
Bacteria_phylum_raw<-LCA_level(Bacteria_raw, "phylum", 1:18)
Bacteria_phylum_raw<-subset(Bacteria_phylum_raw, rownames(Bacteria_phylum_raw)!="root")

### AL vs RO
Bacteria_phylum_dds_AL_RO<-DESeqDataSetFromMatrix(countData= Bacteria_phylum_raw, colData= my_colData[1:3], design=~soil)
Bacteria_phylum_dds_AL_RO<-DESeq(Bacteria_phylum_dds_AL_RO, fitType="parametric", betaPrior=T)
DESeq2_b_phylum_AL_RO<-results(Bacteria_phylum_dds_AL_RO, contrast=c("soil", "AL","RO"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_b_phylum_AL_RO_filtered = subset(DESeq2_f_phylum_AL_RO, padj < sig)

### C vs B
Bacteria_phylum_dds_C_B<-DESeqDataSetFromMatrix(countData= Bacteria_phylum_raw, colData= my_colData[1:3], design=~genotype)
Bacteria_phylum_dds_C_B<-DESeq(Bacteria_phylum_dds_C_B, fitType="parametric", betaPrior=T)
DESeq2_b_phylum_C_B<-results(Bacteria_phylum_dds_C_B, contrast=c("genotype", "C","B"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_b_phylum_C_B_filtered = subset(DESeq2_f_phylum_C_B, padj < sig)

#### Natural vs Control 
Bacteria_phylum_dds_N_C<-DESeqDataSetFromMatrix(countData= Bacteria_phylum_raw, colData= my_colData[1:4], design=~cond)
Bacteria_phylum_dds_N_C<-DESeq(Bacteria_phylum_dds_N_C, fitType="parametric", betaPrior=T)
DESeq2_b_phylum_N_C<-results(Bacteria_phylum_dds_N_C, contrast=c("cond", "N","CONT"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_b_phylum_N_C_filtered = subset(DESeq2_f_phylum_N_C, padj < sig & abs(log2FoldChange)>fc)

#### Family level ####
Bacteria_family_raw<-LCA_level(Bacteria_raw, "family", 1:18)
Bacteria_family_raw<-subset(Bacteria_family_raw, rownames(Bacteria_family_raw)!="root")

### AL vs RO
Bacteria_family_dds_AL_RO<-DESeqDataSetFromMatrix(countData= Bacteria_family_raw, colData= my_colData[1:3], design=~soil)
Bacteria_family_dds_AL_RO<-DESeq(Bacteria_family_dds_AL_RO, fitType="parametric", betaPrior=T)

DESeq2_b_family_AL_RO<-results(Bacteria_family_dds_AL_RO, contrast=c("soil", "AL","RO"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_b_family_AL_RO_filtered = subset(DESeq2_b_family_AL_RO, padj < sig)

DESeq2_b_family_AL_RO$Significant <- ifelse(rownames(DESeq2_b_family_AL_RO) %in% rownames(DESeq2_b_family_AL_RO_filtered) , "Yes", "No")

S_C_DET_bact <- ggplot(data = as.data.frame(DESeq2_b_family_AL_RO), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 0.4) +
  ggtitle("AL  vs RO")+
  scale_x_log10() +
  geom_hline(yintercept = 0, color="firebrick", linetype=2)+
  scale_color_manual(values=c("gray20", pale_red), guide=FALSE) +
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(padj <0.05, rownames(DESeq2_b_family_AL_RO), ""), size=-log(padj,10)-0),show.legend = FALSE, segment.size = 0, force=0.2, max.iter = 100000)+
  scale_size(range=c(0,2))

S_C_DET_bact<-S_C_DET_bact + geom_point(data = as.data.frame(DESeq2_b_phylum_AL_RO_filtered), aes(x = baseMean, y = log2FoldChange), shape=16, color=light_blue)+
  geom_text_repel(data = as.data.frame(DESeq2_b_phylum_AL_RO_filtered), aes(x = baseMean, y = log2FoldChange, label = rownames(DESeq2_b_phylum_AL_RO_filtered), size=-log(padj,10)), colour = light_blue, show.legend = FALSE, segment.size = 0, force=0.2, max.iter = 100000)

### C vs B
Bacteria_family_dds_C_B<-DESeqDataSetFromMatrix(countData= Bacteria_family_raw, colData= my_colData[1:3], design=~genotype)
Bacteria_family_dds_C_B<-DESeq(Bacteria_family_dds_C_B, fitType="parametric", betaPrior=T)

DESeq2_b_family_C_B<-results(Bacteria_family_dds_C_B, contrast=c("genotype", "C","B"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_b_family_C_B_filtered = subset(DESeq2_b_family_C_B, padj < sig)

DESeq2_b_family_C_B$Significant <- ifelse(rownames(DESeq2_b_family_C_B) %in% rownames(DESeq2_b_family_C_B_filtered) , "Yes", "No")

S_R_DET_bact <- ggplot(data = as.data.frame(DESeq2_b_family_C_B), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 0.4) +
  ggtitle("'Cuore di Bue' vs 'Battito'")+
  scale_x_log10() +
  geom_hline(yintercept = 0, color="firebrick", linetype=2)+
  scale_color_manual(values=c("gray20", pale_red), guide=FALSE) +
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(padj <0.05, rownames(DESeq2_b_family_C_B), ""), size=-log(padj,10)-1),show.legend = FALSE, segment.size = 0, force=0.5, max.iter = 500000)+
  scale_size(range=c(1,2.3))

S_R_DET_bact<-S_R_DET_bact + geom_point(data = as.data.frame(DESeq2_b_phylum_C_B_filtered), aes(x = baseMean, y = log2FoldChange), shape=16, color=light_blue)+
  geom_text_repel(data = as.data.frame(DESeq2_b_phylum_C_B_filtered), aes(x = baseMean, y = log2FoldChange, label = rownames(DESeq2_b_phylum_C_B_filtered), size=-log(padj,10)), colour = light_blue, show.legend = FALSE, segment.size = 0, force=1.2, max.iter = 500000)

#### Native vs Control 
Bacteria_family_dds_N_C<-DESeqDataSetFromMatrix(countData= Bacteria_family_raw, colData= my_colData[1:4], design=~cond)
Bacteria_family_dds_N_C<-DESeq(Bacteria_family_dds_N_C, fitType="parametric", betaPrior=T)

DESeq2_b_family_N_C<-results(Bacteria_family_dds_N_C, contrast=c("cond", "N","CONT"), independentFiltering = T, cooksCutoff = F, alpha = sig)
DESeq2_b_family_N_C_filtered = subset(DESeq2_b_family_N_C, padj < 0.05 & abs(log2FoldChange)>abs(0))

DESeq2_b_family_N_C$Significant <- ifelse(rownames(DESeq2_b_family_N_C) %in% rownames(DESeq2_b_family_N_C_filtered) , "Yes", "No")

N_C_DET_bact <- ggplot(data = as.data.frame(DESeq2_b_family_N_C), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 0.4) +
  ggtitle("Native vs Control")+
  scale_x_log10() +
  geom_hline(yintercept = 0, color="firebrick", linetype=2)+
  scale_color_manual(values=c("gray20", pale_red), guide=FALSE) +
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(padj <0.05 & abs(log2FoldChange)>abs(0), rownames(DESeq2_b_family_N_C), ""), size=-log(padj,10)-1),show.legend = FALSE, segment.size = 0, force=0.5, max.iter = 500000)+
  scale_size(range=c(1,2.3))

N_C_DET_bact<-N_C_DET_bact + geom_point(data = as.data.frame(DESeq2_b_phylum_N_C_filtered), aes(x = baseMean, y = log2FoldChange), shape=16, color=light_blue)+
  geom_text_repel(data = as.data.frame(DESeq2_b_phylum_N_C_filtered), aes(x = baseMean, y = log2FoldChange, label = rownames(DESeq2_b_phylum_N_C_filtered), size=-log(padj,10)), colour = light_blue, show.legend = FALSE, segment.size = 0, force=1.2, max.iter = 500000)

#### Export Plots ####
save(S_C_DET_bact, S_R_DET_bact, N_C_DET_bact, DESeq2_b_family_AL_RO, DESeq2_b_family_N_C, DESeq2_b_family_C_B, DESeq2_b_phylum_AL_RO_filtered, DESeq2_b_phylum_C_B_filtered, DESeq2_b_phylum_N_C_filtered, file = "DEseq2_bact.RData")
