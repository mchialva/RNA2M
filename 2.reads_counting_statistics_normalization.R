################################################
######## Anlyse TAXONER64 taxonomy outputs #####
################################################
# This script takes in input taxoner64 results and generates raw counts matrix with taxonomical annotation and statistics
library(CHNOSZ)
library(parallel)
library(ShortRead)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(reshape)
library(Biobase)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(agricolae)
library(car)
library(vegan)
library(devtools)

setwd("/home/user/")

# load custom functions
source_url("https://github.com/mchialva/myNGS_tools/blob/master/TaxR.R?raw=TRUE")
# load plotting styles
source("~/Plotting_styles.R")

# paths to sortmerna and taxoner outputs
taxoner64_output<-"/home/user/5.Taxoner64/"
Sortmerna_output<-"/home/user/4.Sortmerna/"

# Import local NCBI taxonomy
taxdir <- "/home/user/taxonomy/"
nodes <- getnodes(taxdir=taxdir)
names=getnames(taxdir=taxdir)

# use all cores (threads) available
cores<-detectCores()

#### Import TAXONER64 outputs (megan.txt files) ####
# generate taxoner64 outputs names
taxoner64_outputs<-list.files(taxoner64_output, pattern = "megan.txt$", full.names = T, recursive=T)
Sortmerna_unmapped<-list.files(Sortmerna_output, pattern="_rRNA_free.fastq$", full.names=T)

# retrive output files and generate a list
outputs_samples=lapply(taxoner64_outputs, function(x) read.delim(x, header=F))
names(outputs_samples) <- paste(str_replace(substr(taxoner64_outputs, 54, str_length(taxoner64_outputs)-72), "-", "_"), 1:length(outputs_samples), sep = "")

###### Stats ######
# Generate mapping statistics
mapped<-lapply(outputs_samples, function(x) length(x[[1]]))
filtered<-lapply(Sortmerna_unmapped, function(x) length(readFastq(x)))
names(filtered) <- paste(str_replace(substr(Sortmerna_unmapped, 54, str_length(Sortmerna_unmapped)-50), "-", "_"), 1:length(outputs_samples), sep = "")

stats<-data.frame(cbind(t(data.frame(filtered)),t(data.frame(mapped))))
stats$unassigned<-(stats$X1)-(stats$X2)
stats<-cbind(stats, (stats$X2/stats$X1)*100)
colnames(stats)<-c("filtered", "assigned", "unmapped", "percent.mapped")

# Read redundancy statistics
library(ShortRead)
uniquef <- occurrenceFilter(withSread=TRUE)
stats$filtered.unique<-unlist(lapply(Sortmerna_unmapped, function(x) length(readFastq(x)[uniquef(readFastq(x))])))
stats$redundant.reads<-(stats$filtered) - (stats$filtered.unique)
stats$percent.unique<-100-(stats$redundant.reads/stats$filtered*100)

### import mapping table ###
my_colData<-read.delim("my_colData.tsv", header = T, row.names = 1)

###### Summarize reads counts ######
# generate reads count for each library
counts<-mclapply(outputs_samples, function(x) ddply(x,"V2",summarize,Freq=length(V1)), mc.cores=cores, mc.preschedule=T, mc.cleanup=T)
# merge all list object (libraries) in a single data.frame
counts_df<- Reduce(function(...) merge(..., by="V2", all=TRUE), counts)
rownames(counts_df)<-counts_df$V2
counts_df<-counts_df[,-1]
colnames(counts_df)<- Samples_name

#### Reads count Filtering ####
#### Remove tax IDs with less than 5 reads in at least 5 samples ####
counts_df_filtered<-filter_reads(counts_df, 3, 5)

#### annotate DB with LCA ranks ####
data.annot<-data.frame(counts_df_filtered, mclapply(data.frame(counts_df_filtered$taxid), function(x) getrank(x, taxdir, nodes=nodes), mc.cores=cores, mc.preschedule=T, mc.cleanup=T))
colnames(data.annot)[length(data.annot)]<-"rank"
data.annot<-data.frame(data.annot, rank.annotation=mclapply(data.frame(data.annot$taxid), function(x) sciname(x, taxdir, names=names), mc.cores=cores, mc.preschedule=T, mc.cleanup=T))
colnames(data.annot)[length(data.annot)]<-"rank.annotation"

#### Assign complete taxonomy to each reads ####
#ranks<-unique(nodes$rank) ## -- ALL ranks
# set selected ranks to annotate
ncols<-length(data.annot)
ranks<-c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  for (i in 1:length(ranks))
  {
    seq<-as.character(ranks)[i]
    data.annot<-data.frame(data.annot, mcmapply(data.annot$taxid, FUN=function(x)
      sciname(parent(x, taxdir, rank=seq, nodes=nodes),
              taxdir, names=names),
      mc.cores=cores, mc.preschedule=T, mc.cleanup=T))
    colnames(data.annot)[i+ncols]<-seq
  }
save(data.annot, file="data.annot")
load("data.annot")

#### Remove residual host categories: kingdom=Viridiplantae & Metazoa ####
data_host_filtered<-filter(data.annot, kingdom!="Viridiplantae"&kingdom!="Metazoa")
colnames(data_host_filtered)[1:18]<-str_replace(colnames(data_host_filtered)[1:18], "AV", "RO")
#### Write complete filtered counts table & annotation ####
save(data_host_filtered, file="data_host_filtered")
write.table(data_host_filtered, file="data_host_filtered.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

#### Filter reads count at superkingdom rank and plot statistics ####
superkingdom<-t(data.frame(lapply(data_host_filtered[,1:18], function(x) tapply(x, data_host_filtered$superkingdom, sum, na.rm=TRUE))))
#### by condition
superkingdom_abundance_by<-merge(superkingdom, my_colData, by="row.names", all.x=T)
superkingdom_abundance_byCond<-data.frame(lapply(superkingdom_abundance_by[,2:4], function(x) tapply(x, superkingdom_abundance_by$condition, mean, na.rm=TRUE)))
superkingdom_abundance_byCond<-data.frame(melt(superkingdom_abundance_byCond/rowSums(superkingdom_abundance_byCond)), value.names=rownames(superkingdom_abundance_byCond))

superkingdom_abundance_byCond$variable = with(superkingdom_abundance_byCond, factor(variable, levels = rev(levels(variable))))
superkingdom_abundance_byCond$value.names = factor(superkingdom_abundance_byCond$value.names, levels=factor(c("CAL", "BAL", "CRO", "BRO","CCONT", "BCONT")))

superkingdom_plot_byCond<-ggplot(superkingdom_abundance_byCond,aes(value.names,value,fill=variable))+
  geom_bar(stat="identity")+
  #c("#B4B4B4B3","#7fc97f","#1464C8B3","royalblue4", "forestgreen","firebrick", "yellow", "purple")
  scale_fill_manual(values=rev(brewer.pal(3,"Set2")), labels=rev(str_replace(unique(superkingdom_abundance_byCond$variable), "root", "unclassified")))+
  theme_bw()+ylab("Reads mean Relative Abundance")+xlab("Sample")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x=element_text(size=4, angle=90,hjust=1,vjust=0.5),
        legend.title=element_blank(),
        legend.key = element_blank(),
        axis.line.x = element_line(size=0.2, color = "black"),
        axis.line.y = element_line(size=0.2, color = "black"),
        panel.grid=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 5, colour ="black"),
        axis.title.y = element_text(size = 6, colour ="black"),
        axis.text.y = element_text(size=4, angle=0),
        axis.ticks=element_line(size=0.2),
        legend.text=element_text(size=4),
        legend.position="bottom",
        panel.spacing = unit(x = c(0.2, 0, -0.7, 0), "line"),
        legend.key.size = unit(0.4, "cm"))

#### Export plot ####
saveRDS(superkingdom_plot_byCond, file = "./superkingdom.RDS")

#### DESeq2 normalization ####
library(DESeq2)
## All Bacteria and Fungi
Taxa<-data_host_filtered[1:18]
rownames(Taxa)<-data_host_filtered$taxid
Taxa[is.na(Taxa)]<-0

Taxa_dds<-DESeqDataSetFromMatrix(countData= Taxa, colData= my_colData[1:3], design=~condition)
Taxa_dds<-DESeq(Taxa_dds, fitType="parametric", betaPrior=TRUE)
Taxa_norm<-cbind(counts(Taxa_dds, normalized=TRUE), data_host_filtered[19:length(data_host_filtered)])

## Bacteria
Bacteria_raw<-filter(data_host_filtered, superkingdom=="Bacteria")

Bacteria_raw_matrix<-Bacteria_raw[1:18]
rownames(Bacteria_raw_matrix)<-Bacteria_raw$taxid
Bacteria_raw_matrix[is.na(Bacteria_raw_matrix)]<-0

Bacteria_dds<-DESeqDataSetFromMatrix(countData= Bacteria_raw_matrix, colData= my_colData[1:3], design=~condition)
dds <- estimateSizeFactors(Bacteria_dds)
dds <- estimateDispersions(dds)
head(dispersions(dds))
plotDispEsts(dds)
Bacteria_dds<-DESeq(Bacteria_dds, fitType="parametric", betaPrior=TRUE)
Bacteria_norm<-cbind(counts(Bacteria_dds, normalized=TRUE), Bacteria_raw[19:length(Bacteria_raw)])

## Fungi
Fungi_raw<-filter(data_host_filtered, kingdom=="Fungi")

Fungi_raw_matrix<-Fungi_raw[1:18]
rownames(Fungi_raw_matrix)<-Fungi_raw$taxid
Fungi_raw_matrix[is.na(Fungi_raw_matrix)]<-0

Fungi_dds<-DESeqDataSetFromMatrix(countData= Fungi_raw_matrix, colData= my_colData[1:3], design=~condition)
dds <- estimateSizeFactors(Fungi_dds)
dds <- estimateDispersions(dds)
head(dispersions(dds))
plotDispEsts(dds)
Fungi_dds<-DESeq(Fungi_dds, fitType="parametric", betaPrior=TRUE)
Fungi_norm<-cbind(counts(Fungi_dds, normalized=TRUE), Fungi_raw[19:length(Fungi_raw)])

### Add superkingdom information to stats ###
stats$Bacteria<-colSums(filter(data_host_filtered, superkingdom=="Bacteria")[,c(1:18)], na.rm = TRUE)
stats$Eukaryota<-colSums(filter(data_host_filtered, superkingdom=="Eukaryota")[,c(1:18)], na.rm = TRUE)
stats$host_residual<-colSums(filter(data.annot, kingdom=="Viridiplantae")[,c(1:18)], na.rm = TRUE)
stats$contaminants<-colSums(filter(data.annot, kingdom=="Metazoa")[,c(1:18)], na.rm = TRUE)
stats$unassigned<-colSums(filter(data.annot, superkingdom=="root")[,c(1:18)], na.rm = TRUE)
stats$low_expressed <-colSums(anti_join(counts_df,counts_df_filtered), na.rm = TRUE)

### Write stat file ###
save(stats, file="stats")
### Write Normalized values files ###
save(Taxa_norm, file="Taxa_norm")
save(Bacteria_norm, file="Bacteria_norm")
save(Fungi_norm, file="Fungi_norm")
### Write Sample metadata ###
save(my_colData, file="my_colData")