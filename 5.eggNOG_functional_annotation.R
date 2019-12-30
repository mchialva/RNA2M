################################################
##### Anlyse TAXONER64 transcripts outputs #####
################################################
library(parallel)
library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(DESeq2)
library(gplots)
library(clusterProfiler)
#library(rentrez)
#library(KEGGREST)
#library(CHNOSZ)

# use all cores (threads) available
cores<-detectCores()

setwd("/home/user/")

# load custom functions
source_url("https://github.com/mchialva/myNGS_tools/blob/master/TaxR.R?raw=TRUE")
# load plotting styles
source("~/Plotting_styles.R")

#### load Sample metadata ####
load("~my_colData")
# load previously-generated statistics
load("~stats")
# load taxa table
load("~data_host_filtered")

#### Import TAXONER64 outputs (Taxonomy.txt files) ####
taxoner64_output<-"/home/user/5.Taxoner64/"
taxoner64_outputs<-list.files(taxoner64_output, pattern = "Taxonomy.txt$", full.names = T, recursive=T)

# retrive output files and put them in a list
outputs_samples=lapply(taxoner64_outputs, function(x) read.delim(x, header=F))
names(outputs_samples) <- rownames(my_colData)

#### Merge with previously analyzed taxonomy ####
data_host_filtered<-filter(data_host_filtered, rank.annotation!="Eukaryota"&rank.annotation!="root"&rank.annotation!="cellular organisms")
outputs_samples_filtered<-lapply(outputs_samples, function(x) merge(x, data_host_filtered[,19:23], by.x="V2", by.y="taxid"))

#### Select gi to retrieve and add coordinates of pseudoreads ####
reads_metadata<-data.frame(do.call(rbind, lapply(outputs_samples_filtered, function(x) x[,1:length(x)])))
reads_metadata$start<-ifelse(ifelse(reads_metadata$V5!=1, reads_metadata$V5-50, 1)<0, 1, ifelse(reads_metadata$V5!=1, reads_metadata$V5-50, 1))
reads_metadata$stop<-reads_metadata$start+150
write(paste(">", unique(reads_metadata$V3), sep=""), file = "MT_GIs.txt", sep="\n")

## select from fasta ##
# cat nt.fasta volumes from taxoner DB
taxoner_DB<-"/home/mchialva/DB/taxoner/bowtie2_update/"
system(paste("cat ", taxoner_DB, "*.fasta > nt_gi_taxid.fasta", sep=""))

# linearize nt_gi_taxid.fasta
cat(paste("awk '/^>/ {printf(\"\n%s\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\n\");}'  nt_gi_taxid.fasta > nt_gi_taxid_mono.fasta"))

# grep gi in nt database
system(paste("grep --no-group-separator -w -A1 -F -f ", paste(getwd(), "/MT_GIs.txt", sep=""), " nt_gi_taxid_mono.fasta > ", paste(getwd(), "/MT_GIs.fasta", sep=""), sep=""))

# cut fasta header to remove taxonomy
system("cut -d\";\" -f 1 MT_GIs.fasta > MT_GIs_cut.fasta")

# retrive GIs sequence length to better refine pseudoreads coordinates
system(cat("awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' MT_GIs_cut.fasta> MT_GIs_fasta_length.fa | sed ':a;N;$!ba;s/;\n/ /g' | sed 's/;/\t/g' | awk '{print $1\"\t\"$3}'| sed 's/>//g'"))

# add ";" at the end of header line to further parse length in table
system("sed --in-place '/>/s/.*/&;/' MT_GIs_fasta_length.txt")

#remove space at the end of header line (;)
system("sed --in-place ':a;N;$!ba;s/;\n/ /g' MT_GIs_fasta_length.txt")

# remove fasta header symbol ">"
system("sed --in-place 's/>//g' MT_GIs_fasta_length.txt")

# Import GIs sequence length
GIs_length <- read.delim("MT_GIs_fasta_length.txt", sep=" ", header=F)
colnames(GIs_length)<- c("gi", "length")
# Generate pseudoreads coordinates on reference and write .bed file
reads_metadata<-merge(reads_metadata, GIs_length, by.x="V3", by.y="gi")
reads_metadata$stop_adj<-ifelse(reads_metadata$stop>reads_metadata$length, reads_metadata$length, reads_metadata$stop)
reads_metadata$pseudoreads_length<-reads_metadata$stop_adj-reads_metadata$start
reads_metadata$mancanti<-reads_metadata$pseudoreads_length-150
reads_metadata$start_adj<-ifelse(reads_metadata$start+reads_metadata$mancanti>0,reads_metadata$start+reads_metadata$mancanti, 1)
reads_metadata$pseudoreads_final_length<-reads_metadata$stop_adj-reads_metadata$start_adj
# filter for pseudoreads length > 100 bp
reads_metadata<-filter(reads_metadata, pseudoreads_final_length>100)
# generate .bed file
reads_metadata_bed<-as.data.frame(reads_metadata[,c(1,17,15,3)])
write.table(reads_metadata_bed, file = "MT_coord.bed", sep="\t", row.names = F, col.names=F, quote=F)

# extact pseudoreads coordinates using .bed file.
system("bedtools getfasta -fi MT_GIs_cut.fasta -bed MT_coord.bed -name -fo MT_pseudoreads_150.fasta")
# remove coordinates from .fasta header
system("sed --in-place 's/::/;/g' MT_pseudoreads_150.fasta")
system("cut -d\";\" -f 1 MT_pseudoreads_150.fasta > MT_pseudoreads_150_cut.fasta")

#### BLASTX with DIAMOND on eggNOG database ####
# index eggNOG database for DIAMOND
#diamond makedb --in eggnog4.proteins.all.fa -d eggNOG -p 20 -c 1
system("diamond blastx -d /home/user/MT_pseudoreads_150_cut.fasta -k1 -e0.001 -p40 -b15 -c1 -o /home/user/MT_diamond_pseudoreads")
# import DIAMOND results in .tab format.
DIAMOND<- read.delim("MT_diamond_pseudoreads", header=F)
#Remove HSPs (high-scoring segment pairs), V9 & V10
DIAMOND<-unique(DIAMOND[, c(1, 2)])
# eggNOG database annotation files
eggNOG_wd<-"/home/user/eggNOG/"
NOG.annotations<- read.delim(paste(eggNOG_wd, "NOG.annotations.tsv", sep=""), header=F)
COG_functional_categories<-read.delim(paste(eggNOG_wd, "COG_functional_categories.tsv", sep=""), header=F)

# set-up eggNOG database
# download NOG.members.tsv and NOG.annotation.tsv
# merge eggNOG IDs retrieved from DIAMOND blastx with NOG.members running this multicore awk script
awk_cmd<-paste("\nsplit -n l/", cores, " ", file.path(getwd(), "MT_diamond_pseudoreads"),
              " --additional-suffix='.reads'\n",
              "for r in *.reads;\ndo awk 'NR==FNR{r[$2]=1;next} {n=0;for(i in a){if($6~i){print i, $1,$2,$5}}} n' r", file.path(getwd(), "NOG.members.tsv")," > ${r}.out & done\n",
              "wait\n", "cat *.out > outputfile", sep="")
# write bash script
cat(file = file.path(getwd(), "awk_parallel.sh"), awk_cmd, append = F)

# import eggNOG2COG and remove duplicates
eggNOG2COG<-unique(read.delim("/home/user/outputfile", header=F, sep=" "))

# annotate reads IDs
reads2COG<-merge(DIAMOND, eggNOG2COG, by.x="V2", by.y="V1", all.x=T)
# merge with libraries
outputs_samples_filtered_annotated<-mclapply(outputs_samples_filtered, function(x) merge(x, reads2COG, by="V1"), mc.cores=cores, mc.preschedule=T, mc.cleanup=T)

###### Summarize reads counts per COG ######
# generate reads count for each sample
COG_counts<-mclapply(outputs_samples_filtered_annotated, function(x) ddply(x,"V3.y",summarize,Freq=length(V1)), mc.cores=cores, mc.preschedule=T, mc.cleanup=T)
# merge all list object (libraries) in a single data.frame
COG_counts_df<- Reduce(function(...) merge(..., by="V3.y", all=TRUE), COG_counts)
COG_counts_df[is.na(COG_counts_df)]<-0
COG_counts_df$V3.y<-as.character(COG_counts_df$V3.y)
COG_counts_df[length(COG_counts_df$V3.y),1]<-"unassigned"
rownames(COG_counts_df)<-COG_counts_df$V3.y
COG_counts_df<-COG_counts_df[,-1]
colnames(COG_counts_df)<- rownames(my_colData)

### Reads count Filtering ###
### Remove functions with less than 5 reads in at least 5 samples ###
COG_counts_df_filtered<-filter_reads(COG_counts_df, 3, 5)
rownames(COG_counts_df_filtered)<-COG_counts_df_filtered$taxid
COG_counts_df_filtered<-COG_counts_df_filtered[,-19]
COG_counts_df_filtered<-COG_counts_df_filtered[-dim(COG_counts_df_filtered)[1],]

#### add summarized statistics to stats ####
stats$Sample<-substr(rownames(stats), 1, str_locate(rownames(stats), "R")+1)

stats$eggNOG_mapped<-colSums(COG_counts_df, na.rm = TRUE)
stats$eggNOG_unassigned<-t(filter(counts_df_fc, COG_category=="unassigned")[2:19])
stats$eggNOG_assigned<-stats$eggNOG_mapped-stats$eggNOG_unassigned

final_stats<-stats[,c(11,1,4,6,7)]
final_stats$percent.mapped<-round(final_stats$percent.mapped, 2)
final_stats$percent.eggNOG.mapped<-round(stats$eggNOG_mapped/(stats$assigned)*100, 2)
write.table(final_stats, "final_stats")

#### NMDS on COGs ####
#### NMDS - Soil type ####
set.seed(42)
function_ord <- metaMDS(t(COG_counts_df_filtered), distance = "bray", trymax=1000)
stressplot(function_ord)
plot(function_ord, type = "t")

NMDS = data.frame(MDS1 = function_ord$points[,1], MDS2 = function_ord$points[,2], group= my_colData$soil)
NMDS.mean=aggregate(NMDS[,1:2],list(group=my_colData$soil), mean)
NMDS$genotype=as.factor(substr(rownames(NMDS), 1, 1))
# Generate ellipses
ellipses<-ordiellipse(function_ord, my_colData$soil, display = "sites", kind = "se", conf = 0.95, label = T)
df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ellipses[[g]]$cov,ellipses[[g]]$center,ellipses[[g]]$scale))),group=g))}
# Plot NMDS
Soil_NMDS<-ggplot(NMDS) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=0.3, linetype=1)+
  geom_point(aes(x=MDS1, y=MDS2, color=group, shape=genotype), size=1)+
  scale_colour_manual(values=c(very_dark_green, dark_brown, red))+
  scale_shape_manual(values=c(16,17), labels=c("'Battito'", "'Cuore di Bue'")) +
  guides(colour=guide_legend(nrow=3,byrow=T, override.aes=list(size=1.5, shape=15)), shape=guide_legend(nrow=2,override.aes=list(size=1.5)))+ 
  main_theme

saveRDS(Soil_NMDS, file = "./Soil_NMDS.RDS")

# By genotypes
ellipses_genotype<-ordiellipse(function_ord, my_colData$genotype, display = "sites", kind = "se", conf = 0.95, label = T)
df_ell_genotype <- data.frame()
for(g in levels(NMDS$genotype)){
  df_ell_genotype <- rbind(df_ell_genotype, cbind(as.data.frame(with(NMDS[NMDS$genotype==g,],
                                                                     veganCovEllipse(ellipses_genotype[[g]]$cov,ellipses_genotype[[g]]$center,ellipses_genotype[[g]]$scale))),group=g))}
Soil_NMDS_genotype<-ggplot(NMDS) +
  geom_path(data=df_ell_genotype, aes(x=NMDS1, y=NMDS2,colour=group), size=0.3, linetype=1)+
  geom_point(aes(x=MDS1, y=MDS2, color=genotype, shape=group), size=1)+
  scale_colour_manual(values=c(pale_blue, grey, red), labels=c("'Battito'", "'Cuore di Bue'"))+
  scale_shape_manual(values=c(16,22,17), labels=c("AL", "CONT", "RO")) +
  guides(colour=guide_legend(nrow=3,byrow=T, override.aes=list(size=1.5, shape=15)), shape=guide_legend(nrow=2,override.aes=list(size=1.5)))+ 
  main_theme

saveRDS(Soil_NMDS_genotype, file = "./Soil_NMDS_genotype.RDS")

#### PERMANOVA ####
#COG genes
dis <- vegdist(t(COG_counts_df_filtered), method="bray")
beta_disp<-with(my_colData, betadisper(dis,condition)) 
beta_disp
### Plot the groups and distances to centroids on the first two PCoA axes
plot(beta_disp, hull = FALSE, ellipse = TRUE)
boxplot(beta_disp)
permutest(beta_disp)
TukeyHSD(beta_disp)
set.seed(1)
adonis<-adonis(t(COG_counts_df_filtered) ~ genotype*soil, data=my_colData, permutations=9999, method="bray")
adonis$aov.tab$expl_variance<-(adonis$aov.tab$SumsOfSqs/adonis$aov.tab[dim(adonis$aov.tab)[1],2])*100
adonis

sig<-0.05
###### Differentially expressed functions (DEFs) ######
###### DEFs on COGs ######
### N vs CONT
COG_dds_N_CONT<-DESeqDataSetFromMatrix(countData= COG_counts_df_filtered, colData= my_colData[1:4], design=~cond)
COG_dds_N_CONT<-DESeq(COG_dds_N_CONT, fitType="parametric", betaPrior=T)

COG_N_CONT<-results(COG_dds_N_CONT, contrast=c("cond", "N","CONT"), independentFiltering = T, cooksCutoff = F, alpha = sig)
COG_N_CONT_filtered = subset(COG_N_CONT, padj < sig)

COG_N_CONT<-merge(as.data.frame(COG_N_CONT), NOG.annotations, by.x="row.names", by.y="V2", all.x=T)

COG_N_CONT$Significant <- ifelse(COG_N_CONT$Row.names %in% rownames(COG_N_CONT_filtered) , "Yes", "No")

COG_N_C_DEF <- ggplot(data = as.data.frame(COG_N_CONT), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 1) +
  ggtitle("Native vs Control")+
  scale_x_log10() +
  scale_color_manual(values=c("black", "firebrick")) +
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(padj <0.01, as.character(COG_N_CONT$V6), ""), size=-log(padj,10)-1),show.legend = FALSE, segment.size = 0, force=0.4, max.iter = 10000, nudge_y = -0.4)

### C vs B
COG_dds_C_B<-DESeqDataSetFromMatrix(countData= COG_counts_df_filtered, colData= my_colData[1:4], design=~genotype)
COG_dds_C_B<-DESeq(COG_dds_C_B, fitType="parametric", betaPrior=T)

COG_C_B<-results(COG_dds_C_B, contrast=c("genotype", "C","B"), independentFiltering = T, cooksCutoff = F, alpha = sig)
COG_C_B_filtered = subset(COG_C_B, padj < sig)

COG_C_B<-merge(as.data.frame(COG_C_B), NOG.annotations, by.x="row.names", by.y="V2", all.x=T)

COG_C_B$Significant <- ifelse(COG_C_B$Row.names %in% rownames(COG_C_B_filtered) , "Yes", "No")

COG_C_B_DEF <- ggplot(data = as.data.frame(COG_C_B), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 2) +
  ggtitle("'Cuore di Bue' vs 'Battito'")+
  scale_x_log10() +
  scale_color_manual(values=c("black", "firebrick"))+
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(padj <sig, as.character(COG_C_B$V6), ""), size=-log(padj,10)-1),show.legend = FALSE, segment.size = 0, force=0.4, max.iter = 10000, nudge_y = -0.4)

### AL vs RO
COG_dds_AL_RO<-DESeqDataSetFromMatrix(countData= COG_counts_df_filtered, colData= my_colData[1:4], design=~soil)
COG_dds_AL_RO<-DESeq(COG_dds_AL_RO, fitType="parametric", betaPrior=T)

COG_AL_RO<-results(COG_dds_AL_RO, contrast=c("soil", "AL","RO"), independentFiltering = T, cooksCutoff = F, alpha = sig)
COG_AL_RO_filtered = subset(COG_AL_RO, padj < sig)

COG_AL_RO<-merge(as.data.frame(COG_AL_RO), NOG.annotations, by.x="row.names", by.y="V2", all.x=T)
COG_AL_RO$Significant <- ifelse(COG_AL_RO$Row.names %in% rownames(COG_AL_RO_filtered) , "Yes", "No")

COG_AL_RO_DEF <- ggplot(data = as.data.frame(COG_AL_RO), aes(x = baseMean, y = log2FoldChange, color = Significant), shape=16) +
  geom_point(size = 2) +
  ggtitle("AL vs RO")+
  scale_x_log10() +
  scale_color_manual(values=c("black", "firebrick")) +
  labs(x = "Mean abundance", y = expression(paste(log[2], "fold-change")))+theme_bw()+
  geom_text_repel(aes(x = baseMean, y = log2FoldChange, label = ifelse(padj <sig, as.character(COG_AL_RO$V6), ""), size=-log(padj,10)-1),show.legend = FALSE, segment.size = 0, force=0.4, max.iter = 10000, nudge_y = -0.4)

####### DE-COGs boxplots ########
# C versus B
COG_C_B_annot<-merge(COG_C_B, COG_functional_categories, by.x="V5", by.y="V1", all.x=T, all.y=F)
COG_C_B_annot$category<-paste(COG_C_B_annot$V2, " [", COG_C_B_annot$V5, "]", sep="")

C_B_boxplot_COG_by_cat<-ggplot(filter(COG_C_B_annot, V5%in%c("S", "G", "O", "U", "J", "T", "E")), aes(x = category, y = log2FoldChange))+theme_classic()+
  scale_x_discrete()+
  geom_point(position=position_jitter(width = 0.3), alpha=0.7, aes(color=Significant), size=0.5)+
  scale_shape_manual(values=19)+coord_flip()+
  scale_color_manual(values=c("#aec05d", "firebrick"), guide=F)+
  geom_boxplot(fill=NA, outlier.shape=NA, lwd=0.2)+
  geom_hline(yintercept = 0, color="firebrick", linetype=2)+
  labs(x = "COG functions", y = expression(paste(log[2], "fold-change")))+
  ggtitle("Susceptible vs Resistant")

# N versus CONT
COG_N_CONT_annot<-merge(COG_N_CONT, COG_functional_categories, by.x="V5", by.y="V1", all.x=T, all.y=F)
COG_N_CONT_annot$category<-paste(COG_N_CONT_annot$V2, " [", COG_N_CONT_annot$V5, "]", sep="")

N_CONT_boxplot_COG_by_cat<-ggplot(filter(COG_N_CONT_annot, V5%in%c("J", "C", "S", "G")), aes(x = category, y = log2FoldChange))+theme_classic()+
  scale_x_discrete()+
  geom_point(position=position_jitter(width = 0.3), alpha=0.7, aes(color=Significant), size=0.5)+
  scale_shape_manual(values=19)+coord_flip()+
  scale_color_manual(values=c("#aec05d", "firebrick"), guide=F)+
  geom_boxplot(fill=NA, outlier.shape=NA, lwd=0.2)+
  geom_hline(yintercept = 0, color="firebrick", linetype=2)+
  labs(x = "COG functions", y = expression(paste(log[2], "fold-change")))+
  ggtitle("Native vs Control")

save(C_B_boxplot_COG_by_cat, COG_C_B_annot, N_CONT_boxplot_COG_by_cat, COG_N_CONT_annot, N_CONT_boxplot_COG, fc_COG, file="Function_boxplots.RData")

#### DESeq2 normalization and data export ####
## COG genes
COG_dds<-DESeqDataSetFromMatrix(countData= COG_counts_df_filtered, colData= my_colData[1:3], design=~condition)
COG_dds<-DESeq(COG_dds, fitType="parametric", betaPrior=TRUE)
COG_norm<-counts(COG_dds, normalized=TRUE)

## export data
save(COG_counts_df_filtered, file="COG_raw_counts_filtered")
write.table(COG_counts_df_filtered, file="COG_raw_counts_filtered.tsv", quote=F, sep = "\t", row.names = T, col.names = T)
save(COG_norm, file="COG_norm")
