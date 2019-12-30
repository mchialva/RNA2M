library(vegan)
library(adespatial)
library(dplyr)
library(ggplot2)
library(devtools)

setwd("/home/user/")

# load custom functions
source_url("https://github.com/mchialva/myNGS_tools/blob/master/TaxR.R?raw=TRUE")
# load plotting styles
source("~/Plotting_styles.R")

# import mapping table
my_colData<-read.delim("my_colData.tsv", header = T, row.names = 1)

# load normalized counts
load("Taxa_norm")
load("COG_norm")

# load plant host transcriptome normalized counts (from Chialva et al. 2018)
host_norm_counts<- t(read.delim("host_normalized_counts.tsv", header=T, row.names = 1, sep="\t"))

#### VARPART on taxonomy as explanatory variable ####
# all taxa
# Taxa_norm
# only bacteria
taxaAll_bacteria<-filter(Taxa_norm, superkingdom=="Bacteria"&rank.annotation!="root")[,1:18]
# only fungi
taxaAll_fungi<-filter(Taxa_norm, superkingdom=="Eukaryota"&rank.annotation!="root")[,1:18]

###### all taxa at LCA levels ######
# perform Variance Partitioning analysis (all taxa)
Tax_varpart_allTax <- varpart(Y = host_norm_counts, X =~genotype, ~soil, t(Taxa_norm[,1:18]), data = my_colData, transfo ="hel")
# too much collinearity detected.
# perform forward selection on taxa
set.seed(123)
taxa_fwsel_allTax<-forward.sel(host_norm_counts, t(Taxa_norm[1:18]), alpha=0.05, nperm=999)#, adjR2thresh=fun_r2adj)
taxa_sel_allTax<-t(Taxa_norm[1:18])[,taxa_fwsel_allTax$order]
# perform Variance Partitioning using only selected variables
Tax_varpart_allTax <- varpart(Y = host_norm_counts, X =~genotype, ~soil, taxa_sel_allTax, data = my_colData, transfo ="hel")
plot(Tax_varpart_allTax)

# Test individual fractions
taxAll_genotype_rda_1<-rda(host_norm_counts ~genotype+Condition(soil)+Condition(taxa_sel_allTax), data=my_colData)
taxAll_genotype_aov_1<-anova(taxAll_genotype_rda_1, step=1000, perm.max=9999, parallel=10)

taxAll_soil_rda_1<-rda(host_norm_counts ~soil+Condition(genotype)+Condition(taxa_sel_allTax), data=my_colData)
taxAll_soil_aov_1<-anova(taxAll_soil_rda_1, step=1000, perm.max=9999, parallel=10)

taxAll_micro_rda_1<-rda(host_norm_counts ~taxa_sel_allTax+Condition(genotype)+Condition(soil), data=my_colData)
taxAll_micro_aov_1<-anova(taxAll_micro_rda_1, step=1000, perm.max=9999, parallel=10)

#### Fungi Bacteria at LCA levels plus Genotype and Soil factor  ####
# perform Variance Partitioning analysis combining forward-selected bacterial and fungal taxa
Tax_varpart_allTax_sep <- varpart(Y = host_norm_counts, X =~genotype, ~soil, taxa_sel_allTax_bact, taxa_sel_allTax_fung, data = my_colData, transfo ="hel")
plot(Tax_varpart_allTax_sep)

# Test individual fractions
taxAll_genotype_rda_sep<-rda(host_norm_counts ~genotype+Condition(soil)+Condition(taxa_sel_allTax_bact)+Condition(taxa_sel_allTax_fung), data=my_colData)
taxAll_genotype_aov_sep<-anova(taxAll_genotype_rda_sep, step=1000, perm.max=9999, parallel=10)

taxAll_micro_rda_sep<-rda(host_norm_counts ~taxa_sel_allTax_bact+Condition(genotype)+Condition(soil)+Condition(taxa_sel_allTax_fung), data=my_colData)
taxAll_micro_aov_sep<-anova(taxAll_micro_rda_sep, step=1000, perm.max=9999, parallel=10)

#### VARPART on MT COGs as explanatory variable ####
COG_varpart <- varpart(Y = host_norm_counts, X =~genotype, ~soil, t(COG_norm), data = my_colData, transfo ="hel")
COG_varpart
# too much collinearity detected
# perform forward selection on COGs
COG_fwsel<-forward.sel(host_norm_counts, t(COG_norm), alpha=0.05, nperm=999)#adjR2thresh=COG_r2adj,
COG_fwsel
COG_sel<-t(COG_norm)[,COG_fwsel$order]

## Repeat Varpart using selected COGs #
COG_varpart <- varpart(Y = host_norm_counts, X =~genotype, ~soil, COG_sel, data = my_colData, transfo ="hel")
plot(COG_varpart)

# Test Individual Fractions
COG_sel_rda<-rda(host_norm_counts ~COG_sel+Condition(genotype)+Condition(soil), my_colData)
COG_aov<-anova(COG_sel_rda, step=1000, perm.max=9999, parallel=10)

COG_soil_rda<-rda(host_norm_counts ~soil+Condition(genotype)+Condition(COG_sel), data=my_colData)
COG_soil_aov<-anova(COG_soil_rda, step=1000, perm.max=9999, parallel=10)

COG_genotype_rda<-rda(host_norm_counts ~genotype+Condition(soil)+Condition(COG_sel), data=my_colData)
COG_genotype_aov<-anova(COG_genotype_rda, step=1000, perm.max=9999, parallel=10)

#### Plot Fig. 5 ####
tiff(filename = "Fig6_updated.tiff", width = 11.4, height = 8, unit="cm", res=600, compression="lzw")
par(bty = 'n', mar=c(0,0.5,0.2,0.1), mfrow=c(1,2))
#plot(Tax_varpart, Xnames=NA, bg=c(sea_green, pale_blue, pale_red), cex=0.8) # old version
plot(Tax_varpart_allTax, Xnames=NA, bg=c(sea_green, pale_blue, pale_red), cex=0.7)
text("a", cex=1.2, x=-0.7, y=1.2)
text(x = c(-0.5, 1.5, 0.5), y=c(0.8,0.8,-1.75), c("genotype","soil","active microbiota diversity***"), font=3, cex=0.75)
plot(COG_varpart, Xnames=NA, bg=c(sea_green, pale_blue, pale_red), cex=0.7)
text("b", cex=1.2, x=-0.7, y=1.2)
text(x = c(-0.5, 1.5, 0.5), y=c(0.8,0.8,-1.75), c("genotype","soil","functional diversity***"), font=3, cex=0.75)
dev.off()

# Supplementary figure 2
tiff(filename = "FigS2.tiff", width = 8, height = 8, unit="cm", res=600, compression="lzw")
par(bty = 'n', mar=c(0,0.5,0.2,0.1))
plot(Tax_varpart_allTax_sep, Xnames=NA, cex=0.7)
text(x = c(-1.5, -0.9, 1, 1.5), y=c(0.7,1.2,1.2,0.7), c("genotype","soil","bacteria", "fungi"), font=3, cex=0.75)
dev.off()
