library(dplyr)
library(reshape)
library(ggplot2)
library(car)
library(grid)
library(gridExtra)
library(gridGraphics)
library(jpeg)
library(agricolae)
library(multcomp)

setwd("/home/user/")

# load plotting styles
source("Plotting_styles.R")
# load custom functions
source_url("https://github.com/mchialva/myNGS_tools/blob/master/TaxR.R?raw=TRUE")

#### load Sample metadata ####
load("my_colData")

AMF_colonization<-read.delim("AM_colonization.tsv", header=T, dec=",")
colnames(AMF_colonization)<-c("X", "F", "M", "a", "A")

AMF_colonization_melted<-melt(AMF_colonization)
AMF_colonization_melted$genotype<-with(AMF_colonization_melted, substr(X, 1,1))
AMF_colonization_melted$soil<-with(AMF_colonization_melted, substr(X, 2,3))

## test differences ##
with(AMF_colonization, tapply(`F`, X, shapiro.test))
with(AMF_colonization, tapply(`M`, X, shapiro.test))
with(AMF_colonization, tapply(`a`, X, shapiro.test))
with(AMF_colonization, tapply(`A`, X, shapiro.test))
leveneTest(`F`~X, data=AMF_colonization)
model<-lm(`F` ~ X, data = AMF_colonization)
anova(model)
MSerror <- deviance(model)/df.residual(model)

library(lawstat)
with(AMF_colonization, levene.test(`F`, X, location="mean", kruskal.test=T))

library(FSA)
letters<-NULL
params<-as.factor(colnames(AMF_colonization)[2:5])

for (i in 1:length(params))
{
  seq<-params[i]
  formula<-as.formula(paste(seq, "~X"))
x<-data.frame(dunnTest(formula, data=AMF_colonization, method="none", two.sided=T)$res)
colnames(x)[4]<-"p adj"
rownames(x)<- paste(substr(x$Comparison,1, str_locate(x$Comparison, pattern = "-")-2), "-", substr(x$Comparison, str_locate(x$Comparison, pattern = "-")+2, str_length(x$Comparison)), sep="")
x<-x[,-1]
x$`p adj`<-ifelse(x$`p adj`>0.05, FALSE, TRUE)
xletters<-data.frame(multcompLetters(extract_p(x), reversed = T)$Letters)
xletters$trt<-rownames(xletters)
xletters$variable<-seq
colnames(xletters)[1]<-"X"
xletters$soil<-substr(rownames(xletters), 2,3)
xletters$genotype<-substr(rownames(xletters), 1,1)
xletters<-merge(xletters, melt(ddply(AMF_colonization[,c(1,i+1)],.(X), colwise(mean))), by.x="row.names", by.y="X")
letters<-rbind(xletters, letters)
}

colnames(letters)[4]<-"variable"

colonization<-ggplot(AMF_colonization_melted, aes(x = soil, y = value, fill=genotype)) +
  geom_bar(position=position_dodge(), stat="summary", fun.y = "mean", lwd=0.2) +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="gray10", width=0.2, position=position_dodge(0.9), na.rm = T)+
  scale_fill_manual(values = c(grey, sea_green), labels=c("'Battito'", "'Cuore di Bue'"))+
  scale_y_continuous(limits=c(0, 105))+
  geom_text(data = letters, aes(x=soil, y=(sqrt(value)+2)*5, label=X), position=position_dodge(width = 0.9), size = 2)+
  facet_grid(.~ variable) +
  theme_bw() +
  labs(x = "Soils", y = "%")+
  theme(legend.position="right")

load("Fungi_norm")
Fungi_phylum<-LCA_level(Fungi_norm, "phylum", 1:18)
Glomeromycota<-t(Fungi_phylum[row.names(Fungi_phylum)!="root",])

Glomeromycota_by<-merge(Glomeromycota, my_colData, by="row.names", all.x=T)
Glomeromycota_by$soil = factor(Glomeromycota_by$soil, levels=factor(c("CONT", "AL", "RO")))

with(Glomeromycota_by, tapply(Glomeromycota, condition, shapiro.test))
leveneTest(Glomeromycota~condition, data=Glomeromycota_by)
model<-lm(Glomeromycota ~ condition, data = Glomeromycota_by)
anova(model)
MSerror <- deviance(model)/df.residual(model)
comparison_fung <- HSD.test(Glomeromycota_by$Glomeromycota, Glomeromycota_by$condition, df.residual(model), MSerror, group=T, alpha=0.05)

Glomeromycota_boxplot<-ggplot(Glomeromycota_by, aes(x = soil, y = Glomeromycota, fill=genotype)) +
  geom_boxplot(lwd=0.2, outlier.size=0.5) +
  scale_fill_manual(values = c(grey, sea_green), labels=c("'Battito'", "'Cuore di Bue'"))+
  theme_bw()+
  labs(x = "Soils", y = "normalized reads")

#### load micrographs ####
images <- list.files(pattern=".jpg", all.files=T, full.names=T)
list_of_images = lapply(images, function(x) readJPEG(x, native = T))

#### Export Plots ####
save(colonization, Glomeromycota_boxplot, list_of_images, file = "AMF_plot.RData")
