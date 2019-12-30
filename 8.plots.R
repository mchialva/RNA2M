library(ggplot2)
library(gridExtra)
library(gridGraphics)
library(cowplot)

setwd("/home/user/")
# load plotting styles
source("~/Plotting_styles.R")

# load generated graphs and data (previous scripts)
superkingdom<-readRDS("superkingdom.RDS")
load("b_phylum_plot.RData")
load("f_phylum_plot.RData")
load("DEseq2_bact.RData")
load("DEseq2_fung.RData")
load("AMF_plot.RData")
load("Function_boxplots.RData")
Soil_NMDS<-readRDS("Soil_NMDS.RDS")
Soil_NMDS_genotype<-readRDS("Soil_NMDS_genotype.RDS")

# load gglegend function
gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}

#### Fig 1 Diversity Overview ####
tiff(filename = "Fig1.tiff", width = 12, height = 15, unit="cm", res=600, compression="lzw")
grid.arrange(arrangeGrob(superkingdom+labs(title="")+barplot_theme+xlab("")+theme(legend.position='bottom', legend.key.size = unit(0.5, "cm"),legend.text=element_text(size=7), plot.margin=unit(c(0.1,0.5,0.1,0.05), "cm"))+coord_flip()+theme(axis.title.x = element_text(size = 7, vjust= -0.5, colour ="black"), axis.text.x = element_text(size=6, angle = 0, hjust = 0.5)),
                         plot_grid(b_phylum_plot+barplot_theme+theme(axis.text.x = element_text(size=7),legend.position='bottom', legend.key.size = unit(0.5, "cm"),legend.text=element_text(size=7)),
                                   f_phylum_plot+barplot_theme+no_y_theme+theme(axis.text.x = element_text(size=7), legend.position='bottom',legend.key.size = unit(0.5, "cm"),legend.text=element_text(size=7)),
                                   ncol=2, align="h"),
                         nrow=2, heights=c(0.9,2)))
grid.text(c("a", "b", "c"), c(0.02,0.02,0.52), c(0.99,0.68,0.68), gp=gpar(col="black", cex=1),
          just=c("left"))
dev.off()

#### Fig. 2 ####
tiff(filename = "Fig2_all.tiff", width = 11.6, height = 17.4, unit="cm", res=600, compression="lzw")
grid.arrange(arrangeGrob(S_R_DET_bact+MA_theme_large+geom_rug(col=rgb(.5,0,0,alpha=.2)), S_R_DET_fung+MA_theme_large+geom_rug(col=rgb(.5,0,0,alpha=.2)), ncol=2, widths=c(1,1)),
             arrangeGrob(N_C_DET_bact+MA_theme_large+geom_rug(col=rgb(.5,0,0,alpha=.2)), N_C_DET_fung+MA_theme_large+geom_rug(col=rgb(.5,0,0,alpha=.2)), ncol=2, widths=c(1,1)),
             arrangeGrob(S_C_DET_bact+MA_theme_large+geom_rug(col=rgb(.5,0,0,alpha=.2)), S_C_DET_fung+MA_theme_large+geom_rug(col=rgb(.5,0,0,alpha=.2)), ncol=2, widths=c(1,1)),
             nrow=3, heights=c(1,1,1))
grid.text(c("a", "b", "c", "d", "e", "f"), c(0.02,0.52,0.02,0.52,0.02,0.52), c(0.99,0.99,0.66, 0.66, 0.33, 0.33), gp=gpar(col="black", cex=1),
          just=c("left"))
dev.off()

#### Fig 3 ####
tiff(filename = "Fig3_amf.tiff", width = 18, height = 12, unit="cm", res=600, compression="lzw")
grid.arrange(arrangeGrob(colonization+boxplot_theme+theme(legend.position="none", plot.margin=unit(c(0.1,0.2,0.5,0.2), "cm")),
                         Glomeromycota_boxplot+boxplot_theme+theme(legend.position="none", plot.margin=unit(c(0.2,0.2,0.5,0.2), "cm")), widths=c(1.5,1), ncol=2),
             arrangeGrob(gglegend(colonization+theme(legend.direction="horizontal", legend.text = element_text(size = 7), legend.key.size = unit(5, "mm"), legend.margin=margin(0,0,0,0))+guides(fill=guide_legend(title="")))),
             arrangeGrob(rasterGrob(list_of_images[[1]], width = unit(4.9*1.5, "cm"), height = unit(1.974*1.5, "cm")), rasterGrob(list_of_images[[2]],  width = unit(4.9*1.5, "cm"), height = unit(1.974*1.5, "cm")), ncol=2),
             nrow=3, heights=c(1,0.05,0.5))
grid.text(c("a", "b", "c", "d"), c(0.01,0.62, 0.01, 0.5), c(0.98,0.98,0.28, 0.28), gp=gpar(col="black", cex=1),
          just=c("left"))
dev.off()

#### Fig 4 functions (COGs) ####
tiff(filename = "Fig4.tiff", width = 21, height = 12, unit="cm", res=600, compression="lzw")
grid.arrange(arrangeGrob(Soil_NMDS_genotype+guides(colour=guide_legend(nrow=1,byrow=T, override.aes=list(size=1.5, shape=15)), shape=guide_legend(nrow=1,override.aes=list(size=1.5)))+main_theme+theme(legend.text=element_text(size=6), legend.position=c(0,0.1), legend.spacing = unit(-0.5, "cm")),
                         Soil_NMDS+guides(colour=guide_legend(nrow=1,byrow=T, override.aes=list(size=1.5, shape=15)), shape=guide_legend(nrow=1,override.aes=list(size=1.5)))+main_theme+theme(legend.text=element_text(size=6), legend.position=c(0,0.1), legend.spacing = unit(-0.5, "cm")), nrow = 2),
             arrangeGrob(C_B_boxplot_COG_by_cat+main_theme+theme(axis.text.y = element_text(size=6), axis.text.x = element_text(size=6), axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6)), N_CONT_boxplot_COG_by_cat+main_theme+theme(axis.text.y = element_text(size=6), axis.text.x = element_text(size=6), axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6)), heights = c(1, 0.65)), 
             ncol=2, widths=c(0.45,1))
grid.text(c("a", "b", "c", "d"), c(0.02,0.3, 0.02,0.3), c(0.97, 0.97, 0.5,0.5), gp=gpar(col="black", cex=1),
          just=c("left"))
dev.off()
