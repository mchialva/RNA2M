# colours
dark_green <-"#32C864B3"
very_dark_green<-"#32C864B3"
sea_green <-"#2E815AB3"
grey <-"#B4B4B4B3"
dark_brown <-"#654321B3"
red <-"#C80000B3"
orange <-"#FF8200B3"
pale_blue <-"#1464C8B3"
light_blue<-"#00AFBB"
pale_red<-"#FC4E07"

# ggplot2 themes
main_theme<-theme(panel.grid=element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  axis.title.x = element_text(size = 6, vjust= -0.5, colour ="black"),
                  axis.title.y = element_text(size = 6, vjust= -0.5, colour ="black"),
                  axis.text.x = element_text(size=6, angle=0, colour ="black"),
                  axis.text.y = element_text(size=6, angle=0, colour ="black"),
                  axis.ticks=element_line(size=0.2),
                  legend.text=element_text(size=4),
                  legend.title=element_blank(),
                  axis.line.x = element_line(size=0.2, color = "black"),
                  axis.line.y = element_line(size=0.2, color = "black"),
                  legend.key = element_blank(),
                  legend.position=c(0.88,0.25),
                  legend.key.size = unit(0.2, "cm"),
                  legend.spacing=unit(-0.1,"cm"),
                  text=element_text(family="sans"),
                  plot.title=element_text(hjust = 0.5, size=6))

boxplot_theme<-theme(legend.key = element_blank(),
                     legend.position=c(0.2,0.18),
                     legend.key.size = unit(0.4, "cm"),
                     legend.spacing=unit(-0.1,"cm"),
                     legend.title=element_blank(),
                     text=element_text(family="sans"),
                     legend.text=element_text(size=5),
                     axis.line.x = element_line(size=0.2, color = "black"),
                     axis.line.y = element_line(size=0.2, color = "black"),
                     axis.text.x = element_text(size=6, angle=0, colour ="black"),
                     axis.text.y = element_text(size=6, angle=0, colour ="black"),
                     axis.ticks=element_line(size=0.2),
                     axis.title.x = element_text(size = 6, vjust= -0.5, colour ="black"),
                     axis.title.y = element_text(size = 6, vjust= -0.5, colour ="black"),
                     panel.grid=element_blank(),
                     strip.text =element_text(size = 6),
                     panel.background = element_blank(),
                     plot.margin = unit(c(0.1,0,0,0), "cm"))

barplot_theme<-theme(legend.key = element_blank(),
                     #legend.position="bottom",
                     legend.position=c(0.4,-0.35),
                     legend.box = "horizontal",
                     axis.text.x = element_text(size=4, angle=90,hjust=1,vjust=0.5),
                     axis.text.y = element_text(size=6, angle=0, colour ="black"),
                     axis.line.x = element_line(size=0.2, color = "black"),
                     axis.line.y = element_line(size=0.2, color = "black"),
                     legend.title=element_blank(),
                     panel.grid=element_blank(),
                     panel.background = element_blank(),
                     plot.margin = unit(c(0.1,0,1.1,0.1), "cm"),
                     legend.text=element_text(size=3.5),
                     axis.ticks=element_line(size=0.2),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size = 6, vjust= -0.5, colour ="black"),
                     text=element_text(family="sans"),
                     legend.key.size = unit(0.3, "cm"),
                     legend.spacing=unit(0,"cm"),
                     legend.margin =margin(0,0,0,0,"cm"),
                     panel.border = element_blank(),
                     plot.title=element_text(hjust = 0.5, size=6))
                     
MA_theme<-theme(panel.grid=element_blank(),
                  panel.background = element_blank(),
                  axis.title.x = element_text(size = 6, vjust= -0.5, colour ="black"),
                  axis.title.y = element_text(size = 6, vjust= -0.5, colour ="black",),
                  axis.text.x = element_text(size=6, angle=0, colour ="black"),
                  axis.text.y = element_text(size=6, angle=0, colour ="black"),
                  axis.ticks=element_line(size=0.2),
                  #element_text(size=0),
                  axis.line.x = element_line(size=0.2, color = "black"),
                  axis.line.y = element_line(size=0.2, color = "black"),
                  text=element_text(family="sans"),
                  plot.margin = unit(c(0,0.1,0,0.1), "cm"),
                  plot.title=element_text(hjust = 0.5, size=4, vjust =-1))

MA_theme_large<-MA_theme+theme(
                  plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
                  plot.title=element_text(hjust = 0.5, size=6))

bubbles_theme<-theme(panel.grid=element_blank(),
                  panel.background = element_blank(),
                  axis.title.x = element_text(size = 6, colour ="black"),
                  axis.title.y = element_text(size = 6, colour ="black"),
                  axis.text.x = element_text(size=6, angle=0),
                  axis.text.y = element_text(size=6, angle=0),
                  axis.ticks=element_line(size=0.2),
                  legend.text=element_text(size=4.5),
                  legend.title=element_text(size=5),
                  axis.line.x = element_line(size=0.2, color = "black"),
                  axis.line.y = element_line(size=0.2, color = "black"),
                  legend.key = element_blank(),
                  text=element_text(family="sans"),
                  legend.margin=margin(0,0,0,0,"cm"),
                  legend.key.size = unit(0.2, "cm"),
                  plot.margin = unit(x = c(0, 0, 0, 0), units = "mm"),
                  strip.text.x=element_text(size = 6))

no_y_theme<-theme(axis.line.y = element_line(size=0.2, color=NA),
                  axis.text.y = element_text(size=6, color = NA),
                  axis.title.y = element_text(size=6, color = NA),
                  axis.ticks.y = element_line(size=0.2, color=NA),
                  axis.line.x = element_line(size=0.2, color = "black"),
                  plot.margin = unit(c(0.1,0.5,1.1,-0.5), "cm"))

no_x_theme<-theme(axis.line.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.ticks.x = element_blank())