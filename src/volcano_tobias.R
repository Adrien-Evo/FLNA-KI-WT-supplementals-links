
library(ggplot2)
#library(dplyr)
library(yaml)
library(DESeq2)
library(biomaRt)
library(ggrepel)

radio = yaml.load_file(file.path("/home/af/Desktop/Differential_analyses_peaks/radiofile.yml"))
rootfolder = radio$rootfolder

changelimit = 0.1
pvallimit = 100

topgenes = read.table(file.path(rootfolder,radio$tobias_data), stringsAsFactor = FALSE, h=TRUE)
topgenes$diffexpressed <- rep("NO",length(topgenes$WT_KI_change))

topgenes$diffexpressed[topgenes$WT_KI_change > changelimit & -log10(topgenes$WT_KI_pvalue) > pvallimit] <- "UP"
topgenes$diffexpressed[topgenes$WT_KI_change < -changelimit & -log10(topgenes$WT_KI_pvalue) > pvallimit] <- "DOWN"

topgenes$labels <- rep(NA,length(topgenes$WT_KI_change))
topgenes$labels[topgenes$diffexpressed != "NO"] <- topgenes$name[topgenes$diffexpressed != "NO"]

###REnaming of FOS gene

topgenes$labels[intersect(grep("FOS",topgenes$name),grep("::",topgenes$name))] <- NA
#Adjusting JUN::JUNB set 

topgenes$diffexpressed <- as.factor(topgenes$diffexpressed)

p <-ggplot(data=topgenes, aes(x=WT_KI_change, y=-log10(WT_KI_pvalue),label=labels)) + 
    geom_point(aes(colour = topgenes$diffexpressed))+ 
    theme_classic() + 
    geom_text_repel(aes(colour = topgenes$diffexpressed),show.legend = FALSE) + 
    scale_colour_manual(values=c("red","grey","blue"),labels=c("Higher scores in KI", "No changes", "Higher scores in WT")) +
    geom_vline(xintercept=c(-changelimit, changelimit),linetype = "dashed",  alpha=0.4) +
    geom_hline(yintercept=pvallimit, linetype = "dashed", alpha=0.4)  +
    theme(legend.position = c(0.87,0.2),legend.key.size = unit(0.4, "cm")) + 
    labs(x ="Change", y = expression("Significance (-Log[10])")) + 
    expand_limits(x=c(-0.2,0.2)) + 
    guides(colour=guide_legend(title=NULL))



ggsave(file.path(rootfolder,radio$plot,"volcano_tobias.pdf"),p,width=15, height=15,units="cm",dpi=100)

