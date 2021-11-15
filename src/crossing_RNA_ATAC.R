
library(ggplot2)
library(dplyr)
library(yaml)
library(DESeq2)
library(biomaRt)
library(VennDiagram)
library(forcats)

radio = yaml.load_file(file.path("/home/af/Desktop/Differential_analyses_peaks/radiofile.yml"))
rootfolder = radio$rootfolder

topgenes = read.table(file.path(rootfolder,radio$topgenes_RNA), stringsAsFactor = FALSE)

homer_hugo = read.table(file.path(rootfolder,radio$homer_hugo), stringsAsFactor = FALSE)


homer = read.table(file.path(rootfolder,radio$homer),sep = "\t", stringsAsFactor = FALSE, h = T)

  
  # recode the annotation to remove mention of the transcript
  homer <-homer[!is.na(homer$Annotation),]
  homer$Annotation[grep("promoter", homer$Annotation)] <- "Promoter-TSS"
  homer$Annotation[grep("TTS", homer$Annotation)] <- "TTS"
  homer$Annotation[grep("intron", homer$Annotation)] <- "Intron"
  homer$Annotation[grep("exon", homer$Annotation)] <- "Exon"
  homer$Annotation <- factor(homer$Annotation, levels =rev(c('Promoter-TSS','Exon','Intron','TTS','Intergenic')))

head(homer)
homer$ID <-substr(homer$PeakID,1,2)
homer <- homer %>% mutate(ID=factor(ID, levels=c("WT","KI")))

# Changer les couleurs manuellement
ggplot(data=homer, aes(x=Annotation, fill=ID)) +
geom_bar(color="black", position=position_dodge())+
    theme_classic() + scale_fill_manual(values=c('#E69F00','#999999'))
# Couleurs personnalisées

ggsave(file.path(rootfolder,radio$plot,"Homer_annotation_373.png"))


ggplot(data=homer, aes(x=ID, fill=Annotation)) +
geom_bar(color="black", position="fill")+
    theme_classic() + scale_fill_brewer(palette="Set1")+ 
theme(axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.ticks.x=element_blank())


ggsave(file.path(rootfolder,radio$plot,"Homer_annotation_stacked_373.png") )

ggplot(data=homer, aes(x=Annotation, fill=Annotation)) +
geom_bar(color="black", position="stack")+
    theme_classic() + scale_fill_brewer(palette="Set1") + 
theme(axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.ticks.x=element_blank())


#########Checking stats
####Testing differences in Promoter-TSS
table(homer$Annotation[which(homer$ID =="WT")])
table(homer$Annotation[which(homer$ID =="KI")])


######Using all peaks from occupancy set
WT <- read.table(file.path(rootfolder,radio$WT_homer),sep = "\t", stringsAsFactor = FALSE, h = T)
  # recode the annotation to remove mention of the transcript
  WT <-WT[!is.na(WT$Annotation),]
  WT$Annotation[grep("promoter", WT$Annotation)] <- "Promoter-TSS"
  WT$Annotation[grep("TTS", WT$Annotation)] <- "TTS"
  WT$Annotation[grep("intron", WT$Annotation)] <- "Intron"
  WT$Annotation[grep("exon", WT$Annotation)] <- "Exon"
  WT$Annotation <- factor(WT$Annotation, levels =rev(c('Promoter-TSS','Exon','Intron','TTS','Intergenic')))
WT$ID <- "WT"

KI <- read.table(file.path(rootfolder,radio$KI_homer),sep = "\t", stringsAsFactor = FALSE, h = T)
  # recode the annotation to remove mention of the transcript
  KI <-KI[!is.na(KI$Annotation),]
  KI$Annotation[grep("promoter", KI$Annotation)] <- "Promoter-TSS"
  KI$Annotation[grep("TTS", KI$Annotation)] <- "TTS"
  KI$Annotation[grep("intron", KI$Annotation)] <- "Intron"
  KI$Annotation[grep("exon", KI$Annotation)] <- "Exon"
  KI$Annotation <- factor(KI$Annotation, levels =rev(c('Promoter-TSS','Exon','Intron','TTS','Intergenic')))
KI$ID <- "KI"

KIWT <- rbind(KI,WT)
KIWT <- KIWT %>% mutate(ID=factor(ID, levels=c("WT","KI")))
ggplot(data=KIWT, aes(x=ID, fill=Annotation)) +
geom_bar(color="black", position="fill")+
    theme_classic() + scale_fill_brewer(palette="Set1") +
theme(axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.ticks.x=element_blank())
    
ggsave(file.path(rootfolder,radio$plot,"Homer_annotation_KI_8135_WT_7562.pdf"),width=10, height=10,units="cm",dpi=100)

####Testing differences in Promoter-TSS
table(KIWT$Annotation[which(KIWT$ID =="WT")])
table(KIWT$Annotation[which(KIWT$ID =="KI")])



ggplot(data=homer[which(homer$ID=='WT'),], aes(x="",fill=Annotation)) + 
geom_bar(color="black",width = 1) +
coord_polar("y") +
  scale_fill_brewer(palette="Set1")+
theme_void()
ggsave(file.path(rootfolder,radio$plot,"Homer_annotation_WT.png"))


ggplot(data=homer[which(homer$ID=='KI'),], aes(x="",fill=Annotation)) + 
geom_bar(color="black",width = 1) +
coord_polar("y") +
  scale_fill_brewer(palette="Set1")+
theme(axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(file.path(rootfolder,radio$plot,"Homer_annotation_KI.png"))

rna <- read.table(file.path(rootfolder,radio$intersect_rna_atac),sep = "\t", stringsAsFactor = FALSE, h = T)

  rna <-rna[!is.na(rna$Annotation),]
  rna$Annotation[grep("promoter", rna$Annotation)] <- "Promoter-TSS"
  rna$Annotation[grep("TTS", rna$Annotation)] <- "TTS"
  rna$Annotation[grep("intron", rna$Annotation)] <- "Intron"
  rna$Annotation[grep("exon", rna$Annotation)] <- "Exon"
   rna$Annotation <- factor(rna$Annotation, levels =rev(c('Promoter-TSS','Exon','Intron','TTS','Intergenic')))
rna$ID <- "RNA"

ggplot(data=rna, aes(x=ID, fill=Annotation)) +
geom_bar(color="black", position="fill")+
    theme_classic() + scale_fill_brewer(palette="Set1")+
theme(axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(file.path(rootfolder,radio$plot,"RNA_Homer_annotation.pdf") ,width=10, height=10,units="cm",dpi=100)
##Test 
homer$ID <- "Differential"

KIWTDIFF <- rbind(KIWT,homer)

KIWTDIFF <- KIWTDIFF %>% mutate(ID=factor(ID, levels=c("WT","KI","Differential")))
ggplot(data=KIWTDIFF, aes(x=ID, fill=Annotation)) +
geom_bar(color="black", position="fill")+
    theme_classic() + scale_fill_brewer(palette="Set1")+
theme(axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(file.path(rootfolder,radio$plot,"Homer_annotation_KI_8135_WT_7562_Diff_373.pdf") ,width=10, height=10,units="cm",dpi=100)


###VennDiagram for RNA 

vv <- list(diff = homer$PeakID,RNA=c(paste("gene",1:523-length(rna$PeakID)),rna$PeakID))

venn.diagram(vv, filename=file.path(rootfolder,radio$plot,"Venn_RNA_Diff.pdf"),fill = c("#999999", "#009E73"),category.names = c("RNA-Seq" , "ATAC-Seq"),width=10, height=10,units="cm",dpi=100)
write.table(vv, filename="test.txt")







