
suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(DOSE))
suppressMessages(library(enrichplot))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library(magrittr))
suppressMessages(library(rWikiPathways))
library(biomaRt)
library(forcats)
library(yaml)

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


background = read.table(file.path(rootfolder,radio$background),sep = "\t", stringsAsFactor = FALSE, h = T)

  # recode the annotation to remove mention of the transcript
  background <-background[!is.na(background$Annotation),]
  background$Annotation[grep("promoter", background$Annotation)] <- "Promoter-TSS"
  background$Annotation[grep("TTS", background$Annotation)] <- "TTS"
  background$Annotation[grep("intron", background$Annotation)] <- "Intron"
  background$Annotation[grep("exon", background$Annotation)] <- "Exon"
  background$Annotation <- factor(background$Annotation, levels =rev(c('Promoter-TSS','Exon','Intron','TTS','Intergenic')))


#human <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
#rat <- useMart("ENSEMBL_MART_ENSEMBL", dataset="rnorvegicus_gene_ensembl")
human = readRDS(file.path(rootfolder,radio$humanbiomart))
rat = readRDS(file.path(rootfolder,radio$ratbiomart))


viz_ORA <- function(enrichR_output, database, GO = FALSE){
  
  # --------- Dotplot ----------
  dotplot_plot = dotplot(enrichR_output, showCategory=30) + ggtitle("dotplot for ORA")
  try(ggsave(file.path(rootfolder,radio$plot,paste0("GO_",database,"_dotplot_ora.png")),dotplot_plot,width = 10, height = 6))
  
  # ----- Enrichment map or goplot ------
  if(GO == FALSE){
    emapplot_plot = emapplot(enrichR_output)
 try(ggsave(file.path(rootfolder,radio$plot,paste0("GO_",database,"_enrichment_map_ora.png")),emapplot_plot,width = 8, height = 8))
    
  }else{
    goplot_plot = goplot(enrichR_output)
    try(ggsave(file.path(rootfolder,radio$plot,paste0("GO_",database,"_enrichment_map_ora_GO.png")),goplot_plot,width = 6, height = 6))
  }

  # --------- Heatmap ---------- 
  heatplot_plot = heatplot(enrichR_output)
  try(ggsave(file.path(rootfolder,radio$plot,paste0("GO_",database,"_heatmap_ora.png")),heatplot_plot,width = 24, height = 12))

  # ---- Network w/t genes ----- 

  cnetplot_plot = cnetplot(enrichR_output)
  try(ggsave(file.path(rootfolder,radio$plot,paste0("GO_",database,"_cnetplot_ora.png")),cnetplot_plot,width = 12, height = 12))
}



#######Conversion from RAT to human
#This leads to errors with too many extra convertion ( Interleukine over representation
#homerhuman = getLDS(attributes = c("ensembl_gene_id", "entrezgene_id","external_gene_name"), 
#                 filters = "ensembl_gene_id", 
#                 values = homer$Entrez.ID , 
#                 mart = rat, 
#                 attributesL = c("ensembl_gene_id"), 
#                 martL = human, 
#                 uniqueRows=T)

#names(homerhuman) <- c("rat","entrez","name","human")

 homerhuman_list <- bitr(toupper(homer$Gene.Name),fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID","SYMBOL"),OrgDb = "org.Hs.eg.db")



#backgroundhuman = getLDS(attributes = c("ensembl_gene_id", "entrezgene_id","external_gene_name"), 
#                 filters = "ensembl_gene_id", 
#                 values = background$Entrez.ID , 
#                 mart = rat, 
#                 attributesL = c("ensembl_gene_id"), 
#                 martL = human, 
#                 uniqueRows=T)
#names(backgroundhuman) <- c("rat","entrez","name","human")
background_list <- bitr(toupper(background$Gene.Name),fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID","SYMBOL"),OrgDb = "org.Hs.eg.db")





egoMF <- enrichGO(gene=homerhuman_list$ENTREZID,
                  OrgDb = 'org.Hs.eg.db', #the rat or human database
                  ont="MF", #BP pour biological process (on peut aussi avoir cellular component (CC) et molecular function (MF))
                  pvalueCutoff = 0.05,
		  pAdjustMethod = "none",
                  universe = background_list$ENTREZID,
                  qvalueCutoff = 0.1,
                  readable = T,#False
)

egoCC <- enrichGO(gene=homerhuman_list$ENTREZID,
                  OrgDb = 'org.Hs.eg.db', #the rat or human database
                  ont="CC", #BP pour biological process (on peut aussi avoir cellular component (CC) et molecular function (MF))
                  pvalueCutoff = 0.05,
		  pAdjustMethod = "bonferroni",
                  universe = background_list$ENTREZID,
                  qvalueCutoff = 0.1,
                  readable = T,#False
)

egoBP <- enrichGO(gene=homerhuman_list$ENTREZID,
                  OrgDb = 'org.Hs.eg.db', #the rat or human database
                  ont="BP", #BP pour biological process (on peut aussi avoir cellular component (CC) et molecular function (MF))
                  pvalueCutoff = 0.05,
		  pAdjustMethod = "none",
                  universe = background_list$ENTREZID,
                  qvalueCutoff = 1,
                  readable = T,#False
)

  dotplot_plot = dotplot(egoBP, showCategory=5) + ggtitle("dotplot for ORA")

ggplot(head(egoBP), # you can replace the numbers to the row number of pathway of your interest
             aes(x = GeneRatio, y = Description)) + 
             geom_point(aes(size = GeneRatio, color = p.adjust)) +
             theme_bw(base_size = 14) +
             scale_colour_gradient(limits=c(0, 0.10), low="red") +
             ylab(NULL) +
             ggtitle("GO pathway enrichment")



viz_ORA(egoBP,"BP",GO=TRUE)
viz_ORA(egoMF,"MF",GO=TRUE)
viz_ORA(egoCC,"CC",GO=TRUE)
#subset = homer[intersect(which(homer$Gene.Type == "protein_coding"), which(homer$Annotation != "Intergenic")),]
subset = homer[which(homer$Annotation != "Intergenic"),]
subsethuman_list <- bitr(toupper(subset$Gene.Name),fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID","SYMBOL"),OrgDb = "org.Hs.eg.db")



#Not much here : in subset nothing too significativ
subsetBP <- enrichGO(gene=subsethuman_list$ENTREZID,
                  OrgDb = 'org.Hs.eg.db', #the rat or human database
                  ont="BP", #BP pour biological process (on peut aussi avoir cellular component (CC) et molecular function (MF))
                  pvalueCutoff = 0.9,
		  pAdjustMethod = "BH",
                  universe = background_list$ENTREZID,
                  qvalueCutoff = 1,
                  readable = T,#False
)
subsetCC <- enrichGO(gene=subsethuman_list$ENTREZID,
                  OrgDb = 'org.Hs.eg.db', #the rat or human database
                  ont="CC", #BP pour biological process (on peut aussi avoir cellular component (CC) et molecular function (MF))
                  pvalueCutoff = 0.05,
		  pAdjustMethod = "bonferroni",
                  universe = background_list$ENTREZID,
                  qvalueCutoff = 0.2,
                  readable = T,#False
)
subsetMF <- enrichGO(gene=subsethuman_list$ENTREZID,
                  OrgDb = 'org.Hs.eg.db', #the rat or human database
                  ont="MF", #BP pour biological process (on peut aussi avoir cellular component (CC) et molecular function (MF))
                  pvalueCutoff = 1,
		  pAdjustMethod = "bonferroni",
                  universe = background_list$ENTREZID,
                  qvalueCutoff = 0.2,
                  readable = T,#False
)

  write.table(subsetBP,file.path(rootfolder,"data","processed","subsetBP.csv"),quote=F,col.names=T, row.names = F,sep="\t")
toplot = subsetBP@result
frac = sapply(subsetBP$GeneRatio, function(x) eval(parse(text=x)))
sel = toplot[c(2,3,6,8,9,11,13,14,15,16),]

sel = sel[order(sel$GeneRatio),]


 chocho <- mutate(sel,Description=factor(Description, levels=Description))

pp <- ggplot(chocho,# you can replace the numbers to the row number of pathway of your interest
             aes(x = GeneRatio, y = Description)) + 
             geom_point(aes(size = 50, color = pvalue)) +
             theme_bw(base_size = 14) +
             scale_colour_gradient(limits=c(0.002, 0.01), low="red") +
             ylab(NULL) +
             ggtitle("") +
 theme(axis.text.y = element_text(size = 20),axis.text.x = element_text(size = 20),legend.text=element_text(size=20), legend.title = element_text(size = 20)) + guides(size = "none")

ggsave(file.path(rootfolder,radio$plot,paste0("GO_subset_BP_dotplot_ora.pdf")), plot = pp,width = 12, height = 8)



######################################RNA
######################################RNA

changelimit = 0.09
pvallimit = 90

topgenes = read.table(file.path(rootfolder,radio$tobias_data), stringsAsFactor = FALSE, h=TRUE)
topgenes$diffexpressed <- rep("NO",length(topgenes$WT_KI_change))

topgenes$diffexpressed[topgenes$WT_KI_change > changelimit & -log10(topgenes$WT_KI_pvalue) > pvallimit] <- "UP"
topgenes$diffexpressed[topgenes$WT_KI_change < -changelimit & -log10(topgenes$WT_KI_pvalue) > pvallimit] <- "DOWN"
gene = toupper(topgenes$name[union(which(topgenes$diffexpressed == "UP"),which(topgenes$diffexpressed == "DOWN"))])




wp2gene <- read.gmt(file.path(rootfolder,radio$wikipathway))

wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
 tfbs <- bitr(toupper(gene),fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID","SYMBOL"),OrgDb = "org.Hs.eg.db")

ora_wp <- enricher(tfbs$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

ora_wp <- setReadable(ora_wp, org.Hs.eg.db, keyType = "ENTREZID")
ora_wp <- ora_wp@result
ora_wp <- ora_wp[c(1:5),]
ora_wp <- ora_wp[with(ora_wp,order(Count)),]

ora_wp <- mutate(ora_wp,Description=factor(Description, levels=Description))

ora_wp
pp <- ggplot(ora_wp,# you can replace the numbers to the row number of pathway of your interest
             aes(x = GeneRatio, y = Description)) + 
             geom_point(aes(size = 50, color = -log10(p.adjust))) +
             theme_bw(base_size = 14) +
             scale_colour_gradient(limits=c(0, 5), low="red") +
             ylab(NULL) +
             ggtitle("") + labs(color = expression(paste("-Log"[10],"(pValue)")))+
 theme(axis.text.y = element_text(size = 25),axis.text.x = element_text(size = 20),axis.title.x = element_text(size = 20),legend.text=element_text(size=20), legend.title = element_text(size = 20)) + guides(size = "none")

ggsave(file.path(rootfolder,radio$plot,paste0("Wikipathway_tfbs.pdf")), plot = pp,width = 14, height = 8)


panther = read.table(file.path(rootfolder,radio$panther),sep="\t", stringsAsFactor = FALSE, h=TRUE)

panther <- panther[c(1:5),]
panther$Count = panther$GeneRatio
names(panther) <-c("Description","GeneRatio", "pvalue", "p.adjust","Old.P.value","OldAdjusted.P.value", "Odds.Ratio","Combined.Score","Genes")

panther$Description = gsub(" P[0-9]*","", panther$Description)
panther$Description = gsub(" Homo sapiens","", panther$Description)

tt = strsplit(panther$GeneRatio,"/")
gg = c()
for(i in tt){gg = c(as.numeric(i[1])/as.numeric(i[2]),gg)}
panther$Ratio <- rev(gg) # here reverse because it does it reverse, dont know why
pp <- ggplot(panther,# you can replace the numbers to the row number of pathway of your interest
             aes(x = GeneRatio, y = Description)) + 
             geom_point(aes(size = Ratio, color = -log10(p.adjust))) +
             theme_bw(base_size = 14) +
             scale_colour_gradient(limits=c(0, 5), low="red") +
             ylab(NULL) +
             ggtitle("") + labs(color =  expression(paste("-Log"[10],"(pValue)")))+
 theme(axis.text.y = element_text(size = 25),axis.title.x = element_text(size = 20),axis.text.x = element_text(size = 20),legend.text=element_text(size=20), legend.title = element_text(size = 20)) + guides(size = "none") + scale_size(range = c(10, 20))
pp
ggsave(file.path(rootfolder,radio$plot,paste0("PantherDB_tfbs.pdf")), plot = pp,width = 14, height = 8)


####GO MF
gomf = read.table(file.path(rootfolder,radio$GO_MF_diff),sep="\t", stringsAsFactor = FALSE, h=TRUE)

names(gomf) <-c("Description","GeneRatio", "pvalue", "p.adjust","Old.P.value","OldAdjusted.P.value", "Odds.Ratio","Combined.Score","Genes")
gomf$Count = gomf$GeneRatio
gomf <- gomf[c(1:5),]


tt = strsplit(gomf$GeneRatio,"/")
gg = c()
for(i in tt){gg = c(as.numeric(i[1])/as.numeric(i[2]),gg)}
gomf$Ratio <- rev(gg) # here reverse because it does it reverse, dont know why

gomf <- gomf[with(gomf,order(Ratio, decreasing = TRUE)),]
gomf$Description = gsub(" \\(GO:[0-9]*)","",gomf$Description)
gomf <- mutate(gomf,Description=factor(Description, levels=Description))

pp <- ggplot(gomf,# you can replace the numbers to the row number of pathway of your interest
             aes(x = GeneRatio, y = Description)) + 
             geom_point(aes(size = Ratio, color = pvalue)) +
             theme_bw(base_size = 14) +
             scale_colour_gradient(limits=c(0, 0.05), low="red") +
             ylab(NULL) +
             ggtitle("") + labs(color =  "pValue")+
 theme(axis.text.y = element_text(size = 20),axis.text.x = element_text(size = 20),legend.text=element_text(size=25), legend.title = element_text(size = 20)) + guides(size = "none") + scale_size(range = c(10, 20))
pp
ggsave(file.path(rootfolder,radio$plot,paste0("Go_MF_Diff.pdf")), plot = pp,width = 14, height = 8)



wp2gene = read.table(file.path(rootfolder,radio$wiki_diff),sep="\t", stringsAsFactor = FALSE, h=TRUE)

names(wp2gene) <-c("Description","GeneRatio", "pvalue", "p.adjust","Old.P.value","OldAdjusted.P.value", "Odds.Ratio","Combined.Score","Genes")
wp2gene$Count = wp2gene$GeneRatio
wp2gene <- wp2gene[c(1:5),]


tt = strsplit(wp2gene$GeneRatio,"/")
gg = c()
for(i in tt){gg = c(as.numeric(i[1])/as.numeric(i[2]),gg)}
wp2gene$Ratio <- rev(gg) # here reverse because it does it reverse, dont know why
wp2gene <- wp2gene[with(wp2gene,order(pvalue, decreasing = TRUE)),]
wp2gene <- mutate(wp2gene,Description=factor(Description, levels=Description))


pp <- ggplot(wp2gene,# you can replace the numbers to the row number of pathway of your interest
             aes(x = GeneRatio, y = Description)) + 
             geom_point(aes(size = Ratio, color = pvalue)) +
             theme_bw(base_size = 14) +
             scale_colour_gradient(limits=c(0, 0.2), low="red") +
             ylab(NULL) +
             ggtitle("") + labs(color = "pValue")+
 theme(axis.text.y = element_text(size = 20),axis.text.x = element_text(size = 20),legend.text=element_text(size=25), legend.title = element_text(size = 20)) + guides(size = "none") + 
scale_size(range = c(10, 20))

ggsave(file.path(rootfolder,radio$plot,paste0("Wikipathway_Diff.pdf")), plot = pp,width = 14, height = 8)





