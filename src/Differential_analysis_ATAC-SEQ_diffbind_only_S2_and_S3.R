###Works with the DiffBind environnement

library(DiffBind)
library(ggplot2)
#library(dplyr)
library(yaml)
library(DESeq2)


rootfolder = "/home/af/Desktop/Differential_analyses_peaks/"
radio = yaml.load_file(file.path(rootfolder,"radiofile.yml"))

sample=read.csv(file.path(rootfolder,radio$samplesheets),h=T,sep="\t", stringsAsFactor = FALSE)

sample$Condition <- paste0(sample$Tissue,"_",sample$Factor)
## REMOVING samples S1
##
sample <- sample[-c(which(sample$SampleID =="KI_1"),which(sample$SampleID =="WT_1")),]



###Doing the CM vs IPS analysis


###QC selection
DBsample=dba(sampleSheet=sample)

##Plot folder
setwd(file.path(rootfolder,radio$plotDiffBindQC))

#######
####### OCCUPANCY ANALYSIS
#######

png("Clustering_Occupancy.png", width= 840, height = 840)
plot(DBsample)
dev.off()



###Overlap
KI_overlap <- dba.overlap(DBsample,DBsample$masks$KI ,mode=DBA_OLAP_RATE)
png("Overlap_KI.png", width= 840, height = 840)
plot(KI_overlap,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets for KI ')
dev.off()


WT_overlap <- dba.overlap(DBsample,DBsample$masks$WT ,mode=DBA_OLAP_RATE)

png("Overlap_WT.png", width= 840, height = 840)
plot(WT_overlap,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets for WT ')
dev.off()

# ===>Needs complete overlap

DBsample_simple <- dba.peakset(DBsample, consensus=c(DBA_TREATMENT), minOverlap=1)



png("Venn_KI_WT_consensus peakset.png", width= 840, height = 840)

dba.plotVenn(DBsample_simple,DBsample_simple$masks$Consensus)

dev.off()

venn = dba.overlap(DBsample_simple, DBsample_simple$masks$Consensus)

write_report_occupancy <- function(report,reportname){

df <- data.frame(seqnames=seqnames(report),
    starts=start(report)-1,
    ends=end(report),
    names=c(paste0(reportname,"_peak",names(report))),
    scores=c(rep("0", length(report))),
    strands=c(rep(".", length(report)))
            )
  write.table(df,file.path(rootfolder,radio$dataDiffBindQC,paste0(reportname,"_Occupancy.csv")),quote=F,row.names=F,col.names=F,sep="\t")
}


write_report_occupancy(venn$onlyA,"KI")
write_report_occupancy(venn$onlyB,"WT")
write_report_occupancy(venn$inAll,"common")



#################### Differential
################## Affinity analysis
##################################################################
DBcount <- dba.count(DBsample)

saveRDS(DBcount, file = file.path(rootfolder,radio$dataDiffBindQC,"dba.count.rds"))


png("PCA_TREATMENT.png", width= 840, height = 840)

dba.plotPCA(DBcount,attributes=DBA_TREATMENT,title="KI vs WT")

dev.off()

png("PCA_REPLICATE.png", width= 840, height = 840)

dba.plotPCA(DBcount,attributes=DBA_REPLICATE,title="Replicate")

dev.off()


png("Clustering_Affinity.png", width= 840, height = 840)
plot(DBcount)
dev.off()


contrast <- dba.contrast(DBcount, categories=DBA_TREATMENT,minMembers=2)
##Contrast all methods : EDGER 
KIvsWT<- dba.analyze(contrast, method=c(DBA_EDGER))


png("MA_plot_KI_WT_EDGER.png", width= 840, height = 840)

dba.plotMA(KIvsWT)
dev.off()

png("volcano_plot_KI_WT_pval_EDGER.png", width= 840, height = 840)
dba.plotVolcano(KIvsWT, bUsePval=TRUE, th = 0.01, fold=1)
dev.off()


png("box_plot_KI_WT_FDR_EDGER.png", width= 840, height = 840)
dba.plotBox(KIvsWT, bUsePval=FALSE, th = 0.01, fold=1)
dev.off()


KIvsWT.report <- dba.report(KIvsWT, th = 1,bCalled=TRUE,method=DBA_EDGER)
##Reordering df

KIvsWT.report <- sort(KIvsWT.report)
###########Function to write the report. Here fold = fold_cond1 -fold_cond2

#report = KIvsWT.report
#reportname = "KI_vs_WT"
#cond1 = "KI"
#cond2 = "WT"
#FDR = 0.01
#fold = 1
write_report <- function(report,reportname,cond1,cond2,FDR,fold){
  df <- data.frame(chr=seqnames(report),
    starts=start(report)-1,
    ends=end(report),
    names=c(paste0("peak",1:length(report))),
    scores=c(rep("0", length(report))),
    strands=c(rep(".", length(report))),
    conc = report$Conc,
    conc_cond1=mcols(report)[,2],
    conc_cond2=mcols(report)[,3],
    fold =report$Fold,
    pval = report$"p-value",
    FDR = report$FDR
            )
colnames(df)[8] <- paste0(cond1,"_conc")
colnames(df)[9] <- paste0(cond2,"_conc")

write.table(df,file.path(rootfolder,radio$dataDiffBindQC,paste0(reportname,"_full.csv")),quote=F,col.names=T, row.names = F,sep="\t")

  temp = as.character(df$names)
  temp[which(df$fold <= 0)] <- paste0(cond2,"_",temp[which(df$fold <= 0)])
  temp[which(df$fold > 0)] <- paste0(cond1,"_",  temp[which(df$fold > 0)])

  df$names=temp
  ###Selecting based on asb(fold change) > 1 and FDR < 1%
  df_select <- df[which(abs(df$fold) >= fold),]

  df_select <- df_select[which(df_select$FDR <= FDR),]
write.table(df_select[,c(1,2,3,4)],file.path(rootfolder,radio$dataDiffBindQC,paste0(reportname,"_fdr-",FDR,"_fold-",fold,".csv")),quote=F,row.names=F,col.names=F,sep="\t")

}

write_report(KIvsWT.report,"KI_vs_WT_EDGER","KI","WT",0.01,1)

##Contrast all methods : DBA_DESEQ2 WONT works since it doesnt find any peaks

KIvsWT<- dba.analyze(contrast, method=c(DBA_DESEQ2))

png("MA_plot_KI_WT_DESEQ2.png", width= 840, height = 840)

dba.plotMA(KIvsWT)
dev.off()

png("volcano_plot_KI_WT_pval_DESEQ2.png", width= 840, height = 840)
dba.plotVolcano(KIvsWT, bUsePval=TRUE, th = 0.01, fold=1)
dev.off()


png("box_plot_KI_WT_FDR_DESEQ2.png", width= 840, height = 840)
dba.plotBox(KIvsWT, bUsePval=FALSE, th = 0.01, fold=1)
dev.off()


KIvsWT.report <- dba.report(KIvsWT, th = 1,bCalled=TRUE, method=DBA_DESEQ2)
##Reordering df

KIvsWT.report <- sort(KIvsWT.report)

write_report(KIvsWT.report,"KI_vs_WT_DESEQ2","KI","WT",0.01,1)






