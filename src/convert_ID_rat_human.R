annotated <- read.table("/home/af/Desktop/Differential_analyses_peaks/data/processed/motif2gene_mapping.txt",sep="\t", stringsAsFactors = FALSE)

human <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
rat <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="rnorvegicus_gene_ensembl")



genesV2 = getLDS(attributes = c("ensembl_gene_id"), 
                 filters = "ensembl_gene_id", 
                 values = annotated$V2 , 
                 mart = human, 
                 attributesL = c("ensembl_gene_id"), 
                 martL = rat, 
                 uniqueRows=T)



merged = merge(genesV2, annotated, by.y=c("V2"), by.x="Gene.stable.ID")


write.table(merged[,c(3,2)],"test.txt", quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)




