setwd("..\\WGBS\\result\\realdata\\res")
install.packages("BiocManager")
install.packages("clusterProfiler")
BiocManager::install("DO.db")
BiocManager::install(version = "3.13")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db",force = TRUE)
BiocManager::install("org.Bt.eg.db")
BiocManager::install("org.Ss.eg.db",force = TRUE)
BiocManager::install("biomaRt")
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Bt.eg.db)
library(org.Ss.eg.db)
library("biomaRt")
library(clusterProfiler)


######### enrich #############
speEnrich<- function(species,mapper,ensemblData,orgdb,organismName){
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset=ensemblData)
  input <- read.table(paste(species,"\\gene\\",mapper,"_gene_DMR_Res.txt",sep=""),header = F, sep = "\t", quote = "\"", dec = ".")
  inputGene<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol', 'entrezgene_id'), 
                    filters = 'ensembl_gene_id', values = input$V1, mart = ensembl)
  gene <- inputGene$entrezgene_id
  gene <- gene[-which(is.na(gene))]
  allResKEGG <- enrichKEGG(gene,organism=organismName,keyType = 'kegg',pvalueCutoff=0.05, pAdjustMethod = 'BH')
  resKEGG <- allResKEGG@result
  resKEGG <- resKEGG[which(resKEGG$pvalue<0.05),]
  write.table(resKEGG,paste("enrichKEGG\\",species,"_",mapper,"_enrichKEGG.txt",sep=""),quote = F,col.names = F,row.names = F,sep = "\t")
}
speEnrich("pig","bismarkbwt2","sscrofa_gene_ensembl","org.Ss.eg.db","ssc")
speEnrich("pig","bsmap","sscrofa_gene_ensembl","org.Ss.eg.db","ssc")
speEnrich("pig","bwameth","sscrofa_gene_ensembl","org.Ss.eg.db","ssc")
speEnrich("pig","walt","sscrofa_gene_ensembl","org.Ss.eg.db","ssc")

speEnrich("mouse","bismarkbwt2","mmusculus_gene_ensembl","org.Mm.eg.db","mmu")
speEnrich("mouse","bsmap","mmusculus_gene_ensembl","org.Mm.eg.db","mmu")
speEnrich("mouse","bwameth","mmusculus_gene_ensembl","org.Mm.eg.db","mmu")
speEnrich("mouse","walt","mmusculus_gene_ensembl","org.Mm.eg.db","mmu")

speEnrich("human","bismarkbwt2","hsapiens_gene_ensembl","org.Hs.eg.db","hsa")
speEnrich("human","bsmap","hsapiens_gene_ensembl","org.Hs.eg.db","hsa")
speEnrich("human","bwameth","hsapiens_gene_ensembl","org.Hs.eg.db","hsa")
speEnrich("human","walt","hsapiens_gene_ensembl","org.Hs.eg.db","hsa")

speEnrich("cattle","bismarkbwt2","btaurus_gene_ensembl","org.Bt.eg.db","bta")
speEnrich("cattle","bsmap","btaurus_gene_ensembl","org.Bt.eg.db","bta")
speEnrich("cattle","bwameth","btaurus_gene_ensembl","org.Bt.eg.db","bta")
speEnrich("cattle","walt","btaurus_gene_ensembl","org.Bt.eg.db","bta")





