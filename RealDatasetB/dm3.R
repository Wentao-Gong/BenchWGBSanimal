library(data.table)
library(DSS)
require(bsseq)

dmAnalysis <- function(file1.1,file1.2,file1.3,file2.1,file2.2,file2.3,dmlResFileName,dmrResFileName){
  dat1.1 <- fread(file1.1, header=T,stringsAsFactors = FALSE)
  dat1.2 <- fread(file1.2, header=T,stringsAsFactors = FALSE)
  dat1.3 <- fread(file1.3, header=T,stringsAsFactors = FALSE)
  dat2.1 <- fread(file2.1, header=T,stringsAsFactors = FALSE)
  dat2.2 <- fread(file2.2, header=T,stringsAsFactors = FALSE)
  dat2.3 <- fread(file2.3, header=T,stringsAsFactors = FALSE)
  BSobj <- makeBSseqData( list(dat1.1,dat1.2,dat1.3, dat2.1,dat2.2,dat2.3),c("first1","first2","first3","second1","second2","second3") )
  dmlTest.sm <- DMLtest(BSobj,group1 = c("first1","first2","first3"),group2 = c("second1","second2","second3"),smoothing = TRUE)
  dmls <- callDML(dmlTest.sm, p.threshold=0.05)
  dmrs <- callDMR(dmlTest.sm, p.threshold=0.05)
  write.table(dmls,dmlResFileName,quote = F,col.names = T,row.names = F)
  write.table(dmrs,dmrResFileName,quote = F,col.names = T,row.names = F)
}

## human: Differential methylation analysis
setwd("../meth/human2")
dmAnalysis("cpgSite/bwameth_SRR6373923_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR6825466_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR6825471_depth10_CpGofDSS.txt",
           "cpgSite/bwameth_SRR6818517_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR6373926_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR6373932_depth10_CpGofDSS.txt",
           "dssRes2/bwameth_depth10_DML.txt","dssRes2/bwameth_depth10_DMR.txt")
dmAnalysis("cpgSite/bsmap_SRR6373923_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR6825466_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR6825471_depth10_CpGofDSS.txt",
           "cpgSite/bsmap_SRR6818517_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR6373926_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR6373932_depth10_CpGofDSS.txt",
           "dssRes2/bsmap_depth10_DML.txt","dssRes2/bsmap_depth10_DMR.txt")
dmAnalysis("cpgSite/bismarkbwt2_SRR6373923_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR6825466_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR6825471_depth10_CpGofDSS.txt",
           "cpgSite/bismarkbwt2_SRR6818517_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR6373926_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR6373932_depth10_CpGofDSS.txt",
          "dssRes2/bismarkbwt2_depth10_DML.txt","dssRes2/bismarkbwt2_depth10_DMR.txt")
dmAnalysis("cpgSite/walt_SRR6373923_depth10_CpGofDSS.txt","cpgSite/walt_SRR6825466_depth10_CpGofDSS.txt","cpgSite/walt_SRR6825471_depth10_CpGofDSS.txt",
           "cpgSite/walt_SRR6818517_depth10_CpGofDSS.txt","cpgSite/walt_SRR6373926_depth10_CpGofDSS.txt","cpgSite/walt_SRR6373932_depth10_CpGofDSS.txt",
           "dssRes2/walt_depth10_DML.txt","dssRes2/walt_depth10_DMR.txt")

## cattle: Differential methylation analysis
setwd("../meth/cattle2")
dmAnalysis("cpgSite/bwameth_SRR7528450_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR7528456_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR7528458_depth10_CpGofDSS.txt",
           "cpgSite/bwameth_SRR7528459_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR7528464_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR7528465_depth10_CpGofDSS.txt",
           "dssRes2/bwameth_depth10_DML.txt","dssRes2/bwameth_depth10_DMR.txt")
dmAnalysis("cpgSite/bsmap_SRR7528450_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR7528456_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR7528458_depth10_CpGofDSS.txt",
           "cpgSite/bsmap_SRR7528459_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR7528464_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR7528465_depth10_CpGofDSS.txt",
           "dssRes2/bsmap_depth10_DML.txt","dssRes2/bsmap_depth10_DMR.txt")
dmAnalysis("cpgSite/bismarkbwt2_SRR7528450_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR7528456_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR7528458_depth10_CpGofDSS.txt",
           "cpgSite/bismarkbwt2_SRR7528459_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR7528464_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR7528465_depth10_CpGofDSS.txt",
           "dssRes2/bismarkbwt2_depth10_DML.txt","dssRes2/bismarkbwt2_depth10_DMR.txt")
dmAnalysis("cpgSite/walt_SRR7528450_depth10_CpGofDSS.txt","cpgSite/walt_SRR7528456_depth10_CpGofDSS.txt","cpgSite/walt_SRR7528458_depth10_CpGofDSS.txt",
           "cpgSite/walt_SRR7528459_depth10_CpGofDSS.txt","cpgSite/walt_SRR7528464_depth10_CpGofDSS.txt","cpgSite/walt_SRR7528465_depth10_CpGofDSS.txt",
           "dssRes2/walt_depth10_DML.txt","dssRes2/walt_depth10_DMR.txt")

## pig: Differential methylation analysis
setwd("../meth/pig2")
dmAnalysis("cpgSite/bwameth_SRR7812176_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR7812178_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR7812179_depth10_CpGofDSS.txt",
           "cpgSite/bwameth_SRR7812199_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR7812200_depth10_CpGofDSS.txt","cpgSite/bwameth_SRR7812210_depth10_CpGofDSS.txt",
           "dssRes2/bwameth_depth10_DML.txt","dssRes2/bwameth_depth10_DMR.txt")
dmAnalysis("cpgSite/bsmap_SRR7812176_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR7812178_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR7812179_depth10_CpGofDSS.txt",
           "cpgSite/bsmap_SRR7812199_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR7812200_depth10_CpGofDSS.txt","cpgSite/bsmap_SRR7812210_depth10_CpGofDSS.txt",
           "dssRes2/bsmap_depth10_DML.txt","dssRes2/bsmap_depth10_DMR.txt")
dmAnalysis("cpgSite/bismarkbwt2_SRR7812176_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR7812178_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR7812179_depth10_CpGofDSS.txt",
           "cpgSite/bismarkbwt2_SRR7812199_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR7812200_depth10_CpGofDSS.txt","cpgSite/bismarkbwt2_SRR7812210_depth10_CpGofDSS.txt",
           "dssRes2/bismarkbwt2_depth10_DML.txt","dssRes2/bismarkbwt2_depth10_DMR.txt")
dmAnalysis("cpgSite/walt_SRR7812176_depth10_CpGofDSS.txt","cpgSite/walt_SRR7812178_depth10_CpGofDSS.txt","cpgSite/walt_SRR7812179_depth10_CpGofDSS.txt",
           "cpgSite/walt_SRR7812199_depth10_CpGofDSS.txt","cpgSite/walt_SRR7812200_depth10_CpGofDSS.txt","cpgSite/walt_SRR7812210_depth10_CpGofDSS.txt",
           "dssRes2/walt_depth10_DML.txt","dssRes2/walt_depth10_DMR.txt")
