setwd("..\\result\\depth5")

library(ggplot2)
library(RColorBrewer)
################### Figure 4a: The class of unsatisfactory mapped reads #######################

errRes <- function(error){
  speciesRes<- function(species,species2,error,allReadsNum){
    report1 <- read.csv(paste("depth5Err",error,"\\",species,"1\\readId\\nonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    report2 <- read.csv(paste("depth5Err",error,"\\",species,"2\\readId\\nonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    report3 <- read.csv(paste("depth5Err",error,"\\",species,"3\\readId\\nonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    res <- data.frame(Mapper=c("Bismark-bwt2-e2e","Bismark-bwt2-e2e","Bismark-bwt2-e2e",
                               "BSBolt","BSBolt","BSBolt",
                               "BSMAP","BSMAP","BSMAP",
                               "Bwa-meth","Bwa-meth","Bwa-meth",
                               "Walt","Walt","Walt"),
                      Class=c("Error mapped reads","Multiple mapped reads","Discarded reads",
                              "Error mapped reads","Multiple mapped reads","Discarded reads",
                              "Error mapped reads","Multiple mapped reads","Discarded reads",
                              "Error mapped reads","Multiple mapped reads","Discarded reads",
                              "Error mapped reads","Multiple mapped reads","Discarded reads"),
                      Mean=c(mean(report1[report1$V2=="aErrorMapNum",3]/allReadsNum,report2[report2$V2=="aErrorMapNum",3]/allReadsNum,report2[report2$V2=="aErrorMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="aMultMapNum",3]/allReadsNum,report2[report2$V2=="aMultMapNum",3]/allReadsNum,report2[report2$V2=="aMultMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="aNonMapNum",3]/allReadsNum,report2[report2$V2=="aNonMapNum",3]/allReadsNum,report2[report2$V2=="aNonMapNum",3]/allReadsNum)*100,
                             
                             mean(report1[report1$V2=="eErrorMapNum",3]/allReadsNum,report2[report2$V2=="eErrorMapNum",3]/allReadsNum,report2[report2$V2=="eErrorMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="eMultMapNum",3]/allReadsNum,report2[report2$V2=="eMultMapNum",3]/allReadsNum,report2[report2$V2=="eMultMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="eNonMapNum",3]/allReadsNum,report2[report2$V2=="eNonMapNum",3]/allReadsNum,report2[report2$V2=="eNonMapNum",3]/allReadsNum)*100,
                             
                             mean(report1[report1$V2=="bErrorMapNum",3]/allReadsNum,report2[report2$V2=="bErrorMapNum",3]/allReadsNum,report2[report2$V2=="bErrorMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="bMultMapNum",3]/allReadsNum,report2[report2$V2=="bMultMapNum",3]/allReadsNum,report2[report2$V2=="bMultMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="bNonMapNum",3]/allReadsNum,report2[report2$V2=="bNonMapNum",3]/allReadsNum,report2[report2$V2=="bNonMapNum",3]/allReadsNum)*100,
                             
                             mean(report1[report1$V2=="cErrorMapNum",3]/allReadsNum,report2[report2$V2=="cErrorMapNum",3]/allReadsNum,report2[report2$V2=="cErrorMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="cMultMapNum",3]/allReadsNum,report2[report2$V2=="cMultMapNum",3]/allReadsNum,report2[report2$V2=="cMultMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="cNonMapNum",3]/allReadsNum,report2[report2$V2=="cNonMapNum",3]/allReadsNum,report2[report2$V2=="cNonMapNum",3]/allReadsNum)*100,
                             
                             mean(report1[report1$V2=="dErrorMapNum",3]/allReadsNum,report2[report2$V2=="dErrorMapNum",3]/allReadsNum,report2[report2$V2=="dErrorMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="dMultMapNum",3]/allReadsNum,report2[report2$V2=="dMultMapNum",3]/allReadsNum,report2[report2$V2=="dMultMapNum",3]/allReadsNum)*100,
                             mean(report1[report1$V2=="dNonMapNum",3]/allReadsNum,report2[report2$V2=="dNonMapNum",3]/allReadsNum,report2[report2$V2=="dNonMapNum",3]/allReadsNum)*100),
                      
                      Sd=c(sd(c(report1[report1$V2=="aErrorMapNum",3]/allReadsNum,report2[report2$V2=="aErrorMapNum",3]/allReadsNum,report2[report2$V2=="aErrorMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="aMultMapNum",3]/allReadsNum,report2[report2$V2=="aMultMapNum",3]/allReadsNum,report2[report2$V2=="aMultMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="aNonMapNum",3]/allReadsNum,report2[report2$V2=="aNonMapNum",3]/allReadsNum,report2[report2$V2=="aNonMapNum",3]/allReadsNum)),
                           
                           sd(c(report1[report1$V2=="eErrorMapNum",3]/allReadsNum,report2[report2$V2=="eErrorMapNum",3]/allReadsNum,report2[report2$V2=="eErrorMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="eMultMapNum",3]/allReadsNum,report2[report2$V2=="eMultMapNum",3]/allReadsNum,report2[report2$V2=="eMultMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="eNonMapNum",3]/allReadsNum,report2[report2$V2=="eNonMapNum",3]/allReadsNum,report2[report2$V2=="eNonMapNum",3]/allReadsNum)),
                           
                           sd(c(report1[report1$V2=="bErrorMapNum",3]/allReadsNum,report2[report2$V2=="bErrorMapNum",3]/allReadsNum,report2[report2$V2=="bErrorMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="bMultMapNum",3]/allReadsNum,report2[report2$V2=="bMultMapNum",3]/allReadsNum,report2[report2$V2=="bMultMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="bNonMapNum",3]/allReadsNum,report2[report2$V2=="bNonMapNum",3]/allReadsNum,report2[report2$V2=="bNonMapNum",3]/allReadsNum)),
                           
                           sd(c(report1[report1$V2=="cErrorMapNum",3]/allReadsNum,report2[report2$V2=="cErrorMapNum",3]/allReadsNum,report2[report2$V2=="cErrorMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="cMultMapNum",3]/allReadsNum,report2[report2$V2=="cMultMapNum",3]/allReadsNum,report2[report2$V2=="cMultMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="cNonMapNum",3]/allReadsNum,report2[report2$V2=="cNonMapNum",3]/allReadsNum,report2[report2$V2=="cNonMapNum",3]/allReadsNum)),
                           
                           sd(c(report1[report1$V2=="dErrorMapNum",3]/allReadsNum,report2[report2$V2=="dErrorMapNum",3]/allReadsNum,report2[report2$V2=="dErrorMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="dMultMapNum",3]/allReadsNum,report2[report2$V2=="dMultMapNum",3]/allReadsNum,report2[report2$V2=="dMultMapNum",3]/allReadsNum)),
                           sd(c(report1[report1$V2=="dNonMapNum",3]/allReadsNum,report2[report2$V2=="dNonMapNum",3]/allReadsNum,report2[report2$V2=="dNonMapNum",3]/allReadsNum))))
    res$SpeMap <- paste(species2,": ",res$Mapper,sep = "")
    res$Class <- factor(res$Class,levels = c("Discarded reads","Multiple mapped reads","Error mapped reads"))
    res$Species <- species2
    return(res)
  }
  human<- speciesRes("human","Human",error,106976204)
  cattle<- speciesRes("cattle","Cattle",error,89014078)
  pig<- speciesRes("pig","Pig",error,83397080)
  all <- rbind(pig,cattle,human)
  all$Mapper <- factor(all$Mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","BSBolt","Bwa-meth"))
  all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$Mapper)]))
  return(all)
}
err0<- errRes("0")
err0$`Sequencing error rate` <- "Sequencing error rate: 0"
err1<- errRes("1")
err1$`Sequencing error rate` <- "Sequencing error rate: 1.00%"
allErr <- rbind(err0,err1)

allErr$SpeMap <- factor(allErr$SpeMap,levels =c("Human: Bwa-meth","Human: BSBolt","Human: BSMAP","Human: Bismark-bwt2-e2e","Human: Walt",
                                                "Cattle: Bwa-meth","Cattle: BSBolt","Cattle: BSMAP","Cattle: Bismark-bwt2-e2e","Cattle: Walt",
                                                "Pig: Bwa-meth","Pig: BSBolt","Pig: BSMAP","Pig: Bismark-bwt2-e2e","Pig: Walt"))
tiff("depth5-figure/err01.tiff",width = 3000,height =1700,res = 300)
ggplot(allErr, aes(SpeMap, weight = Mean, fill =Class)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_manual(values = c( "darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,15),breaks=seq(0,15,3)) +
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title = element_text(size=13,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",face = "bold") , 
        legend.text=element_text(colour ="black",face = "bold"),
        axis.text.y = element_text(size=13,color = "black",face = "bold"),
        axis.text.x = element_text(size=12,angle = 62,hjust = 1,colour ="black",face = "bold"),
        strip.text = element_text(colour ="black",face = "bold",size=11))+
  facet_wrap( .~ `Sequencing error rate`,ncol = 2)
dev.off()


################### Figure 4b: unsatisfactory mapped reads in repeat #########################################

repErr<- function(error){
  repeatFactor <- function(species,species2,error){
    meanFunc <- function(species,className,className2,error){
      allReport1 <- read.csv(paste("depth5Err",error,"\\",species,"1\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      allReport2 <- read.csv(paste("depth5Err",error,"\\",species,"2\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      allReport3 <- read.csv(paste("depth5Err",error,"\\",species,"3\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      
      commonReport1 <- read.csv(paste("depth5Err",error,"\\",species,"1\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      commonReport2 <- read.csv(paste("depth5Err",error,"\\",species,"2\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      commonReport3 <- read.csv(paste("depth5Err",error,"\\",species,"3\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      
      data1 <- data.frame(mapper=c(className2),
                          region=c("Repetitive region","Non-repetitiven region"),
                          mean=c(mean(c(commonReport1[which(commonReport1$V1=="repeat"),2]/allReport1[which(allReport1$V1=="repeat"),2],
                                        commonReport2[which(commonReport2$V1=="repeat"),2]/allReport2[which(allReport2$V1=="repeat"),2],
                                        commonReport3[which(commonReport3$V1=="repeat"),2]/allReport3[which(allReport3$V1=="repeat"),2])),
                                 mean(c(commonReport1[which(commonReport1$V1=="nonRepeat"),2]/allReport1[which(allReport1$V1=="nonRepeat"),2],
                                        commonReport2[which(commonReport2$V1=="nonRepeat"),2]/allReport2[which(allReport2$V1=="nonRepeat"),2],
                                        commonReport3[which(commonReport3$V1=="nonRepeat"),2]/allReport3[which(allReport3$V1=="nonRepeat"),2]))),
                          sd=c(sd(c(commonReport1[which(commonReport1$V1=="repeat"),2]/allReport1[which(allReport1$V1=="repeat"),2],
                                    commonReport2[which(commonReport2$V1=="repeat"),2]/allReport2[which(allReport2$V1=="repeat"),2],
                                    commonReport3[which(commonReport3$V1=="repeat"),2]/allReport3[which(allReport3$V1=="repeat"),2])),
                               sd(c(commonReport1[which(commonReport1$V1=="nonRepeat"),2]/allReport1[which(allReport1$V1=="nonRepeat"),2],
                                    commonReport2[which(commonReport2$V1=="nonRepeat"),2]/allReport2[which(allReport2$V1=="nonRepeat"),2],
                                    commonReport3[which(commonReport3$V1=="nonRepeat"),2]/allReport3[which(allReport3$V1=="nonRepeat"),2]))))
      data1$mean <- data1$mean*100
      return(data1)
    }
    
    waltData <- meanFunc(species,"walt","Walt",error)
    bsmapData <- meanFunc(species,"bsmap","BSMAP",error)
    bismarkbwt2Data <- meanFunc(species,"bismarkbwt2","Bismark-bwt2-e2e",error)
    bsboltData <- meanFunc(species,"bsbolt","BSBolt",error)
    bwamethData <- meanFunc(species,"bwameth","Bwa-meth",error)
    res <- rbind(waltData,bsmapData,bismarkbwt2Data,bsboltData,bwamethData)
    res$region <- factor(res$region,levels = c("Repetitive region","Non-repetitiven region"))
    res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
    return(res)
  }
  
  human <- repeatFactor("human","Human",error)
  pig <- repeatFactor("pig","Pig",error)
  cattle <- repeatFactor("cattle","Cattle",error)
  all <- rbind(pig,cattle, human)
  names(all)[2] <- "Region"
  all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","BSBolt","Bwa-meth"))
  all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
  return(all)  
}
err0<- repErr("0")
err0$`Sequencing error rate` <- "Sequencing error rate: 0"
err1<- repErr("1")
err1$`Sequencing error rate` <- "Sequencing error rate: 1.00%"
allErr <- rbind(err0,err1)

allErr$SpeMap <- factor(allErr$SpeMap,levels =c("Human: Bwa-meth","Human: BSBolt","Human: BSMAP","Human: Bismark-bwt2-e2e","Human: Walt",
                                                "Cattle: Bwa-meth","Cattle: BSBolt","Cattle: BSMAP","Cattle: Bismark-bwt2-e2e","Cattle: Walt",
                                                "Pig: Bwa-meth","Pig: BSBolt","Pig: BSMAP","Pig: Bismark-bwt2-e2e","Pig: Walt"))
tiff("depth5-figure/repErr01.tiff",width = 3000,height =1700,res = 300)
ggplot(allErr, aes(SpeMap, weight = mean, fill =Region)) +
  geom_bar(width = .9, position = 'dodge') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_manual(values = c( "tan2","peachpuff2"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,16),breaks=seq(0,16,4)) +
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title = element_text(size=13,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",face = "bold") , 
        legend.text=element_text(colour ="black",face = "bold"),
        axis.text.y = element_text(size=13,color = "black",face = "bold"),
        axis.text.x = element_text(size=12,angle = 62,hjust = 1,colour ="black",face = "bold"),
        strip.text = element_text(colour ="black",face = "bold",size=11))+
  facet_wrap( .~ `Sequencing error rate`,ncol = 2)
dev.off()

################### Figure 4c: unsatisfactory mapped reads in cgi #######################################
cgiErr<- function(error){
  cgiFactor <- function(species,species2,error){
    meanFunc <- function(species,className,className2,error){
      allReport1 <- read.csv(paste("depth5Err",error,"\\",species,"1\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      allReport2 <- read.csv(paste("depth5Err",error,"\\",species,"2\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      allReport3 <- read.csv(paste("depth5Err",error,"\\",species,"3\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      
      commonReport1 <- read.csv(paste("depth5Err",error,"\\",species,"1\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      commonReport2 <- read.csv(paste("depth5Err",error,"\\",species,"2\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      commonReport3 <- read.csv(paste("depth5Err",error,"\\",species,"3\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      
      data1 <- data.frame(mapper=c(className2),
                          region=c("CGI","Non-CGI"),
                          mean=c(mean(c(commonReport1[which(commonReport1$V1=="cgiNum"),2]/allReport1[which(allReport1$V1=="cgiNum"),2],
                                        commonReport2[which(commonReport2$V1=="cgiNum"),2]/allReport2[which(allReport2$V1=="cgiNum"),2],
                                        commonReport3[which(commonReport3$V1=="cgiNum"),2]/allReport3[which(allReport3$V1=="cgiNum"),2])),
                                 mean(c(commonReport1[which(commonReport1$V1=="nonCgiNum"),2]/allReport1[which(allReport1$V1=="nonCgiNum"),2],
                                        commonReport2[which(commonReport2$V1=="nonCgiNum"),2]/allReport2[which(allReport2$V1=="nonCgiNum"),2],
                                        commonReport3[which(commonReport3$V1=="nonCgiNum"),2]/allReport3[which(allReport3$V1=="nonCgiNum"),2]))),
                          sd=c(sd(c(commonReport1[which(commonReport1$V1=="cgiNum"),2]/allReport1[which(allReport1$V1=="cgiNum"),2],
                                    commonReport2[which(commonReport2$V1=="cgiNum"),2]/allReport2[which(allReport2$V1=="cgiNum"),2],
                                    commonReport3[which(commonReport3$V1=="cgiNum"),2]/allReport3[which(allReport3$V1=="cgiNum"),2])),
                               sd(c(commonReport1[which(commonReport1$V1=="nonCgiNum"),2]/allReport1[which(allReport1$V1=="nonCgiNum"),2],
                                    commonReport2[which(commonReport2$V1=="nonCgiNum"),2]/allReport2[which(allReport2$V1=="nonCgiNum"),2],
                                    commonReport3[which(commonReport3$V1=="nonCgiNum"),2]/allReport3[which(allReport3$V1=="nonCgiNum"),2]))))
      data1$mean <- data1$mean*100
      return(data1)
    }
    waltData <- meanFunc(species,"walt","Walt",error)
    bsmapData <- meanFunc(species,"bsmap","BSMAP",error)
    bismarkbwt2Data <- meanFunc(species,"bismarkbwt2","Bismark-bwt2-e2e",error)
    bsboltData <- meanFunc(species,"bsbolt","BSBolt",error)
    bwamethData <- meanFunc(species,"bwameth","Bwa-meth",error)
    res <- rbind(waltData,bsmapData,bismarkbwt2Data,bsboltData,bwamethData)
    res$region <- factor(res$region,levels = c("CGI","Non-CGI"))
    res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
    return(res)
  }
  human <- cgiFactor("human","Human",error)
  pig <- cgiFactor("pig","Pig",error)
  cattle <- cgiFactor("cattle","Cattle",error)
  all <- rbind(pig,cattle, human)
  names(all)[2] <- "Region"
  all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","BSBolt","Bwa-meth"))
  all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
  return(all)
}
cgiErr0<- cgiErr("0")
cgiErr0$`Sequencing error rate` <- "Sequencing error rate: 0"
cgiErr1<- cgiErr("1")
cgiErr1$`Sequencing error rate` <- "Sequencing error rate: 1.00%"
allCgiErr <- rbind(cgiErr0,cgiErr1)


allCgiErr$SpeMap <- factor(allCgiErr$SpeMap,levels  =c("Human: Bwa-meth","Human: BSBolt","Human: BSMAP","Human: Bismark-bwt2-e2e","Human: Walt",
                                                       "Cattle: Bwa-meth","Cattle: BSBolt","Cattle: BSMAP","Cattle: Bismark-bwt2-e2e","Cattle: Walt",
                                                       "Pig: Bwa-meth","Pig: BSBolt","Pig: BSMAP","Pig: Bismark-bwt2-e2e","Pig: Walt"))
tiff("depth5-figure/cgiErr01.tiff",width = 3000,height =1700,res = 300)
ggplot(allCgiErr, aes(SpeMap, weight = mean, fill =Region)) +
  geom_bar(width = .9, position = 'dodge') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_manual(values = c( "paleturquoise4","azure3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,24),breaks=seq(0,24,4)) +
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title = element_text(size=13,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",face = "bold") , 
        legend.text=element_text(colour ="black",face = "bold"),
        axis.text.y = element_text(size=13,color = "black",face = "bold"),
        axis.text.x = element_text(size=12,angle = 62,hjust = 1,colour ="black",face = "bold"),
        strip.text = element_text(colour ="black",face = "bold",size=11))+
  facet_wrap( .~ `Sequencing error rate`,ncol = 2)
dev.off()
