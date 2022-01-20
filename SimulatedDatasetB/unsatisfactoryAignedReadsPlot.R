setwd("..\\result\\depth5")
library(ggplot2)
library(RColorBrewer)
library(ggplot2)
library(RColorBrewer)
################### The class of unsatisfactory mapped reads #######################

sinError<- function(error){
  mapClassOfIncorrectMapRead <- function(species,singleReadsNum,species2,error){
    resRbind <- function(species,mapperName,mapperName2,singleReadsNum,error){
      report1 <- read.csv(paste("depth5Error",error,"\\",species,"1\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      report2 <- read.csv(paste("depth5Error",error,"\\",species,"2\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      report3 <- read.csv(paste("depth5Error",error,"\\",species,"3\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      data1 <- data.frame(mapper=character(),class=character(),value=numeric())
      data1 <- rbind(data1,data.frame(mapper=mapperName2,class="Error mapped reads",
                                      mean=c(mean(c(report1[1,3]/singleReadsNum,report2[1,3]/singleReadsNum,report3[1,3]/singleReadsNum))*100),
                                      sd=c(sd(c(report1[1,3]/singleReadsNum,report2[1,3]/singleReadsNum,report3[1,3]/singleReadsNum)))),
                     data.frame(mapper=mapperName2,class="Multiple mapped reads",
                                mean=c(mean(c(report1[2,3]/singleReadsNum,report2[2,3]/singleReadsNum,report3[2,3]/singleReadsNum))*100),
                                sd=c(sd(c(report1[2,3]/singleReadsNum,report2[2,3]/singleReadsNum,report3[2,3]/singleReadsNum)))),
                     data.frame(mapper=mapperName2,class="Discarded reads",
                                mean=c(mean(c(report1[3,3]/singleReadsNum,report2[3,3]/singleReadsNum,report3[3,3]/singleReadsNum))*100),
                                sd=c(sd(c(report1[3,3]/singleReadsNum,report2[3,3]/singleReadsNum,report3[3,3]/singleReadsNum)))))
      return(data1)
    }
    
    bismarkbwt2Data<- resRbind(species,"bismarkbwt2","Bismark-bwt2-e2e",singleReadsNum,error)
    bsmapData<- resRbind(species,"bsmap","BSMAP",singleReadsNum,error)
    bwamethData<- resRbind(species,"bwameth","Bwa-meth",singleReadsNum,error)
    waltData<- resRbind(species,"walt","Walt",singleReadsNum,error)
    res <- rbind(bismarkbwt2Data,bsmapData,bwamethData,waltData)
    res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
    res$class <- factor(res$class,levels = c("Discarded reads","Multiple mapped reads","Error mapped reads"))
    res$Species <- species2
    return(res)
  }
  human <- mapClassOfIncorrectMapRead("human",106976204,"Human",error)
  cattle <- mapClassOfIncorrectMapRead("cattle",89014078,"Cattle",error)
  pig <- mapClassOfIncorrectMapRead("pig",83397080,"Pig",error)
  all <- rbind(pig,cattle,human)
  names(all)[2] <- "Class"
  all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
  all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
  return(all)
}

err0<- sinError("0")
err0$`Sequencing error rate` <- "Sequencing error rate: 0"
err1<- sinError("1")
err1$`Sequencing error rate` <- "Sequencing error rate: 1.00%"
allErr <- rbind(err0,err1)


allErr$SpeMap <- factor(allErr$SpeMap,levels =c("Human: Bwa-meth","Human: BSMAP","Human: Bismark-bwt2-e2e","Human: Walt",
                                                "Cattle: Bwa-meth","Cattle: BSMAP","Cattle: Bismark-bwt2-e2e","Cattle: Walt",
                                                "Pig: Bwa-meth","Pig: BSMAP","Pig: Bismark-bwt2-e2e","Pig: Walt"))
tiff("depth5-figure/err01.tiff",width = 2700,height =1700,res = 300)
ggplot(allErr, aes(SpeMap, weight = mean, fill =Class)) +
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


################### unsatisfactory mapped reads in repeat #########################################
repErr<- function(error){
  repeatFactor <- function(species,species2,error){
    meanFunc <- function(species,className,className2,error){
      allReport1 <- read.csv(paste("depth5Error",error,"\\",species,"1\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      allReport2 <- read.csv(paste("depth5Error",error,"\\",species,"2\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      allReport3 <- read.csv(paste("depth5Error",error,"\\",species,"3\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      
      commonReport1 <- read.csv(paste("depth5Error",error,"\\",species,"1\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      commonReport2 <- read.csv(paste("depth5Error",error,"\\",species,"2\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      commonReport3 <- read.csv(paste("depth5Error",error,"\\",species,"3\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      
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
    bwamethData <- meanFunc(species,"bwameth","Bwa-meth",error)
    res <- rbind(waltData,bsmapData,bismarkbwt2Data,bwamethData)
    res$region <- factor(res$region,levels = c("Repetitive region","Non-repetitiven region"))
    res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
    return(res)
  }
  
  human <- repeatFactor("human","Human",error)
  pig <- repeatFactor("pig","Pig",error)
  cattle <- repeatFactor("cattle","Cattle",error)
  all <- rbind(pig,cattle, human)
  names(all)[2] <- "Region"
  all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
  all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
  return(all)  
}
err0<- repErr("0")
err0$`Sequencing error rate` <- "Sequencing error rate: 0"
err1<- repErr("1")
err1$`Sequencing error rate` <- "Sequencing error rate: 1.00%"
allErr <- rbind(err0,err1)


allErr$SpeMap <- factor(allErr$SpeMap,levels =c("Human: Bwa-meth","Human: BSMAP","Human: Bismark-bwt2-e2e","Human: Walt",
                                                "Cattle: Bwa-meth","Cattle: BSMAP","Cattle: Bismark-bwt2-e2e","Cattle: Walt",
                                                "Pig: Bwa-meth","Pig: BSMAP","Pig: Bismark-bwt2-e2e","Pig: Walt"))
tiff("depth5-figure/repErr01.tiff",width = 2700,height =1700,res = 300)
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

################### unsatisfactory mapped reads in cgi #######################################
cgiErr<- function(error){
  cgiFactor <- function(species,species2,error){
    meanFunc <- function(species,className,className2,error){
      allReport1 <- read.csv(paste("depth5Error",error,"\\",species,"1\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      allReport2 <- read.csv(paste("depth5Error",error,"\\",species,"2\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      allReport3 <- read.csv(paste("depth5Error",error,"\\",species,"3\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      
      commonReport1 <- read.csv(paste("depth5Error",error,"\\",species,"1\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      commonReport2 <- read.csv(paste("depth5Error",error,"\\",species,"2\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      commonReport3 <- read.csv(paste("depth5Error",error,"\\",species,"3\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
      
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
    bwamethData <- meanFunc(species,"bwameth","Bwa-meth",error)
    res <- rbind(waltData,bsmapData,bismarkbwt2Data,bwamethData)
    res$region <- factor(res$region,levels = c("CGI","Non-CGI"))
    res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
    return(res)
  }
  human <- cgiFactor("human","Human",error)
  pig <- cgiFactor("pig","Pig",error)
  cattle <- cgiFactor("cattle","Cattle",error)
  all <- rbind(pig,cattle, human)
  names(all)[2] <- "Region"
  all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
  all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
  return(all)
}
cgiErr0<- cgiErr("0")
cgiErr0$`Sequencing error rate` <- "Sequencing error rate: 0"
cgiErr1<- cgiErr("1")
cgiErr1$`Sequencing error rate` <- "Sequencing error rate: 1.00%"
allCgiErr <- rbind(cgiErr0,cgiErr1)

allCgiErr$SpeMap <- factor(allCgiErr$SpeMap,levels =c("Human: Bwa-meth","Human: BSMAP","Human: Bismark-bwt2-e2e","Human: Walt",
                                                      "Cattle: Bwa-meth","Cattle: BSMAP","Cattle: Bismark-bwt2-e2e","Cattle: Walt",
                                                      "Pig: Bwa-meth","Pig: BSMAP","Pig: Bismark-bwt2-e2e","Pig: Walt"))
tiff("depth5-figure/cgiErr01.tiff",width = 2700,height =1700,res = 300)
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
