setwd("..\\result\\depth5")
library(ggplot2)
library(RColorBrewer)
################### Error rate (0)£ºthe class of incorrect mapped reads ################
mapClassOfIncorrectMapRead <- function(species,singleReadsNum,species2){
  resRbind <- function(species,mapperName,mapperName2,singleReadsNum){
    report1 <- read.csv(paste(species,"1\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    report2 <- read.csv(paste(species,"2\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    report3 <- read.csv(paste(species,"3\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
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
  
  bismarkbwt2Data<- resRbind(species,"bismarkbwt2","Bismark-bwt2-e2e",singleReadsNum)
  bsmapData<- resRbind(species,"bsmap","BSMAP",singleReadsNum)
  bwamethData<- resRbind(species,"bwameth","Bwa-meth",singleReadsNum)
  waltData<- resRbind(species,"walt","Walt",singleReadsNum)
  res <- rbind(bismarkbwt2Data,bsmapData,bwamethData,waltData)
  res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
  res$class <- factor(res$class,levels = c("Discarded reads","Multiple mapped reads","Error mapped reads"))
  return(res)
  
}
human <- mapClassOfIncorrectMapRead("human",106976204,"Human")
mouse <- mapClassOfIncorrectMapRead("mouse",91029060,"Mouse")
cattle <- mapClassOfIncorrectMapRead("cattle",89014078,"Cattle")
pig <- mapClassOfIncorrectMapRead("pig",83397080,"Pig")
all <- rbind(pig,cattle, mouse, human)
names(all)[2] <- "Class"
all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))

tiff("firstPart\\second\\unsatisfactoryV2.tiff",width = 1200,height =740,res = 300)
ggplot(all, aes(SpeMap, weight = mean, fill =Class)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,15),breaks=seq(0,15,3)) +
  theme(axis.text.y = element_text(colour =c(rep("orchid4",4),rep("darkseagreen4",4),rep("lightskyblue4",4),rep("palevioletred4",4)),size = 6,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=6,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",size = 6,face = "bold") , legend.text=element_text(colour ="black",size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,colour ="black",face = "bold"))
dev.off()


################### Error rate (0)£ºthe unsatisfactory aligned reads in repetitive sequence ####################
repeatFactor <- function(species,species2){
  meanFunc <- function(species,className,className2){
    allReport1 <- read.csv(paste(species,"1\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    allReport2 <- read.csv(paste(species,"2\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    allReport3 <- read.csv(paste(species,"3\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    
    commonReport1 <- read.csv(paste(species,"1\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    commonReport2 <- read.csv(paste(species,"2\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    commonReport3 <- read.csv(paste(species,"3\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    
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

  waltData <- meanFunc(species,"walt","Walt")
  bsmapData <- meanFunc(species,"bsmap","BSMAP")
  bismarkbwt2Data <- meanFunc(species,"bismarkbwt2","Bismark-bwt2-e2e")
  bwamethData <- meanFunc(species,"bwameth","Bwa-meth")
  res <- rbind(waltData,bsmapData,bismarkbwt2Data,bwamethData)
  res$region <- factor(res$region,levels = c("Repetitive region","Non-repetitiven region"))
  res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
  return(res)
}

human <- repeatFactor("human","Human")
pig <- repeatFactor("pig","Pig")
mouse <- repeatFactor("mouse","Mouse")
cattle <- repeatFactor("cattle","Cattle")
all <- rbind(pig,cattle, mouse, human)
names(all)[2] <- "Region"
### v2
all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("firstPart\\second\\repAndNonRepV2.tiff",width = 1200,height =740,res = 300)
ggplot(all, aes(SpeMap, weight = mean, fill =Region)) +
  geom_bar(width = .9, position = 'dodge') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "tan2","peachpuff2"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,18),breaks=seq(0,18,3)) +
  theme(axis.text.y = element_text(colour =c(rep("orchid4",4),rep("darkseagreen4",4),rep("lightskyblue4",4),rep("palevioletred4",4)),size = 6,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=6,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",size = 6,face = "bold") , legend.text=element_text(colour ="black",size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,colour ="black",face = "bold"))
dev.off()
################### Error rate (0)£ºthe unsatisfactory aligned reads in CGI #####################
cgiFactor <- function(species,species2){
  meanFunc <- function(species,className,className2){
    allReport1 <- read.csv(paste(species,"1\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    allReport2 <- read.csv(paste(species,"2\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    allReport3 <- read.csv(paste(species,"3\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    
    commonReport1 <- read.csv(paste(species,"1\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    commonReport2 <- read.csv(paste(species,"2\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    commonReport3 <- read.csv(paste(species,"3\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    
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
  waltData <- meanFunc(species,"walt","Walt")
  bsmapData <- meanFunc(species,"bsmap","BSMAP")
  bismarkbwt2Data <- meanFunc(species,"bismarkbwt2","Bismark-bwt2-e2e")
  bwamethData <- meanFunc(species,"bwameth","Bwa-meth")
  res <- rbind(waltData,bsmapData,bismarkbwt2Data,bwamethData)
  res$region <- factor(res$region,levels = c("CGI","Non-CGI"))
  res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
  return(res)
}
human <- cgiFactor("human","Human")
pig <- cgiFactor("pig","Pig")
mouse <- cgiFactor("mouse","Mouse")
cattle <- cgiFactor("cattle","Cattle")
all <- rbind(pig,cattle, mouse, human)
names(all)[2] <- "Region"
all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("firstPart\\second\\cgiAndNonCgiV2.tiff",width = 1100,height =740,res = 300)
ggplot(all, aes(SpeMap, weight = mean, fill =Region)) +
  geom_bar(width = .9, position = 'dodge') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise4","azure3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,24),breaks=seq(0,24,4)) +
  theme(axis.text.y = element_text(colour =c(rep("orchid4",4),rep("darkseagreen4",4),rep("lightskyblue4",4),rep("palevioletred4",4)),size = 6,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=6,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",size = 6,face = "bold") , legend.text=element_text(colour ="black",size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,colour ="black",face = "bold"))
dev.off()

################### Error rate (1.00%): the class of incorrect mapped reads ####################
mapClassOfIncorrectMapRead <- function(species,singleReadsNum,species2){
  resRbind <- function(species,mapperName,mapperName2,singleReadsNum){
    report1 <- read.csv(paste(species,"1\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    report2 <- read.csv(paste(species,"2\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    report3 <- read.csv(paste(species,"3\\readId\\",mapperName,"NonRightMapReadMappingClassReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
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
  
  bismarkbwt2Data<- resRbind(species,"bismarkbwt2","Bismark-bwt2-e2e",singleReadsNum)
  bsmapData<- resRbind(species,"bsmap","BSMAP",singleReadsNum)
  bwamethData<- resRbind(species,"bwameth","Bwa-meth",singleReadsNum)
  waltData<- resRbind(species,"walt","Walt",singleReadsNum)
  res <- rbind(bismarkbwt2Data,bsmapData,bwamethData,waltData)
  res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
  res$class <- factor(res$class,levels = c("Discarded reads","Multiple mapped reads","Error mapped reads"))
  return(res)
  
}
human <- mapClassOfIncorrectMapRead("human",106976204,"Human")
mouse <- mapClassOfIncorrectMapRead("mouse",91029060,"Mouse")
cattle <- mapClassOfIncorrectMapRead("cattle",89014078,"Cattle")
pig <- mapClassOfIncorrectMapRead("pig",83397080,"Pig")
all <- rbind(pig,cattle, mouse, human)
names(all)[2] <- "Class"

all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))

tiff("figure\\unMapClassV2.tiff",width = 1200,height =740,res = 300)
ggplot(all, aes(SpeMap, weight = mean, fill =Class)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,15),breaks=seq(0,15,3)) +
  theme(axis.text.y = element_text(colour =c(rep("orchid4",4),rep("darkseagreen4",4),rep("lightskyblue4",4),rep("palevioletred4",4)),size = 6,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=6,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",size = 6,face = "bold") , legend.text=element_text(colour ="black",size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,colour ="black",face = "bold"))
dev.off()

################### Error rate (1.00%): the unsatisfactory aligned reads in repetitive sequence ############################
repeatFactor <- function(species,species2){
  meanFunc <- function(species,className,className2){
    allReport1 <- read.csv(paste(species,"1\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    allReport2 <- read.csv(paste(species,"2\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    allReport3 <- read.csv(paste(species,"3\\readAna\\allReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    
    commonReport1 <- read.csv(paste(species,"1\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    commonReport2 <- read.csv(paste(species,"2\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    commonReport3 <- read.csv(paste(species,"3\\readAna\\",className,"NonRightMapReadIdRepeatAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    
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
  
  waltData <- meanFunc(species,"walt","Walt")
  bsmapData <- meanFunc(species,"bsmap","BSMAP")
  bismarkbwt2Data <- meanFunc(species,"bismarkbwt2","Bismark-bwt2-e2e")
  bwamethData <- meanFunc(species,"bwameth","Bwa-meth")
  res <- rbind(waltData,bsmapData,bismarkbwt2Data,bwamethData)
  res$region <- factor(res$region,levels = c("Repetitive region","Non-repetitiven region"))
  res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
  return(res)
}

human <- repeatFactor("human","Human")
pig <- repeatFactor("pig","Pig")
mouse <- repeatFactor("mouse","Mouse")
cattle <- repeatFactor("cattle","Cattle")
all <- rbind(pig,cattle, mouse, human)
names(all)[2] <- "Region"
all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\repAndNonRepV2.tiff",width = 1200,height =740,res = 300)
ggplot(all, aes(SpeMap, weight = mean, fill =Region)) +
  geom_bar(width = .9, position = 'dodge') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "tan2","peachpuff2"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,18),breaks=seq(0,18,3)) +
  theme(axis.text.y = element_text(colour =c(rep("orchid4",4),rep("darkseagreen4",4),rep("lightskyblue4",4),rep("palevioletred4",4)),size = 6,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=6,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",size = 6,face = "bold") , legend.text=element_text(colour ="black",size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,colour ="black",face = "bold"))
dev.off()

################### Error rate (1.00%): the unsatisfactory aligned reads in CGI ########################
cgiFactor <- function(species,species2){
  meanFunc <- function(species,className,className2){
    allReport1 <- read.csv(paste(species,"1\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    allReport2 <- read.csv(paste(species,"2\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    allReport3 <- read.csv(paste(species,"3\\readAna\\allReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    
    commonReport1 <- read.csv(paste(species,"1\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    commonReport2 <- read.csv(paste(species,"2\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    commonReport3 <- read.csv(paste(species,"3\\readAna\\",className,"NonRightMapReadIdCgiAnnoResStaReport.txt",sep = ""),header = F, sep = "\t", quote = "\"", dec = ".")
    
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
  waltData <- meanFunc(species,"walt","Walt")
  bsmapData <- meanFunc(species,"bsmap","BSMAP")
  bismarkbwt2Data <- meanFunc(species,"bismarkbwt2","Bismark-bwt2-e2e")
  bwamethData <- meanFunc(species,"bwameth","Bwa-meth")
  res <- rbind(waltData,bsmapData,bismarkbwt2Data,bwamethData)
  res$region <- factor(res$region,levels = c("CGI","Non-CGI"))
  res$SpeMap <- paste(species2,": ",res$mapper,sep = "")
  return(res)
}
human <- cgiFactor("human","Human")
pig <- cgiFactor("pig","Pig")
mouse <- cgiFactor("mouse","Mouse")
cattle <- cgiFactor("cattle","Cattle")
all <- rbind(pig,cattle, mouse, human)
names(all)[2] <- "Region"
all$mapper <- factor(all$mapper,levels = c("Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\cgiAndNonCgiV2.tiff",width = 1100,height =740,res = 300)
ggplot(all, aes(SpeMap, weight = mean, fill =Region)) +
  geom_bar(width = .9, position = 'dodge') +
  labs( x="",y = 'Unsatisfactory aligned reads (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise4","azure3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,24),breaks=seq(0,24,4)) +
  theme(axis.text.y = element_text(colour =c(rep("orchid4",4),rep("darkseagreen4",4),rep("lightskyblue4",4),rep("palevioletred4",4)),size = 6,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=6,face = "bold"),
        legend.key.size = unit(12, "pt"),
        legend.title=element_text(colour ="black",size = 6,face = "bold") , legend.text=element_text(colour ="black",size = 6,face = "bold"),
        axis.text.x = element_text(size = 6,colour ="black",face = "bold"))
dev.off()

