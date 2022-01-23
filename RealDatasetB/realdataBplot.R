setwd("..\\WGBS\\result\\realdata")
library("ggplot2")
library(stringr)
library(VennDiagram)
############## CpG: The total number of CpG sites ################################# 
countToalNumCpG <- function(species,sampleName1,sampleName2,sampleName3,sampleName4,sampleName5,sampleName6,depth){
  sampleReport1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName1,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName2,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName3,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName4,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName5,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName6,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  
  me<- function(sampleReport1,sampleReport2,sampleReport3,sampleReport4,sampleReport5,sampleReport6){
    res <- c()
    for(i in 5:19){
      temp<- mean(c(as.numeric(as.character(sampleReport1[i,2])),
                    as.numeric(as.character(sampleReport2[i,2])),
                    as.numeric(as.character(sampleReport3[i,2])),
                    as.numeric(as.character(sampleReport4[i,2])),
                    as.numeric(as.character(sampleReport5[i,2])),
                    as.numeric(as.character(sampleReport6[i,2]))))
      
      res <- c(res,temp)
    }
    return(res)
  }
  
  classStr <- function(sampleReport1){
    res <- c()
    for(i in 5:19){
      temp <- str_sub(as.character(sampleReport1[i,1]),end = -4)
      res <- c(res,temp)
    }
    return(res)
  }
  
  res<- data.frame(class=classStr(sampleReport1),
                   value=c(me(sampleReport1,sampleReport2,sampleReport3,sampleReport4,sampleReport5,sampleReport6)))
  value <- function(class){
    temp<- res[which(res$class==class),2]
    temp <- as.numeric(temp)
    return(temp)
  }
  sumTotal <- value("a")+value("b")+value("c")+value("d")+value("abc")+value("abd")+value("acd")+value("bcd")-
    value("ab")-value("ac")-value("ad")-value("bc")-value("bd")-value("cd")-value("abcd")
  t <- data.frame(class="totalVenn",
                  value=sumTotal)
  
  res <- rbind(res,t)
  res$rate <- res$value/res$value[16]
  return(res)
}

human10 <- countToalNumCpG("human","SRR6373932","SRR6818517","SRR6373926","SRR6373923","SRR6825466","SRR6825471","10")
cattle10 <- countToalNumCpG("cattle","SRR7528450","SRR7528456","SRR7528458","SRR7528459","SRR7528464","SRR7528456","10")
pig10 <- countToalNumCpG("pig","SRR7812176","SRR7812178","SRR7812179","SRR7812199","SRR7812200","SRR7812210","10")

humanRes <- human10[1:4,]
humanRes$Mapper <- c("Bismark-bwt2-e2e","BSMAP","Bwa-meth","Walt")
humanRes$Species <- "Human"
cattleRes <- cattle10[1:4,]
cattleRes$Mapper <- c("Bismark-bwt2-e2e","BSMAP","Bwa-meth","Walt")
cattleRes$Species <- "Cattle"
pigRes <- pig10[1:4,]
pigRes$Mapper <- c("Bismark-bwt2-e2e","BSMAP","Bwa-meth","Walt")
pigRes$Species <- "Pig"

allRes <- rbind(humanRes,cattleRes,pigRes)
rm(human10,cattle10,pig10) 
allRes$Mapper <- factor(allRes$Mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/cpg/numOfCpGMapper.tiff",width = 2200,height =1540,res = 300)
ggplot(allRes, aes(Mapper, weight = value)) +
  geom_bar(color = "black",fill='pink3', width = .8, position = 'stack') +
  labs( x="",y = 'The number of CpG sties')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 47, hjust = 1,colour ="black",face = "bold",size = 13),
        axis.text.y = element_text(colour ="black",face = "bold",size = 13),
        axis.title.y = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        axis.line = element_line(color="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()

############## CpG: The number of shared CpG sites ######################################
countToalNumCpG <- function(species,sampleName1,sampleName2,sampleName3,sampleName4,sampleName5,sampleName6,depth){
  sampleReport1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName1,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName2,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName3,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName4,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName5,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName6,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  
  me<- function(sampleReport1,sampleReport2,sampleReport3,sampleReport4,sampleReport5,sampleReport6){
    res <- c()
    for(i in 5:19){
      temp<- mean(c(as.numeric(as.character(sampleReport1[i,2])),
                    as.numeric(as.character(sampleReport2[i,2])),
                    as.numeric(as.character(sampleReport3[i,2])),
                    as.numeric(as.character(sampleReport4[i,2])),
                    as.numeric(as.character(sampleReport5[i,2])),
                    as.numeric(as.character(sampleReport6[i,2]))))
      
      res <- c(res,temp)
    }
    return(res)
  }
  
  classStr <- function(sampleReport1){
    res <- c()
    for(i in 5:19){
      temp <- str_sub(as.character(sampleReport1[i,1]),end = -4)
      res <- c(res,temp)
    }
    return(res)
  }
  
  res<- data.frame(class=classStr(sampleReport1),
                   value=c(me(sampleReport1,sampleReport2,sampleReport3,sampleReport4,sampleReport5,sampleReport6)))
  value <- function(class){
    temp<- res[which(res$class==class),2]
    temp <- as.numeric(temp)
    return(temp)
  }
  sumTotal <- value("a")+value("b")+value("c")+value("d")+value("abc")+value("abd")+value("acd")+value("bcd")-
    value("ab")-value("ac")-value("ad")-value("bc")-value("bd")-value("cd")-value("abcd")
  t <- data.frame(class="totalVenn",
                  value=sumTotal)
  
  res <- rbind(res,t)
  res$rate <- res$value/res$value[16]
  return(res)
}

human10 <- countToalNumCpG("human","SRR6373932","SRR6818517","SRR6373926","SRR6373923","SRR6825466","SRR6825471","10")
cattle10 <- countToalNumCpG("cattle","SRR7528450","SRR7528456","SRR7528458","SRR7528459","SRR7528464","SRR7528456","10")
pig10 <- countToalNumCpG("pig","SRR7812176","SRR7812178","SRR7812179","SRR7812199","SRR7812200","SRR7812210","10")

singleSpe<- function(data,species){
  tempRes <- data[11:15,]
  tempRes$Mapper <- c("Bwa,BSM,Bis","BSM,Bis,Wal","Bwa,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis,Wal")
  tempRes$Species <- species
  tempRes$SpeMap <- c(paste0(species,": Bwa,BSM,Bis"),
                      paste0(species,": BSM,Bis,Wal"),
                      paste0(species,": Bwa,Bis,Wal"),
                      paste0(species,": Bwa,BSM,Wal"),
                      paste0(species,": Bwa,BSM,Bis,Wal"))
  return(tempRes)
}
humanRes<- singleSpe(human10,"Human")
cattleRes<- singleSpe(cattle10,"Cattle")
pigRes<- singleSpe(pig10,"Pig")
allRes <- rbind(humanRes,cattleRes,pigRes)
rm(human10,cattle10,pig10) 

allRes$Mapper <- factor(allRes$Mapper,levels = c("Bwa,BSM,Bis,Wal","Bwa,Bis,Wal","BSM,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))


tiff("figure/cpg/concordantCpG.tiff",width = 2400,height =1740,res = 300)
ggplot(allRes,aes(x=Mapper,y=value,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'The number of shared CpG sites')+
  scale_y_continuous(expand = c(0,0),limits=c(0,2.0e+07),breaks = seq(0,2.0e+07,5e+06)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") , 
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()

############## CpG：The number of reliable CpG sites #########################
countToalNumCpG <- function(species,sampleName1,sampleName2,sampleName3,sampleName4,sampleName5,sampleName6,depth){
  sampleReport1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName1,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName2,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName3,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName4,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName5,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName6,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  
  me<- function(sampleReport1,sampleReport2,sampleReport3,sampleReport4,sampleReport5,sampleReport6){
    res <- c()
    for(i in 5:19){
      temp<- mean(c(as.numeric(as.character(sampleReport1[i,2])),
                    as.numeric(as.character(sampleReport2[i,2])),
                    as.numeric(as.character(sampleReport3[i,2])),
                    as.numeric(as.character(sampleReport4[i,2])),
                    as.numeric(as.character(sampleReport5[i,2])),
                    as.numeric(as.character(sampleReport6[i,2]))))
      
      res <- c(res,temp)
    }
    return(res)
  }
  
  classStr <- function(sampleReport1){
    res <- c()
    for(i in 5:19){
      temp <- str_sub(as.character(sampleReport1[i,1]),end = -4)
      res <- c(res,temp)
    }
    return(res)
  }
  
  res<- data.frame(class=classStr(sampleReport1),
                   value=c(me(sampleReport1,sampleReport2,sampleReport3,sampleReport4,sampleReport5,sampleReport6)))
  value <- function(class){
    temp<- res[which(res$class==class),2]
    temp <- as.numeric(temp)
    return(temp)
  }
  sumTotal <- value("a")+value("b")+value("c")+value("d")+value("abc")+value("abd")+value("acd")+value("bcd")-
    value("ab")-value("ac")-value("ad")-value("bc")-value("bd")-value("cd")-value("abcd")
  t <- data.frame(class="totalVenn",
                  value=sumTotal)
  
  res <- rbind(res,t)
  res$rate <- res$value/res$value[16]
  return(res)
}

human10 <- countToalNumCpG("human","SRR6373932","SRR6818517","SRR6373926","SRR6373923","SRR6825466","SRR6825471","10")
cattle10 <- countToalNumCpG("cattle","SRR7528450","SRR7528456","SRR7528458","SRR7528459","SRR7528464","SRR7528456","10")
pig10 <- countToalNumCpG("pig","SRR7812176","SRR7812178","SRR7812179","SRR7812199","SRR7812200","SRR7812210","10")

singleSpe<- function(data,species){
  res<- data.frame(value=c(sum(data[c(11,13,14),2])-2*data[15,2], ## Bwa-meth 
                           sum(data[c(11,12,14),2])-2*data[15,2], ## BSMAP
                           sum(data[c(11,12,13),2])-2*data[15,2], ## Bismark-bwt2-e2e
                           sum(data[c(12,13,14),2])-2*data[15,2]), ## Walt
                   Rate=c((sum(data[c(11,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bwa-meth 
                          (sum(data[c(11,12,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## BSMAP
                          (sum(data[c(11,12,13),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bismark-bwt2-e2e
                          (sum(data[c(12,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2])), ## Walt
                   mapper=c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"),
                   Class=c("Covered by three alignment algorithms in Bwa-meth",
                           "Covered by three alignment algorithms in BSMAP",
                           "Covered by three alignment algorithms in Bismark-bwt2-e2e",
                           "Covered by three alignment algorithms in Walt"),
                   Species=species,
                   SpeCla=c(paste0(species,": Covered by three alignment algorithms in Bwa-meth"),
                            paste0(species,": Covered by three alignment algorithms in BSMAP"),
                            paste0(species,": Covered by three alignment algorithms in Bismark-bwt2-e2e"),
                            paste0(species,": Covered by three alignment algorithms in Walt")))
  return(res)
}
humanRes<- singleSpe(human10,"Human")
cattleRes<- singleSpe(cattle10,"Cattle")
pigRes<- singleSpe(pig10,"Pig")

allRes <- rbind(humanRes,cattleRes,pigRes)
rm(human10,cattle10,pig10) 

allRes$mapper <- factor(allRes$mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))

allRes$Class <- factor(allRes$Class,levels = c(c("Covered by three alignment algorithms in Bwa-meth",
                                                 "Covered by three alignment algorithms in BSMAP",
                                                 "Covered by three alignment algorithms in Bismark-bwt2-e2e",
                                                 "Covered by three alignment algorithms in Walt")))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/20220113/reliable_CpG_rate.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=Rate*100,group=Species))+
  labs( x="",y = 'The proportion of accurate CpG sites (%)',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  scale_y_continuous(expand = c(0,0),limits=c(86,99),breaks=seq(87,99,3)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()


tiff("figure/20220113/reliable_CpG_value.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=value,group=Species))+
  labs( x="",y = 'The number of accurate CpG sites',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()

############## CpG: The proportion of discordant and concordant CpG sites ####################################################
singleSpeciesRes <- function(species,species2, depth,sampleName1,sampleName2,sampleName3,sampleName4,sampleName5,sampleName6){
  sampleReport1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName1,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName2,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName3,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")         
  sampleReport4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName4,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName5,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName6,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")         
  res <- data.frame(Species=species,
                    Type=c("Concordant CpG","Discordant CpG"),
                    Depth=c(paste("Depth",depth,sep = "")),
                    SpeDep=c(paste(species2,": Depth",depth,sep = "")),
                    sample1=c(as.numeric(as.character(sampleReport1[20,2]))/as.numeric(as.character(sampleReport1[23,2])),
                              as.numeric(as.character(sampleReport1[22,2]))/as.numeric(as.character(sampleReport1[23,2]))),
                    sample2=c(as.numeric(as.character(sampleReport2[20,2]))/as.numeric(as.character(sampleReport2[23,2])),
                              as.numeric(as.character(sampleReport2[22,2]))/as.numeric(as.character(sampleReport2[23,2]))),
                    sample3=c(as.numeric(as.character(sampleReport3[20,2]))/as.numeric(as.character(sampleReport3[23,2])),
                              as.numeric(as.character(sampleReport3[22,2]))/as.numeric(as.character(sampleReport3[23,2]))),
                    sample4=c(as.numeric(as.character(sampleReport4[20,2]))/as.numeric(as.character(sampleReport4[23,2])),
                              as.numeric(as.character(sampleReport4[22,2]))/as.numeric(as.character(sampleReport4[23,2]))),
                    sample5=c(as.numeric(as.character(sampleReport5[20,2]))/as.numeric(as.character(sampleReport5[23,2])),
                              as.numeric(as.character(sampleReport5[22,2]))/as.numeric(as.character(sampleReport5[23,2]))),
                    sample6=c(as.numeric(as.character(sampleReport6[20,2]))/as.numeric(as.character(sampleReport6[23,2])),
                              as.numeric(as.character(sampleReport6[22,2]))/as.numeric(as.character(sampleReport6[23,2]))),
                    
                    
                    mean=c(mean(c(as.numeric(as.character(sampleReport1[20,2]))/as.numeric(as.character(sampleReport1[23,2])),
                                  as.numeric(as.character(sampleReport2[20,2]))/as.numeric(as.character(sampleReport2[23,2])),
                                  as.numeric(as.character(sampleReport3[20,2]))/as.numeric(as.character(sampleReport3[23,2])),
                                  as.numeric(as.character(sampleReport4[20,2]))/as.numeric(as.character(sampleReport4[23,2])),
                                  as.numeric(as.character(sampleReport5[20,2]))/as.numeric(as.character(sampleReport5[23,2])),
                                  as.numeric(as.character(sampleReport6[20,2]))/as.numeric(as.character(sampleReport6[23,2])))),
                           mean(c(as.numeric(as.character(sampleReport1[22,2]))/as.numeric(as.character(sampleReport1[23,2])),
                                  as.numeric(as.character(sampleReport2[22,2]))/as.numeric(as.character(sampleReport2[23,2])),
                                  as.numeric(as.character(sampleReport3[22,2]))/as.numeric(as.character(sampleReport3[23,2])),
                                  as.numeric(as.character(sampleReport4[22,2]))/as.numeric(as.character(sampleReport4[23,2])),
                                  as.numeric(as.character(sampleReport5[22,2]))/as.numeric(as.character(sampleReport5[23,2])),
                                  as.numeric(as.character(sampleReport6[22,2]))/as.numeric(as.character(sampleReport6[23,2]))))),
                    
                    sd=c(sd(c(as.numeric(as.character(sampleReport1[20,2]))/as.numeric(as.character(sampleReport1[23,2])),
                              as.numeric(as.character(sampleReport2[20,2]))/as.numeric(as.character(sampleReport2[23,2])),
                              as.numeric(as.character(sampleReport3[20,2]))/as.numeric(as.character(sampleReport3[23,2])),
                              as.numeric(as.character(sampleReport4[20,2]))/as.numeric(as.character(sampleReport4[23,2])),
                              as.numeric(as.character(sampleReport5[20,2]))/as.numeric(as.character(sampleReport5[23,2])),
                              as.numeric(as.character(sampleReport6[20,2]))/as.numeric(as.character(sampleReport6[23,2])))),
                         sd(c(as.numeric(as.character(sampleReport1[22,2]))/as.numeric(as.character(sampleReport1[23,2])),
                              as.numeric(as.character(sampleReport2[22,2]))/as.numeric(as.character(sampleReport2[23,2])),
                              as.numeric(as.character(sampleReport3[22,2]))/as.numeric(as.character(sampleReport3[23,2])),
                              as.numeric(as.character(sampleReport4[22,2]))/as.numeric(as.character(sampleReport4[23,2])),
                              as.numeric(as.character(sampleReport5[22,2]))/as.numeric(as.character(sampleReport5[23,2])),
                              as.numeric(as.character(sampleReport6[22,2]))/as.numeric(as.character(sampleReport6[23,2]))))))
  return(res)
}

human10 <- singleSpeciesRes("human","Human","10","SRR6373932","SRR6818517","SRR6373926","SRR6373923","SRR6825466","SRR6825471")
cattle10 <- singleSpeciesRes("cattle","Cattle","10","SRR7528450","SRR7528456","SRR7528458","SRR7528459","SRR7528464","SRR7528456")
pig10 <- singleSpeciesRes("pig","Pig","10","SRR7812176","SRR7812178","SRR7812179","SRR7812199","SRR7812200","SRR7812210")
allDiscordant <- rbind(human10,cattle10,pig10)

allDiscordant[which(allDiscordant$Species=="human"),1] <- "Human"
allDiscordant[which(allDiscordant$Species=="cattle"),1] <- "Cattle"
allDiscordant[which(allDiscordant$Species=="pig"),1] <- "Pig"

allDiscordant$Species <- factor(allDiscordant$Species,levels = c("Human","Cattle","Pig"))
allDiscordant$Type <- factor(allDiscordant$Type,levels = c("Discordant CpG","Concordant CpG"))

tiff("figure/cpg/speciesAndDepth10.tiff",width = 1500,height =1000,res = 300)
p <- ggplot(allDiscordant, aes(Species, weight = mean*100, fill = Type)) +
  geom_bar(color = "black", width = .8, position = 'stack') +
  labs( x="",y = 'Proportion (%)') +
  scale_fill_manual(values = c("pink3","paleturquoise3"))+
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour ="black",size = 11,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line= element_line(color="black"),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        legend.key.size = unit(15, "pt"),
        legend.title=element_text(colour ="black",size = 12,face = "bold") , legend.text=element_text(colour ="black",size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,colour ="black",face = "bold"))
print(p)
dev.off()
############## CpG: The proportion of discordant and concordant CpG sites in cgi and repeat ###########################################################
repAndCgi <- function(species,species2,depth,sampleName1,sampleName2,sampleName3,sampleName4,sampleName5,sampleName6){
  sampleReport1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName1,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName2,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName3,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName4,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName5,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName6,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleCgi1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName1,"resStaOfCgiAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleCgi2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName2,"resStaOfCgiAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleCgi3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName3,"resStaOfCgiAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleCgi4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName4,"resStaOfCgiAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleCgi5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName5,"resStaOfCgiAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleCgi6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName6,"resStaOfCgiAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleRep1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName1,"resStaOfRepAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleRep2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName2,"resStaOfRepAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleRep3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName3,"resStaOfRepAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleRep4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName4,"resStaOfRepAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleRep5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName5,"resStaOfRepAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleRep6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName6,"resStaOfRepAnno.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  
  repResSingle<- function(sampleRep,sampleReport){
    ## non-rep: discordant CpG=fourDetec_nonConsisCpgAnnoResRep+fourDiff_allAnnoResRep
    nonRepDiscordant <- sum(as.numeric(as.character(sampleRep[which(sampleRep$resClass=="."),c(3,4)])))
    ## non-rep: concordant CpG
    nonRepConcordant <- sum(as.numeric(as.character(sampleRep[which(sampleRep$resClass=="."),2])))    
    ## rep: discordant CpG = all_discordant_CpG  -  non_rep_discordant_CpG  
    repDiscordant <- (as.numeric(as.character(sampleReport[22,2]))-sum(as.numeric(as.character(sampleRep[which(sampleRep$resClass=="."),c(3,4)]))))
    ## rep: concordant CpG = all_concordant_CpG - non_rep_concordant_CpG 
    repConcordant <- as.numeric(as.character(sampleReport[20,2]))-as.numeric(as.character(sampleRep[which(sampleRep$resClass=="."),2]))
    temp <- c(nonRepDiscordant,nonRepConcordant,repDiscordant,repConcordant)
    return(temp)
  }
  
  repRes <- data.frame(Species=species2,
                       Region=c("Non-repeatitive Sequence","Non-repeatitive Sequence",
                                "Repeatitive Sequence","Repeatitive Sequence"),
                       Type=c("Discordant CpG","Concordant CpG","Discordant CpG","Concordant CpG"),
                       sample1=repResSingle(sampleRep1,sampleReport1),
                       sample2=repResSingle(sampleRep2,sampleReport2),
                       sample3=repResSingle(sampleRep3,sampleReport3),
                       sample4=repResSingle(sampleRep4,sampleReport4),
                       sample5=repResSingle(sampleRep5,sampleReport5),
                       sample6=repResSingle(sampleRep6,sampleReport6))
  repRes$Mean <- c(mean(as.numeric(repRes[1,4:9])),
                   mean(as.numeric(repRes[2,4:9])),
                   mean(as.numeric(repRes[3,4:9])),
                   mean(as.numeric(repRes[4,4:9])))
  repRes$Sd <- c(sd(as.numeric(repRes[1,4:9])),
                 sd(as.numeric(repRes[2,4:9])),
                 sd(as.numeric(repRes[3,4:9])),
                 sd(as.numeric(repRes[4,4:9])))
  repRes$Rate <- c(repRes$Mean[1]/(repRes$Mean[1]+repRes$Mean[2]),
                   repRes$Mean[2]/(repRes$Mean[1]+repRes$Mean[2]),
                   repRes$Mean[3]/(repRes$Mean[3]+repRes$Mean[4]),
                   repRes$Mean[4]/(repRes$Mean[3]+repRes$Mean[4]))
  
  
  cgiResSingle<- function(sampleCgi,sampleReport){
    ## CGI: Discordant CpG
    cgiDiscordant <- sum(as.numeric(as.character(sampleCgi[which(sampleCgi$resClass=="cpgIsland"),c(3,4)])))
    ## CGI: Concordant CpG
    cgiConcordant <- as.numeric(as.character(sampleCgi[which(sampleCgi$resClass=="cpgIsland"),2]))
    ## non-CGI: Discordant CpG = all_discordant_CpG - CGI_discordant_CpG
    nonCgiDiscordant <- as.numeric(as.character(sampleReport[22,2]))-sum(as.numeric(as.character(sampleCgi[which(sampleCgi$resClass=="cpgIsland"),c(3,4)])))
    ## non-CGI: Concordant CpG = all_concordant_CpG - non_CGI_concordant_CpG 
    nonCgiConcordant <- as.numeric(as.character(sampleReport[20,2]))-as.numeric(as.character(sampleCgi[which(sampleCgi$resClass=="cpgIsland"),2]))
    temp <- c(nonCgiDiscordant,nonCgiConcordant,cgiDiscordant,cgiConcordant)
    return(temp)
  }
  
  cgiRes <- data.frame(Species=species2,
                       Region=c("Non-CGI","Non-CGI","CGI","CGI"),
                       Type=c("Discordant CpG","Concordant CpG","Discordant CpG","Concordant CpG"),
                       sample1=cgiResSingle(sampleCgi1,sampleReport1),
                       sample2=cgiResSingle(sampleCgi2,sampleReport2),
                       sample3=cgiResSingle(sampleCgi3,sampleReport3),
                       sample4=cgiResSingle(sampleCgi4,sampleReport4),
                       sample5=cgiResSingle(sampleCgi5,sampleReport5),
                       sample6=cgiResSingle(sampleCgi6,sampleReport6))
  cgiRes$Mean <- c(mean(as.numeric(cgiRes[1,4:9])),
                   mean(as.numeric(cgiRes[2,4:9])),
                   mean(as.numeric(cgiRes[3,4:9])),
                   mean(as.numeric(cgiRes[4,4:9])))
  cgiRes$Sd <- c(sd(as.numeric(cgiRes[1,4:9])),
                 sd(as.numeric(cgiRes[2,4:9])),
                 sd(as.numeric(cgiRes[3,4:9])),
                 sd(as.numeric(cgiRes[4,4:9])))
  cgiRes$Rate <- c(cgiRes$Mean[1]/(cgiRes$Mean[1]+cgiRes$Mean[2]),
                   cgiRes$Mean[2]/(cgiRes$Mean[1]+cgiRes$Mean[2]),
                   cgiRes$Mean[3]/(cgiRes$Mean[3]+cgiRes$Mean[4]),
                   cgiRes$Mean[4]/(cgiRes$Mean[3]+cgiRes$Mean[4]))
  
  allRes <- rbind(repRes,cgiRes)
  allRes$SpeRegion <- c(paste(species2,": ",allRes$Region,sep = ""))
  return(allRes)
}
human10 <- repAndCgi("human","Human","10","SRR6373932","SRR6818517","SRR6373926","SRR6373923","SRR6825466","SRR6825471")
cattle10 <- repAndCgi("cattle","Cattle","10","SRR7528450","SRR7528456","SRR7528458","SRR7528459","SRR7528464","SRR7528456")
pig10 <- repAndCgi("pig","Pig","10","SRR7812176","SRR7812178","SRR7812179","SRR7812199","SRR7812200","SRR7812210")
allSpecies10 <- rbind(human10,pig10,cattle10)
allSpecies10$SpeRegion <- factor(allSpecies10$SpeRegion,
                                 levels = c("Pig: Non-repeatitive Sequence","Cattle: Non-repeatitive Sequence",
                                            "Human: Non-repeatitive Sequence",
                                            "Pig: Repeatitive Sequence","Cattle: Repeatitive Sequence",
                                            "Human: Repeatitive Sequence",
                                            "Pig: Non-CGI","Cattle: Non-CGI","Human: Non-CGI",
                                            "Pig: CGI","Cattle: CGI","Human: CGI"))
allSpecies10$Type <- factor(allSpecies10$Type,levels =  c("Discordant CpG","Concordant CpG"))
tiff("figure\\cpg\\Dpeth10RepAndCGI.tiff",width = 1950,height =750,res = 300)
ggplot(allSpecies10, aes(SpeRegion, weight = Rate*100, fill = Type)) +
  geom_bar(color = "black", width = .8, position = 'stack') +
  labs( x="",y = 'Proportion (%)') +
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c("pink3","paleturquoise3"))+
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour =c(rep("lightskyblue4",6),rep("palevioletred4",6)),size = 11,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        legend.key.size = unit(15, "pt"),
        legend.title=element_text(colour ="black",size = 12,face = "bold") , legend.text=element_text(colour ="black",size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,colour ="black",face = "bold"))
dev.off()
############## CpG: The proportion of discordant and concordant CpG sites in different methylation levels#######################
methDis <- function(species,species2,depth,sample1,sample2,sample3,sample4,sample5,sample6){
  sample1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sample1,"cpgMehtLevelDistri.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sample2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sample2,"cpgMehtLevelDistri.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sample3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sample3,"cpgMehtLevelDistri.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sample4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sample4,"cpgMehtLevelDistri.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sample5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sample5,"cpgMehtLevelDistri.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sample6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sample6,"cpgMehtLevelDistri.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  
  sam<- function(sampleName){
    # 样本不同甲基化水平的总和,总和是四个软件cpg位点的并集
    total <- as.numeric(colSums(sampleName[,2:4]))
    # 样本不同甲基化水平的一致性总数
    con <- as.numeric(sampleName[which(sampleName$class=="consis"),2:4])
    # 样本不同甲基化水平的非一致性总数
    dis <- as.numeric(colSums(sampleName[which(sampleName$class!="consis"),2:4]))
    # 三种甲基化水平的total,concordant,discordant
    te <- c(total,con,dis)
    return(te)
  }
  
  allSam<- data.frame(Species=species2,
                      Meth=rep(c("Hight Methylation","Intermediate Methylation","Low Methylation"),3),
                      Class=c(rep("total CpG",3),rep("Concordant CpG",3),rep("Discordant CpG",3)),
                      Sample1=sam(sample1),
                      Sample2=sam(sample2),
                      Sample3=sam(sample3),
                      Sample4=sam(sample4),
                      Sample5=sam(sample5),
                      Sample6=sam(sample6))
  allSam$Mean <- rowMeans(allSam[,4:9])
  allSam$Rate <- c(allSam$Mean[1]/allSam$Mean[1],
                   allSam$Mean[2]/allSam$Mean[2],
                   allSam$Mean[3]/allSam$Mean[3],
                   allSam$Mean[4]/allSam$Mean[1],
                   allSam$Mean[5]/allSam$Mean[2],
                   allSam$Mean[6]/allSam$Mean[3],
                   allSam$Mean[7]/allSam$Mean[1],
                   allSam$Mean[8]/allSam$Mean[2],
                   allSam$Mean[9]/allSam$Mean[3])
  allSam$SpeRegion <- c(paste0(species2,": ",allSam$Meth))
  return(allSam)
}
human10 <- methDis("human","Human","10","SRR6373932","SRR6818517","SRR6373926","SRR6373923","SRR6825466","SRR6825471")
cattle10 <- methDis("cattle","Cattle","10","SRR7528450","SRR7528456","SRR7528458","SRR7528459","SRR7528464","SRR7528456")
pig10 <- methDis("pig","Pig","10","SRR7812176","SRR7812178","SRR7812179","SRR7812199","SRR7812200","SRR7812210")
allSpec10 <- rbind(human10,pig10,cattle10)
allSpecies10 <- allSpec10[which(allSpec10$Class!="total CpG"),]
allSpecies10$SpeRegion <- factor(allSpecies10$SpeRegion,
                                 levels = c("Pig: Hight Methylation","Cattle: Hight Methylation",
                                            "Human: Hight Methylation",
                                            "Pig: Intermediate Methylation","Cattle: Intermediate Methylation",
                                            "Human: Intermediate Methylation",
                                            "Pig: Low Methylation","Cattle: Low Methylation",
                                            "Human: Low Methylation"))
allSpecies10$Class <- factor(allSpecies10$Class,levels = c("Discordant CpG","Concordant CpG"))
tiff("figure\\cpg\\Dpeth10MethAll.tiff",width = 1880,height =650,res = 300)
ggplot(allSpecies10, aes(SpeRegion, weight = Rate*100, fill = Class)) +
  geom_bar(color = "black", width = .8, position = 'stack') +
  labs( x="",y = 'Proportion (%)') +
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c("pink3","paleturquoise3"))+
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour =c(rep("orchid4",3),rep("darkseagreen4",3),rep("lightskyblue4",3)),size = 11,face = "bold"),
        #    axis.text.y = element_text(colour =c(rep("orchid4",3),rep("darkseagreen4",3),rep("lightskyblue4",3),rep("palevioletred4",3)),size = 11,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        legend.key.size = unit(15, "pt"),
        legend.title=element_text(colour ="black",size = 12,face = "bold") , legend.text=element_text(colour ="black",size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,colour ="black",face = "bold"))

dev.off()
############## CpG: venn for four alignment algorithms ############################################################################
countToalNumCpG <- function(species,sampleName1,sampleName2,sampleName3,sampleName4,sampleName5,sampleName6,depth){
  sampleReport1 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName1,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport2 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName2,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport3 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName3,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport4 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName4,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport5 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName5,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport6 <- read.csv(paste(species,"\\cpgAnaRes\\dpeth",depth,"_",sampleName6,"_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  
  me<- function(sampleReport1,sampleReport2,sampleReport3,sampleReport4,sampleReport5,sampleReport6){
    res <- c()
    for(i in 5:19){
      temp<- mean(c(as.numeric(as.character(sampleReport1[i,2])),
                    as.numeric(as.character(sampleReport2[i,2])),
                    as.numeric(as.character(sampleReport3[i,2])),
                    as.numeric(as.character(sampleReport4[i,2])),
                    as.numeric(as.character(sampleReport5[i,2])),
                    as.numeric(as.character(sampleReport6[i,2]))))
      
      res <- c(res,temp)
    }
    return(res)
  }
  
  classStr <- function(sampleReport1){
    res <- c()
    for(i in 5:19){
      temp <- str_sub(as.character(sampleReport1[i,1]),end = -4)
      res <- c(res,temp)
    }
    return(res)
  }
  
  res<- data.frame(class=classStr(sampleReport1),
                   value=c(me(sampleReport1,sampleReport2,sampleReport3,sampleReport4,sampleReport5,sampleReport6)))
  value <- function(class){
    temp<- res[which(res$class==class),2]
    temp <- as.numeric(temp)
    return(temp)
  }
  sumTotal <- value("a")+value("b")+value("c")+value("d")+value("abc")+value("abd")+value("acd")+value("bcd")-
    value("ab")-value("ac")-value("ad")-value("bc")-value("bd")-value("cd")-value("abcd")
  #  return(sumTotal)
  tiff(paste("figure/cpg/",species,"_",depth,"_DiscordantCpG.tiff",sep = ""),width = 1850,height =1800,res = 300)
  p <- draw.quad.venn(
    area1=value("c"), area2=value("d"), area3=value("b"), area4=value("a"),
    n12=value("cd"), n13=value("bc"), n23=value("bd"), n14=value("ac"), n24=value("ad") , n34=value("ab"), 
    n123=value("bcd"), n124=value("acd"), n234=value("abd"), n134=value("abc"), 
    n1234=value("abcd") ,category = c('Bwa-meth','Walt','BSMAP','Bismark-bwt2-e2e') ,
    col=c("pink3","darkseagreen3","khaki3","lightskyblue3"), fill=c("pink3","darkseagreen3","khaki3","lightskyblue3") ,
    alpha = 0.7 ,
    cex=1.4,
    cat.cex=1.4,
    fontface = "bold",
    cat.fontface= "bold",
    print.mode="percent",
    sigdigs=4)
  print(p)
  dev.off()
}

countToalNumCpG("human","SRR6373932","SRR6818517","SRR6373926","SRR6373923","SRR6825466","SRR6825471","10")
countToalNumCpG("cattle","SRR7528450","SRR7528456","SRR7528458","SRR7528459","SRR7528464","SRR7528456","10")
countToalNumCpG("pig","SRR7812176","SRR7812178","SRR7812179","SRR7812199","SRR7812200","SRR7812210","10")


############## DMC: The total number of DMC  #######################################################
DMLres<- function(species,species2,depth){
  sampleReport <- read.csv(paste(species,"\\dmcRes\\depth",depth,"\\depth",depth,"report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleRes<- sampleReport[5:8,]
  sampleRes$Mapper <- c("Bismark-bwt2-e2e","BSMAP","Bwa-meth","Walt")
  sampleRes$Species <- species2
  colnames(sampleRes) <- c("Class","Number","Mapper","Species")
  return(sampleRes)  
}
humanRes<- DMLres("human","Human","10")
cattleRes<- DMLres("cattle","Cattle","10")
pigRes<- DMLres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)
allRes$Number <- as.numeric(allRes$Number)
allRes$Mapper <- factor(allRes$Mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/dmc/numOfDmcMapper.tiff",width = 2200,height =1540,res = 300)
ggplot(allRes, aes(Mapper, weight = Number)) +
  geom_bar(color = "black",fill='khaki3', width = .8, position = 'stack') +
  labs( x="",y = 'The number of DMC')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 47, hjust = 1,colour ="black",face = "bold",size = 13),
        axis.text.y = element_text(colour ="black",face = "bold",size = 13),
        axis.title.y = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        axis.line = element_line(color="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 4)
dev.off()
############## DMC: The number of shared DMC ######################

DMLres<- function(species,species2,depth){
  sampleReport <- read.csv(paste(species,"\\dmcRes\\depth",depth,"\\depth",depth,"report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleRes<- sampleReport[15:19,]
  sampleRes$Type <- c("Bwa,BSM,Bis","BSM,Bis,Wal","Bwa,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis,Wal")
  sampleRes$Species <- species2
  sampleRes$SpeMap <- c(paste0(species,": Bwa,BSM,Bis"),
                        paste0(species,": BSM,Bis,Wal"),
                        paste0(species,": Bwa,Bis,Wal"),
                        paste0(species,": Bwa,BSM,Wal"),
                        paste0(species,": Bwa,BSM,Bis,Wal"))
  colnames(sampleRes)[1:2] <- c("Class","Number")
  sampleRes$Number <- as.numeric(sampleRes$Number)
  return(sampleRes)  
}

humanRes<- DMLres("human","Human","10")
cattleRes<- DMLres("cattle","Cattle","10")
pigRes<- DMLres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)
allRes$Type <- factor(allRes$Type,levels = c("Bwa,BSM,Bis,Wal","Bwa,Bis,Wal","BSM,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))

tiff("figure/dmc/concordantDMC.tiff",width = 2400,height =1740,res = 300)
ggplot(allRes,aes(x=Type,y=Number,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'The number of shared DMC')+
  scale_y_continuous(expand = c(0,0),limits=c(0,8.2e+05),breaks = seq(0,8.2e+05,2e+05)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 4)
dev.off()

############## DMC：The number of reliable DMCs ##############################################

DMLres <- function(species,species2,depth){
  sampleReport1 <- read.csv(paste(species,"\\dmcRes\\depth",depth,"\\depth",depth,"report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport <-sampleReport1[-c(1:4),] 
  sampleReport$V2 <- as.numeric(sampleReport$V2)
  singleSpe<- function(data,species){
    res<- data.frame(value=c(sum(data[c(11,13,14),2])-2*data[15,2], ## Bwa-meth 
                             sum(data[c(11,12,14),2])-2*data[15,2], ## BSMAP
                             sum(data[c(11,12,13),2])-2*data[15,2], ## Bismark-bwt2-e2e
                             sum(data[c(12,13,14),2])-2*data[15,2]), ## Walt
                     Rate=c((sum(data[c(11,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bwa-meth 
                            (sum(data[c(11,12,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## BSMAP
                            (sum(data[c(11,12,13),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bismark-bwt2-e2e
                            (sum(data[c(12,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2])), ## Walt
                     mapper=c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"),
                     Class=c("Covered by three alignment algorithms in Bwa-meth",
                             "Covered by three alignment algorithms in BSMAP",
                             "Covered by three alignment algorithms in Bismark-bwt2-e2e",
                             "Covered by three alignment algorithms in Walt"),
                     Species=species,
                     SpeCla=c(paste0(species,": Covered by three alignment algorithms in Bwa-meth"),
                              paste0(species,": Covered by three alignment algorithms in BSMAP"),
                              paste0(species,": Covered by three alignment algorithms in Bismark-bwt2-e2e"),
                              paste0(species,": Covered by three alignment algorithms in Walt")))
    return(res)
  }
  temp<- singleSpe(sampleReport,species2)
  return(temp)
}

humanRes<- DMLres("human","Human","10")
cattleRes<- DMLres("cattle","Cattle","10")
pigRes<- DMLres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)
allRes$mapper <- factor(allRes$mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/20220113/reliable_DMC_rate.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=Rate*100,group=Species))+
  labs( x="",y = 'The proportion of accurate DMCs (%)',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  scale_y_continuous(expand = c(0,0),limits=c(78,98),breaks=seq(78,98,5)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()


tiff("figure/20220113/reliable_DMC_value.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=value,group=Species))+
  labs( x="",y = 'The number of accurate DMCs',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  scale_y_continuous(expand = c(0,0),limits=c(4e+05,1.2e+06),breaks=seq(4e+05,1.2e+06,2e+05)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()
############## DMC: The proportion of discordant and concordant DMC in cgi and repeat ###############
sinSpe <- function(species,species2,depth){
  sampleCgi <- read.csv(paste(species,"\\dmcRes\\depth",depth,"\\depth",depth,"dmcAnnoResStaCGI.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  sampleRep <- read.csv(paste(species,"\\dmcRes\\depth",depth,"\\depth",depth,"dmcAnnoResStaRep.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  
  rate <- function(data,class){
    a <- data[which(data$resClass==class),3]/(data[which(data$resClass==class),3]+data[which(data$resClass==class),2])
    return(a)
  }
  repRes <- data.frame(species=species2,
                       Type=c("Discordant DMC","Concordant DMC","Discordant DMC","Concordant DMC"),
                       region=c("Repeatitive Sequence","Repeatitive Sequence","Non-repeatitive Sequence","Non-repeatitive Sequence"),
                       Value=c(rate(sampleRep,"repeat"),1-rate(sampleRep,"repeat"),
                               rate(sampleRep,"nonRepeat"),1-rate(sampleRep,"nonRepeat")))
  repRes$SpeRegion <- c(paste(species2,": ",repRes$region,sep = ""))
  cgiRes <- data.frame(species=species2,
                       Type=c("Discordant DMC","Concordant DMC",
                              "Discordant DMC","Concordant DMC"),
                       region=c("CGI","CGI","Non-CGI","Non-CGI"),
                       Value=c(rate(sampleCgi,"cgi"),1-rate(sampleCgi,"cgi"),
                               rate(sampleCgi,"nonCgi"),1-rate(sampleCgi,"nonCgi")))
  cgiRes$SpeRegion <- c(paste(species2,": ",cgiRes$region,sep = ""))
  cgiAndRepRes <- rbind(repRes,cgiRes)
  return(cgiAndRepRes)
}

human <- sinSpe("human","Human","10")
cattle <- sinSpe("cattle","Cattle","10")
pig <- sinSpe("pig","Pig","10")
cgiAndRepRes <- rbind(human,cattle,pig)
cgiAndRepRes$SpeRegion <- factor(cgiAndRepRes$SpeRegion,
                                 levels = c("Pig: Non-repeatitive Sequence","Cattle: Non-repeatitive Sequence",
                                            "Human: Non-repeatitive Sequence",
                                            "Pig: Repeatitive Sequence","Cattle: Repeatitive Sequence",
                                            "Human: Repeatitive Sequence",
                                            "Pig: Non-CGI","Cattle: Non-CGI","Human: Non-CGI",
                                            "Pig: CGI","Cattle: CGI","Human: CGI"))
cgiAndRepRes$Type <- factor(cgiAndRepRes$Type,levels = c("Discordant DMC","Concordant DMC"))

tiff("figure\\dmc\\repAndCGIDMC.tiff",width = 1950,height =750,res = 300)
ggplot(cgiAndRepRes, aes(SpeRegion, weight = Value*100, fill = Type)) +
  geom_bar(color = "black", width = .8, position = 'stack') +
  labs( x="",y = 'Proportion (%)') +
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c("khaki3","paleturquoise3"))+
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour =c(rep("lightskyblue4",6),rep("palevioletred4",6)),size = 11,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        legend.key.size = unit(15, "pt"),
        legend.title=element_text(colour ="black",size = 12,face = "bold") , legend.text=element_text(colour ="black",size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,colour ="black",face = "bold"))
dev.off()

############## DMC: venn for four alignment algorithms ###############################################################
veenDisDmc<- function(species,depth){
  sampleReport <- read.csv(paste(species,"\\dmcRes\\depth",depth,"\\depth",depth,"report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  me<- function(sampleReport){
    res <- c()
    for(i in 5:19){
      temp<- as.numeric(as.character(sampleReport[i,2]))
      res <- c(res,temp)
    }
    return(res)
  }
  classStr <- function(sampleReport1){
    res <- c()
    for(i in 5:19){
      temp <- str_sub(as.character(sampleReport1[i,1]),end = -4)
      res <- c(res,temp)
    }
    return(res)
  }
  res<- data.frame(class=classStr(sampleReport),
                   value=c(me(sampleReport)))
  value <- function(class){
    temp<- res[which(res$class==class),2]
    return(temp)
  }
  tiff(paste("figure/dmc/",species,"_DiscordantDmc","_depth",depth,".tiff",sep = ""),width = 1850,height =1800,res = 300)
  p <- draw.quad.venn(
    area1=value("c"), area2=value("d"), area3=value("b"), area4=value("a"),
    n12=value("cd"), n13=value("bc"), n23=value("bd"), n14=value("ac"), n24=value("ad") , n34=value("ab"), 
    n123=value("bcd"), n124=value("acd"), n234=value("abd"), n134=value("abc"), 
    n1234=value("abcd") ,category = c('Bwa-meth','Walt','BSMAP','Bismark-bwt2-e2e') ,
    col=c("pink3","darkseagreen3","khaki3","lightskyblue3"), fill=c("pink3","darkseagreen3","khaki3","lightskyblue3") ,
    alpha = 0.7 ,
    cex=1.4,
    cat.cex=1.4,
    fontface = "bold",
    cat.fontface= "bold",
    print.mode="percent",
    sigdigs=4)
  print(p)
  dev.off()
}

veenDisDmc("human","10")
veenDisDmc("cattle","10")
veenDisDmc("pig","10")

############## DMR: The total length of DMR ################################
DMRres<- function(species,species2,depth){
  sampleReport <- read.csv(paste(species,"\\dmrRes\\depth",depth,"\\",species,"_vennRes.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleRes<- sampleReport[5:8,]
  sampleRes$Mapper <- c("Bismark-bwt2-e2e","BSMAP","Bwa-meth","Walt")
  sampleRes$Species <- species2
  colnames(sampleRes) <- c("Class","Number","Mapper","Species")
  return(sampleRes)  
}

humanRes<- DMRres("human","Human","10")
cattleRes<- DMRres("cattle","Cattle","10")
pigRes<- DMRres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)
allRes$Number <- as.numeric(allRes$Number)
allRes$Mapper <- factor(allRes$Mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/dmr/numOfDmrMapper.tiff",width = 2200,height =1540,res = 300)
ggplot(allRes, aes(Mapper, weight = Number)) +
  geom_bar(color = "black",fill='lightskyblue3', width = .8, position = 'stack') +
  labs( x="",y = 'The length of DMR')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 47, hjust = 1,colour ="black",face = "bold",size = 13),
        axis.text.y = element_text(colour ="black",face = "bold",size = 13),
        axis.title.y = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        axis.line = element_line(color="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 4)
dev.off()
############## DMR: The length of shared DMR ######################################
DMRres<- function(species,species2,depth){
  sampleReport <- read.csv(paste(species,"\\dmrRes\\depth",depth,"\\",species,"_vennRes.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleRes<- sampleReport[15:19,]
  sampleRes$Type <- c("Bwa,BSM,Bis","BSM,Bis,Wal","Bwa,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis,Wal")
  sampleRes$Species <- species2
  sampleRes$SpeMap <- c(paste0(species,": Bwa,BSM,Bis"),
                        paste0(species,": BSM,Bis,Wal"),
                        paste0(species,": Bwa,Bis,Wal"),
                        paste0(species,": Bwa,BSM,Wal"),
                        paste0(species,": Bwa,BSM,Bis,Wal"))
  colnames(sampleRes)[1:2] <- c("Class","Number")
  sampleRes$Number <- as.numeric(sampleRes$Number)
  return(sampleRes)  
}

humanRes<- DMRres("human","Human","10")
cattleRes<- DMRres("cattle","Cattle","10")
pigRes<- DMRres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)

allRes$Type <- factor(allRes$Type,levels = c("Bwa,BSM,Bis,Wal","Bwa,Bis,Wal","BSM,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))


tiff("figure/dmr/concordantDMR.tiff",width = 2400,height =1740,res = 300)
ggplot(allRes,aes(x=Type,y=Number,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'The length of shared DMR')+
  scale_y_continuous(expand = c(0,0),limits=c(0,4.3e+07),breaks = seq(0,4.3e+07,1e+07)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") , legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 4)
dev.off()

############## DMR: The proportion of discordant and concordant DMR in cgi and repeat #################################################
singleSpeciesResDmr <- function(species,species2,depth){
  sampleReport <- read.csv(paste(species,"\\dmrRes\\depth",depth,"\\",species,"_rep_cgi_report.txt",sep=""),header = T,sep = "\t", quote = "\"", dec = ".")
  regionAll<- function(sampleReport){
    sampleReport$type <- as.character(sampleReport$type)
    temp <- rbind(c(as.character(sampleReport[1,1]),"repeat","all",sum(sampleReport$length[which(sampleReport$region=="repeat")])),
                  c(as.character(sampleReport[1,1]),"non_repeat","all",sum(sampleReport$length[which(sampleReport$region=="non_repeat")])),
                  c(as.character(sampleReport[1,1]),"CGI","all",sum(sampleReport$length[which(sampleReport$region=="CGI")])),
                  c(as.character(sampleReport[1,1]),"non_CGI","all",sum(sampleReport$length[which(sampleReport$region=="non_CGI")])))
    temp <- as.data.frame(temp)
    names(temp) <- names(sampleReport)
    sampleReport <- rbind(sampleReport,temp)
    sampleReport$length<- as.numeric(sampleReport$length)
    return(sampleReport)
  }
  sampleReport<- regionAll(sampleReport)
  returnVal<- function(sampleReport,regionName,typeName){
    value <- sampleReport$length[which(sampleReport$region==regionName & sampleReport$type==typeName)]/
      sampleReport$length[which(sampleReport$region==regionName & sampleReport$type=="all")]
    return(value)
  }
  
  res <- data.frame(Species=species2,
                    Region=c("Repeatitive Sequence","Repeatitive Sequence","Non-repeatitive Sequence","Non-repeatitive Sequence",
                             "CGI","CGI","Non-CGI","Non-CGI"),
                    Type=c("Discordant DMR","Concordant DMR","Discordant DMR","Concordant DMR","Discordant DMR","Concordant DMR",
                           "Discordant DMR","Concordant DMR"),
                    mean=c(returnVal(sampleReport,"repeat","discordant"),
                           returnVal(sampleReport,"repeat","concordant"),
                           returnVal(sampleReport,"non_repeat","discordant"),
                           returnVal(sampleReport,"non_repeat","concordant"),
                           returnVal(sampleReport,"CGI","discordant"),
                           returnVal(sampleReport,"CGI","concordant"),
                           returnVal(sampleReport,"non_CGI","discordant"),
                           returnVal(sampleReport,"non_CGI","concordant")))
  res$SpeRegion <- c(paste(res$Species,": ",res$Region,sep = ""))
  return(res)
}

human <- singleSpeciesResDmr("human","Human",10)
cattle <- singleSpeciesResDmr("cattle","Cattle",10)
pig <- singleSpeciesResDmr("pig","Pig",10)
allSpecies <- rbind(human,cattle,pig)

allSpecies$SpeRegion <- factor(allSpecies$SpeRegion,
                               levels = c("Pig: Non-repeatitive Sequence","Cattle: Non-repeatitive Sequence",
                                          "Human: Non-repeatitive Sequence",
                                          "Pig: Repeatitive Sequence","Cattle: Repeatitive Sequence",
                                          "Human: Repeatitive Sequence",
                                          "Pig: Non-CGI","Cattle: Non-CGI","Human: Non-CGI",
                                          "Pig: CGI","Cattle: CGI","Human: CGI"))
allSpecies$Type <- factor(allSpecies$Type,levels = c("Discordant DMR","Concordant DMR"))
tiff("figure\\dmr\\repAndCGI.tiff",width = 1950,height =750,res = 300)
ggplot(allSpecies, aes(SpeRegion, weight = mean*100, fill = Type)) +
  geom_bar(color = "black", width = .8, position = 'stack') +
  labs( x="",y = 'Proportion (%)') +
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c("lightskyblue3","paleturquoise3"))+
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_text(colour =c(rep("lightskyblue4",6),rep("palevioletred4",6)),size = 11,face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size=12,face = "bold"),
        legend.key.size = unit(15, "pt"),
        legend.title=element_text(colour ="black",size = 12,face = "bold") , legend.text=element_text(colour ="black",size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,colour ="black",face = "bold"))
dev.off()

############## DMR：The number of reliable DMRs #############################################
DMRres<- function(species,species2,depth){
  sampleReport1 <- read.csv(paste(species,"\\dmrRes\\depth",depth,"\\",species,"_vennRes.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport <-sampleReport1[-c(1:4),] 
  sampleReport$V2 <- as.numeric(sampleReport$V2)
  singleSpe<- function(data,species){
    res<- data.frame(value=c(sum(data[c(11,13,14),2])-2*data[15,2], ## Bwa-meth 
                             sum(data[c(11,12,14),2])-2*data[15,2], ## BSMAP
                             sum(data[c(11,12,13),2])-2*data[15,2], ## Bismark-bwt2-e2e
                             sum(data[c(12,13,14),2])-2*data[15,2]), ## Walt
                     Rate=c((sum(data[c(11,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bwa-meth 
                            (sum(data[c(11,12,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## BSMAP
                            (sum(data[c(11,12,13),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bismark-bwt2-e2e
                            (sum(data[c(12,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2])), ## Walt
                     mapper=c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"),
                     Class=c("Covered by three alignment algorithms in Bwa-meth",
                             "Covered by three alignment algorithms in BSMAP",
                             "Covered by three alignment algorithms in Bismark-bwt2-e2e",
                             "Covered by three alignment algorithms in Walt"),
                     Species=species,
                     SpeCla=c(paste0(species,": Covered by three alignment algorithms in Bwa-meth"),
                              paste0(species,": Covered by three alignment algorithms in BSMAP"),
                              paste0(species,": Covered by three alignment algorithms in Bismark-bwt2-e2e"),
                              paste0(species,": Covered by three alignment algorithms in Walt")))
    return(res)
  }
  temp<- singleSpe(sampleReport,species2)
  return(temp)
}

humanRes<- DMRres("human","Human","10")
cattleRes<- DMRres("cattle","Cattle","10")
pigRes<- DMRres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)
allRes$mapper <- factor(allRes$mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/20220113/reliable_DMR_rate.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=Rate*100,group=Species))+
  labs( x="",y = 'The proportion of accurate DMRs (%)',shape='Algorithms',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  scale_y_continuous(expand = c(0,0),limits=c(79,97),breaks=seq(81,97,4)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()


tiff("figure/20220113/reliable_DMR_value.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=value,group=Species))+
  labs( x="",y = 'The length of accurate DMRs',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  scale_y_continuous(expand = c(0,0),limits=c(3.5e+07,6.0e+07),breaks=seq(3.5e+07,6.0e+07,0.5e+07)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()
############## DMR: Venn for four alignment algorithms ##########################################
dmrVenn <- function(species,depth){
  sampleReport1 <- read.csv(paste(species,"\\dmrRes\\depth",depth,"\\",species,"_vennRes.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  me<- function(sampleReport1){
    res <- c()
    for(i in 5:19){
      temp<- as.numeric(as.character(sampleReport1[i,2]))
      res <- c(res,temp)
    }
    return(res)
  }
  res<- data.frame(class=as.character(sampleReport1$V1[5:19]),
                   value=c(me(sampleReport1)))
  value <- function(class){
    temp<- res[which(res$class==class),2]
    return(temp)
  }
  tiff(paste("figure/dmr/",species,"_ConcordantDmr_depth",depth,".tiff",sep = ""),width = 1850,height =1800,res = 300)
  p<- draw.quad.venn(
    area1=value("c"), area2=value("d"), area3=value("b"), area4=value("a"),
    n12=value("cd"), n13=value("bc"), n23=value("bd"), n14=value("ac"), n24=value("ad") , n34=value("ab"), 
    n123=value("bcd"), n124=value("acd"), n234=value("abd"), n134=value("abc"), 
    n1234=value("abcd") ,category = c('Bwa-meth','Walt','BSMAP','Bismark-bwt2-e2e') ,
    col=c("pink3","darkseagreen3","khaki3","lightskyblue3"), fill=c("pink3","darkseagreen3","khaki3","lightskyblue3") ,
    alpha = 0.7 ,
    cex=1.4,
    cat.cex=1.4,
    fontface = "bold",
    cat.fontface= "bold",
    print.mode="percent",
    sigdigs=4)
  print(p)
  dev.off()
}

dmrVenn("human","10")
dmrVenn("cattle","10")
dmrVenn("pig","10")

############## DMR-related genes: The total number of genes ############################

GENEres<- function(species,species2,depth){
  sampleReport <- read.csv(paste(species,"\\gene\\depth",depth,"\\","depth",depth,"_gene_DMR_Venn_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleRes<- sampleReport[5:8,]
  sampleRes$Mapper <- c("Bismark-bwt2-e2e","BSMAP","Bwa-meth","Walt")
  sampleRes$Species <- species2
  colnames(sampleRes) <- c("Class","Number","Mapper","Species")
  return(sampleRes)  
}

humanRes<- GENEres("human","Human","10")
cattleRes<- GENEres("cattle","Cattle","10")
pigRes<- GENEres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)
allRes$Number <- as.numeric(allRes$Number)
allRes$Mapper <- factor(allRes$Mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/gene/numOfgeneMapper.tiff",width = 2200,height =1540,res = 300)
ggplot(allRes, aes(Mapper, weight = Number)) +
  geom_bar(color = "black",fill='darkseagreen3', width = .8, position = 'stack') +
  labs( x="",y = 'The number of Gene')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 47, hjust = 1,colour ="black",face = "bold",size = 13),
        axis.text.y = element_text(colour ="black",face = "bold",size = 13),
        axis.title.y = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        axis.line = element_line(color="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 4)
dev.off()
############## DMR-related genes: The number of shared gene ##########################################################
GENEres<- function(species,species2,depth){
  sampleReport <- read.csv(paste(species,"\\gene\\depth",depth,"\\","depth",depth,"_gene_DMR_Venn_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleRes<- sampleReport[15:19,]
  sampleRes$Type <- c("Bwa,BSM,Bis","BSM,Bis,Wal","Bwa,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis,Wal")
  sampleRes$Species <- species2
  sampleRes$SpeMap <- c(paste0(species,": Bwa,BSM,Bis"),
                        paste0(species,": BSM,Bis,Wal"),
                        paste0(species,": Bwa,Bis,Wal"),
                        paste0(species,": Bwa,BSM,Wal"),
                        paste0(species,": Bwa,BSM,Bis,Wal"))
  colnames(sampleRes)[1:2] <- c("Class","Number")
  sampleRes$Number <- as.numeric(sampleRes$Number)
  return(sampleRes) 
}

humanRes<- GENEres("human","Human","10")
cattleRes<- GENEres("cattle","Cattle","10")
pigRes<- GENEres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)

allRes$Type <- factor(allRes$Type,levels = c("Bwa,BSM,Bis,Wal","Bwa,Bis,Wal","BSM,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))

tiff("figure/gene/concordantGENE.tiff",width = 2400,height =1740,res = 300)
ggplot(allRes,aes(x=Type,y=Number,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'The number of shared gene')+
  scale_y_continuous(expand = c(0,0),limits=c(10000,28000),breaks = seq(10000,28000,6000)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") , legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()

############## DMR-related genes：The numnber of reliable genes #################################################
GENEres<- function(species,species2,depth){
  sampleReport1 <- read.csv(paste(species,"\\gene\\depth",depth,"\\","depth",depth,"_gene_DMR_Venn_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport <-sampleReport1[-c(1:4),] 
  sampleReport$V2 <- as.numeric(sampleReport$V2)
  singleSpe<- function(data,species){
    res<- data.frame(value=c(sum(data[c(11,13,14),2])-2*data[15,2], ## Bwa-meth 
                             sum(data[c(11,12,14),2])-2*data[15,2], ## BSMAP
                             sum(data[c(11,12,13),2])-2*data[15,2], ## Bismark-bwt2-e2e
                             sum(data[c(12,13,14),2])-2*data[15,2]), ## Walt
                     Rate=c((sum(data[c(11,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bwa-meth 
                            (sum(data[c(11,12,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## BSMAP
                            (sum(data[c(11,12,13),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bismark-bwt2-e2e
                            (sum(data[c(12,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2])), ## Walt
                     mapper=c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"),
                     Class=c("Covered by three alignment algorithms in Bwa-meth",
                             "Covered by three alignment algorithms in BSMAP",
                             "Covered by three alignment algorithms in Bismark-bwt2-e2e",
                             "Covered by three alignment algorithms in Walt"),
                     Species=species,
                     SpeCla=c(paste0(species,": Covered by three alignment algorithms in Bwa-meth"),
                              paste0(species,": Covered by three alignment algorithms in BSMAP"),
                              paste0(species,": Covered by three alignment algorithms in Bismark-bwt2-e2e"),
                              paste0(species,": Covered by three alignment algorithms in Walt")))
    return(res)
  }
  temp<- singleSpe(sampleReport,species2)
  return(temp)
}

humanRes<- GENEres("human","Human","10")
cattleRes<- GENEres("cattle","Cattle","10")
pigRes<- GENEres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)
allRes$mapper <- factor(allRes$mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/20220113/reliable_gene_rate.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=Rate*100,group=Species))+
  labs( x="",y = 'The proportion of accurate genes (%)',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  scale_y_continuous(expand = c(0,0),limits=c(94,100),breaks=seq(94,100,2)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()


tiff("figure/20220113/reliable_gene_value.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=value,group=Species))+
  labs( x="",y = 'The length of accurate genes',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  #  scale_y_continuous(expand = c(0,0),limits=c(3.5e+07,6.0e+07),breaks=seq(3.5e+07,6.0e+07,0.5e+07)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()
############## DMR-related genes: venn for four alignment algorithms ###############################################################
geneVenn <- function(species,depth){
  sampleReport1 <- read.csv(paste(species,"\\gene\\depth",depth,"\\depth",depth,"_gene_DMR_Venn_report.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  me<- function(sampleReport1){
    res <- c()
    for(i in 5:19){
      temp<- as.numeric(as.character(sampleReport1[i,2]))
      res <- c(res,temp)
    }
    return(res)
  }
  res<- data.frame(class=as.character(sampleReport1$V1[5:19]),
                   value=c(me(sampleReport1)))
  value <- function(class){
    temp<- res[which(res$class==class),2]
    return(temp)
  }
  tiff(paste("figure/gene/",species,"_ConcordantGene_depth",depth,".tiff",sep = ""),width = 1850,height =1800,res = 300)
  p<- draw.quad.venn(
    area1=value("cNum"), area2=value("dNum"), area3=value("bNum"), area4=value("aNum"),
    n12=value("cdNum"), n13=value("bcNum"), n23=value("bdNum"), n14=value("acNum"), n24=value("adNum") , n34=value("abNum"), 
    n123=value("bcdNum"), n124=value("acdNum"), n234=value("abdNum"), n134=value("abcNum"), 
    n1234=value("abcdNum") ,category = c('Bwa-meth','Walt','BSMAP','Bismark-bwt2-e2e') ,
    col=c("pink3","darkseagreen3","khaki3","lightskyblue3"), fill=c("pink3","darkseagreen3","khaki3","lightskyblue3") ,
    alpha = 0.7 ,
    cex=1.4,
    cat.cex=1.4,
    fontface = "bold",
    cat.fontface= "bold",
    print.mode="percent",
    sigdigs=4)
  print(p)
  dev.off()
}

geneVenn("human","10")
geneVenn("cattle","10")
geneVenn("pig","10")


############## KEGG: The total number of KEGG ##############################################

KEGGres<- function(species,species2,depth){
  sampleReport <- read.csv(paste("enrichKEGG\\depth",depth,"_p001\\",species,"_kegg_Venn.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleRes<- sampleReport[5:8,]
  sampleRes$Mapper <- c("Bismark-bwt2-e2e","BSMAP","Bwa-meth","Walt")
  sampleRes$Species <- species2
  colnames(sampleRes) <- c("Class","Number","Mapper","Species")
  return(sampleRes)  
}

humanRes<- KEGGres("human","Human","10")
cattleRes<- KEGGres("cattle","Cattle","10")
pigRes<- KEGGres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)
allRes$Number <- as.numeric(allRes$Number)
allRes$Mapper <- factor(allRes$Mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/kegg/numOfkeggMapper.tiff",width = 2200,height =1540,res = 300)
ggplot(allRes, aes(Mapper, weight = Number)) +
  geom_bar(color = "black",fill='paleturquoise3', width = .8, position = 'stack') +
  labs( x="",y = 'The number of pathway')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 47, hjust = 1,colour ="black",face = "bold",size = 13),
        axis.text.y = element_text(colour ="black",face = "bold",size = 13),
        axis.title.y = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        axis.line = element_line(color="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()
############## KEGG: The number of shared pathway ######################################################
KEGGres<- function(species,species2,depth){
  sampleReport <- read.csv(paste("enrichKEGG\\depth",depth,"_p001\\",species,"_kegg_Venn.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleRes<- sampleReport[15:19,]
  sampleRes$Type <- c("Bwa,BSM,Bis","BSM,Bis,Wal","Bwa,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis,Wal")
  sampleRes$Species <- species2
  sampleRes$SpeMap <- c(paste0(species,": Bwa,BSM,Bis"),
                        paste0(species,": BSM,Bis,Wal"),
                        paste0(species,": Bwa,Bis,Wal"),
                        paste0(species,": Bwa,BSM,Wal"),
                        paste0(species,": Bwa,BSM,Bis,Wal"))
  colnames(sampleRes)[1:2] <- c("Class","Number")
  sampleRes$Number <- as.numeric(sampleRes$Number)
  return(sampleRes) 
}

humanRes<- KEGGres("human","Human","10")
cattleRes<- KEGGres("cattle","Cattle","10")
pigRes<- KEGGres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)

allRes$Type <- factor(allRes$Type,levels = c("Bwa,BSM,Bis,Wal","Bwa,Bis,Wal","BSM,Bis,Wal","Bwa,BSM,Wal","Bwa,BSM,Bis"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))

tiff("figure/kegg/concordantkegg.tiff",width = 2400,height =1740,res = 300)
ggplot(allRes,aes(x=Type,y=Number,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'The number of shared pathway')+
#  scale_y_continuous(expand = c(0,0),limits=c(100,230),breaks = seq(100,230,40)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") , legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 4)
dev.off()

############## KEGG：The number of reliable pathway ############################################
KEGGres<- function(species,species2,depth){
  sampleReport1 <- read.csv(paste("enrichKEGG\\depth",depth,"_p001\\",species,"_kegg_Venn.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  sampleReport <-sampleReport1[-c(1:4),] 
  sampleReport$V2 <- as.numeric(sampleReport$V2)
  singleSpe<- function(data,species){
    res<- data.frame(value=c(sum(data[c(11,13,14),2])-2*data[15,2], ## Bwa-meth 
                             sum(data[c(11,12,14),2])-2*data[15,2], ## BSMAP
                             sum(data[c(11,12,13),2])-2*data[15,2], ## Bismark-bwt2-e2e
                             sum(data[c(12,13,14),2])-2*data[15,2]), ## Walt
                     Rate=c((sum(data[c(11,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bwa-meth 
                            (sum(data[c(11,12,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## BSMAP
                            (sum(data[c(11,12,13),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2]), ## Bismark-bwt2-e2e
                            (sum(data[c(12,13,14),2])-2*data[15,2])/(sum(data[c(11,12,13,14),2])-3*data[15,2])), ## Walt
                     mapper=c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"),
                     Class=c("Covered by three alignment algorithms in Bwa-meth",
                             "Covered by three alignment algorithms in BSMAP",
                             "Covered by three alignment algorithms in Bismark-bwt2-e2e",
                             "Covered by three alignment algorithms in Walt"),
                     Species=species,
                     SpeCla=c(paste0(species,": Covered by three alignment algorithms in Bwa-meth"),
                              paste0(species,": Covered by three alignment algorithms in BSMAP"),
                              paste0(species,": Covered by three alignment algorithms in Bismark-bwt2-e2e"),
                              paste0(species,": Covered by three alignment algorithms in Walt")))
    return(res)
  }
  temp<- singleSpe(sampleReport,species2)
  return(temp)
}

humanRes<- KEGGres("human","Human","10")
cattleRes<- KEGGres("cattle","Cattle","10")
pigRes<- KEGGres("pig","Pig","10")
allRes <- rbind(humanRes,cattleRes,pigRes)

allRes$mapper <- factor(allRes$mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt"))
allRes$Species <- factor(allRes$Species,levels = c("Human","Cattle","Pig"))
tiff("figure/20220113/reliable_pathway_rate.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=Rate*100,group=Species))+
  labs( x="",y = 'The proportion of accurate pathways (%)',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  scale_y_continuous(expand = c(0,0),limits=c(90,100),breaks=seq(91,100,3)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()


tiff("figure/20220113/reliable_pathway_value.tiff",width = 2100,height =1700,res = 300)
ggplot(allRes,aes(x=mapper,y=value,group=Species))+
  labs( x="",y = 'The length of accurate pathways',shape='Algorithms')+
  geom_point(aes(shape=mapper,colour=Species),size =5,stroke = 2) +
  scale_shape_manual(values = c(1,2,3,4))+
  scale_y_continuous(expand = c(0,0),limits=c(75,205),breaks=seq(80,200,30)) +
  theme(axis.title.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 15,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 47, hjust = 1, size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black",face = "bold") ,
        legend.text=element_text(size=13,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        strip.text = element_text(colour ="black",face = "bold",size = 13))+
  facet_wrap( .~ Species,ncol = 3)
dev.off()


############## KEGG: venn for for alignment algorithms #############################################################
keggVenn <- function(species,depth){
  sampleReport1 <- read.csv(paste("enrichKEGG\\depth",depth,"_p001\\",species,"_kegg_Venn.txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
  me<- function(sampleReport1){
    res <- c()
    for(i in 5:19){
      temp<- as.numeric(as.character(sampleReport1[i,2]))
      res <- c(res,temp)
    }
    return(res)
  }
  res<- data.frame(class=as.character(sampleReport1$V1[5:19]),
                   value=c(me(sampleReport1)))
  value <- function(class){
    temp<- res[which(res$class==class),2]
    return(temp)
  }
  tiff(paste("figure/kegg/",species,"_Concordantkegg_depth",depth,".tiff",sep = ""),width = 1850,height =1800,res = 300)
  p<- draw.quad.venn(
    area1=value("cNum"), area2=value("dNum"), area3=value("bNum"), area4=value("aNum"),
    n12=value("cdNum"), n13=value("bcNum"), n23=value("bdNum"), n14=value("acNum"), n24=value("adNum") , n34=value("abNum"), 
    n123=value("bcdNum"), n124=value("acdNum"), n234=value("abdNum"), n134=value("abcNum"), 
    n1234=value("abcdNum") ,category = c('Bwa-meth','Walt','BSMAP','Bismark-bwt2-e2e') ,
    col=c("pink3","darkseagreen3","khaki3","lightskyblue3"), fill=c("pink3","darkseagreen3","khaki3","lightskyblue3") ,
    alpha = 0.7 ,
    cex=1.4,
    cat.cex=1.4,
    fontface = "bold",
    cat.fontface= "bold",
    print.mode="percent",
    sigdigs=4)
  print(p)
  dev.off()
}

keggVenn("human","10")
keggVenn("cattle","10")
keggVenn("pig","10")
############## KEGG: top 30 pathway ################################
enrichKEGGPlot <- function(species,speciesName,depth,class,wid){
  mapperRes30<- function(class,depth,species,mapperName1,mapperName2){
    mapperRes <- read.csv(paste("enrich",class,"\\depth",depth,"_p001\\",species,"_",mapperName1,"_enrich",class,"_depth",depth,".txt",sep=""),header = F,sep = "\t", quote = "\"", dec = ".")
    mapperRes$Mapper <- c(mapperName2)
    names(mapperRes) <- c("ID","Description","GeneRatio","BaRatio","pvalue","P.adjust","qvalue","geneID","Count","Mapper")
    mapperRes <- mapperRes[order(mapperRes$P.adjust),]
    if (length(mapperRes$ID)>30){
      mapperRes <- mapperRes[1:30,]
    }
    mapperRes$Description <- factor(mapperRes$Description,levels =rev(mapperRes$Description) )
    return(mapperRes)
  }
  
  bismark<- mapperRes30(class,depth,species,"bismarkbwt2","Bismark-bwt2-e2e")
  bsmap <- mapperRes30(class,depth,species,"bsmap","BSMAP")
  bwameth <- mapperRes30(class,depth,species,"bwameth","Bwa-meth")
  walt <- mapperRes30(class,depth,species,"walt","Walt")
  go <- rbind(bwameth,bsmap,bismark, walt)
  go <- within(go, Mapper <- factor(go$Mapper,levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt")))
  
  tiff(paste("figure\\kegg\\",species,"_depth",depth,"_enrich",class,".tiff",sep = ""),
       width = wid,height =1800,res = 300)
  p <- ggplot(go,aes(Mapper,Description)) + 
    geom_point(aes(size=Count,color=-1*log(P.adjust)))+
    scale_color_gradient(low = "blue", high = "red")+ 
    scale_size_continuous(range=c(1,3))+
    labs(color=expression(-log[10](P.adjust)),
         size="Gene Number",
         x="-log10(P.adjust)",
         y="",
         title=speciesName,face="bold")+
    theme_bw()+
    guides(color=guide_colourbar(order = 2),size=guide_legend(order=1))+
    theme(axis.text.y = element_text(size=6,face = "bold",colour ="black"),
          axis.text.x =element_text(size = 7,face = "bold",angle = 15,hjust = 0.5,vjust=0.8,colour ="black"),
          axis.title.x = element_text(colour ="black"),plot.title=element_text(size=9,colour ="black"),
          legend.title=element_text(colour ="black") , legend.text=element_text(colour ="black"))
  print(p)
  dev.off()
}

enrichKEGGPlot("human","Human","10","KEGG",1715)
enrichKEGGPlot("cattle","Cattle","10","KEGG",1739)
enrichKEGGPlot("pig","Pig","10","KEGG",1740)


