setwd("..\\WGBS\\result\\bench5")
library(ggplot2)
library(pheatmap)
library(ggpubr)

################## simulated data: five error rate #####################
SingleSpeciesRes<- function(species,species2){
  bench <- read.csv(paste("metaData\\",species,"BenchSimuDataAccuUniConbine.csv",sep = ""),header = TRUE, sep = "\t", quote = "\"", dec = ".")
  bench$cpuTime<- bench$cpusysTime+bench$cpuuserTime
  bench$errorRate<- paste(bench$errorRate,"%",sep = "", collapse = NULL)
  bench$RSS <- bench$RSS/1000000
  bench$cpuTime <- bench$cpuTime/60
  bench$matchedReadsRate <- bench$matchedReads/2000000
  
  bench_mean <- aggregate(bench[,c(6,9,15,20,21)],by=list(bench$errorRate,bench$mapper),mean)
  names(bench_mean)[1:7] <- c("errorRate","mapper","macroAvgPrecisionMean","microAvgPrecisionMean","RSSMean","cpuTimeMean","matchedReadsRateMean")
  bench_sd <- aggregate(bench[,c(6,9,15,20,21)],by=list(bench$errorRate,bench$mapper),sd)
  names(bench_sd)[3:7] <- c("macroAvgPrecisionSd","microAvgPrecisionSd","RSSSd","cpuTimeSd","matchedReadsRateSd")
  bench_mean_sd <- cbind(bench_mean,bench_sd[,3:7])
  bench_mean_sd$mapper <- as.character(bench_mean_sd$mapper)
  bench_mean_sd[bench_mean_sd$mapper=="batmeth2",2] <- "Batmeth2"
  bench_mean_sd[bench_mean_sd$mapper=="bismarkbwt2",2] <- "Bismark-bwt2-e2e"
  bench_mean_sd[bench_mean_sd$mapper=="bismarkhis2",2] <- "Bismark-his2"
  bench_mean_sd[bench_mean_sd$mapper=="bsmap",2] <- "BSMAP"
  bench_mean_sd[bench_mean_sd$mapper=="bsseeker2bt",2] <- "BSSeeker2-bwt"
  bench_mean_sd[bench_mean_sd$mapper=="bsseeker2bt2end",2] <- "BSSeeker2-bwt2-e2e"
  bench_mean_sd[bench_mean_sd$mapper=="bsseeker2bt2loc",2] <- "BSSeeker2-bwt2-local"
  bench_mean_sd[bench_mean_sd$mapper=="bsseeker2soap",2] <- "BSSeeker2-soap2"
  bench_mean_sd[bench_mean_sd$mapper=="bwameth",2] <- "Bwa-meth"
  bench_mean_sd[bench_mean_sd$mapper=="walt",2] <- "Walt"
  bench_mean_sd[bench_mean_sd$errorRate=="0%",1] <- "0"
  bench_mean_sd[bench_mean_sd$errorRate=="0.5%",1] <- "0.50%"
  bench_mean_sd[bench_mean_sd$errorRate=="1%",1] <- "1.00%"
  
  names(bench_mean_sd)[1] <- "Error rate"
  bench_mean_sd$mapper <- factor(bench_mean_sd$mapper,
                                 levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt","Batmeth2",
                                            "BSSeeker2-bwt2-local","BSSeeker2-soap2","Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
  bench_mean_sd$SpeMap <- paste(species2,": ",bench_mean_sd$mapper,sep = "")
  bench_mean_sd$matchedReadsRateMean <- bench_mean_sd$matchedReadsRateMean*100
  bench_mean_sd$macroAvgPrecisionMean <- bench_mean_sd$macroAvgPrecisionMean*100
  return(bench_mean_sd)
}
human<- SingleSpeciesRes("human","Human")
mouse <-SingleSpeciesRes("mouse","Mouse")
cattle <- SingleSpeciesRes("cattle","Cattle")
pig <- SingleSpeciesRes("pig","Pig")
all <- rbind(pig,cattle,mouse,human)
all$`Error rate` <- factor(all$`Error rate`,levels=sort(unique(all$`Error rate`),decreasing = T))

################## uniquely mapped reads ###############################
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\uniqueMap.tiff",width = 1800,height =1640,res = 300)
ggplot(all, aes(SpeMap, weight = matchedReadsRateMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Uniquely mapped read (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,500),breaks=seq(0,500,100)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",4),rep("indianred4",4),
                                             rep("hotpink3",4),rep("salmon3",4),
                                             rep("azure4",4),rep("tan4",4),
                                             rep("orchid4",4),rep("darkseagreen4",4),
                                             rep("lightskyblue4",4),rep("palevioletred4",4)),size = 8,
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        legend.title=element_text(colour ="black",face = "bold") , legend.text=element_text(colour ="black",face = "bold"),
        axis.text.x = element_text(size = 11,colour ="black",face = "bold"),
        axis.title.x = element_text(face = "bold"))
dev.off()

################## mapped precisoin #################################################
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\mapPre.tiff",width = 1800,height =1640,res = 300)
ggplot(all, aes(SpeMap, weight = macroAvgPrecisionMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Mapped precision (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,500),breaks=seq(0,500,100)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",4),rep("indianred4",4),
                                             rep("hotpink3",4),rep("salmon3",4),
                                             rep("azure4",4),rep("tan4",4),
                                             rep("orchid4",4),rep("darkseagreen4",4),
                                             rep("lightskyblue4",4),rep("palevioletred4",4)),size = 8,
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        legend.title=element_text(colour ="black",face = "bold") , legend.text=element_text(colour ="black",face = "bold"),
        axis.text.x = element_text(size = 11,colour ="black",face = "bold"),
        axis.title.x = element_text(face = "bold"))
dev.off()

################## rss ################################
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\RSS.tiff",width = 1800,height =1640,res = 300)
ggplot(all, aes(SpeMap, weight = RSSMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Memory consumption (GB)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,400),breaks=seq(0,400,80)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",4),rep("indianred4",4),
                                             rep("hotpink3",4),rep("salmon3",4),
                                             rep("azure4",4),rep("tan4",4),
                                             rep("orchid4",4),rep("darkseagreen4",4),
                                             rep("lightskyblue4",4),rep("palevioletred4",4)),size = 8,
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        legend.title=element_text(colour ="black",face = "bold") , legend.text=element_text(colour ="black",face = "bold"),
        axis.text.x = element_text(size = 11,colour ="black",face = "bold"),
        axis.title.x = element_text(face = "bold"))
dev.off()
################## runtime ##########
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\runtime.tiff",width = 1800,height =1640,res = 300)
ggplot(all, aes(SpeMap, weight = cpuTimeMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Runtime (min)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,1500),breaks=seq(0,1500,300)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",4),rep("indianred4",4),
                                             rep("hotpink3",4),rep("salmon3",4),
                                             rep("azure4",4),rep("tan4",4),
                                             rep("orchid4",4),rep("darkseagreen4",4),
                                             rep("lightskyblue4",4),rep("palevioletred4",4)),size = 8,
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        legend.title=element_text(colour ="black",face = "bold") , legend.text=element_text(colour ="black",face = "bold"),
        axis.text.x = element_text(size = 11,colour ="black",face = "bold"),
        axis.title.x = element_text(face = "bold"))
dev.off()
################## simulated data: average of five error rate ######################################################
dataRehape <- function(species,species2){
  fileName <- read.csv(paste("metaData\\",species,"BenchSimuDataAccuUniConbine.csv",sep = ""),header = T,sep = "\t", quote = "\"", dec = ".")
  fileName$cpuTime <- fileName$cpusysTime+fileName$cpuuserTime
  fileName$RSS <- fileName$RSS/1000000
  fileName$cpuTime <- fileName$cpuTime/60
  fileName$matchedReadsRate <- fileName$matchedReads/2000000
  #  fileName<- fileName[,c(1,2,3,6,15,20,21,19)]
  fileName_mean <- aggregate(fileName[,c(6,9,15,20,21)],by=list(fileName$mapper,fileName$species),mean)
  names(fileName_mean)[1:7] <- c("mapper","species","macroAvgPrecisionMean","microAvgPrecisionMean","RSSMean","cpuTimeMean","matchedReadsRateMean")
  fileName_sd <- aggregate(fileName[,c(6,9,15,20,21)],by=list(fileName$mapper,fileName$species),sd)
  names(fileName_sd)[3:7] <- c("macroAvgPrecisionSd","microAvgPrecisionSd","RSSSd","cpuTimeSd","matchedReadsRateSd")
  fileName_mean_sd <- cbind(fileName_mean,fileName_sd[,3:7])
  fileName_mean_sd$mapper <- as.character(fileName_mean_sd$mapper)
  oldName <- c("batmeth2","bismarkbwt2","bismarkhis2","bsmap","bsseeker2bt","bsseeker2bt2end","bsseeker2bt2loc",
               "bsseeker2soap","bwameth","walt")
  newName <- c("Batmeth2","Bismark-bwt2-e2e","Bismark-his2","BSMAP","BSSeeker2-bwt","BSSeeker2-bwt2-e2e","BSSeeker2-bwt2-local",
               "BSSeeker2-soap2","Bwa-meth","Walt")
  for(i in 1:10){
    fileName_mean_sd$mapper[which(fileName_mean_sd$mapper==oldName[i])] <- c(newName[i])
  }
  fileName_mean_sd$species <- species2
  names(fileName_mean_sd)[2] <- "Species"
  return(fileName_mean_sd)
}

human <- dataRehape("human","Human")
mouse <- dataRehape("mouse","Mouse")
pig <- dataRehape("pig","Pig")
cattle <- dataRehape("cattle","Cattle")
dataFrame <- rbind(human,mouse,pig,cattle)
dataFrame$mapper <- factor(dataFrame$mapper,levels =c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt", "Batmeth2","BSSeeker2-bwt2-local",
                                                      "BSSeeker2-soap2","Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
dataFrame$matchedReadsRateMean <- dataFrame$matchedReadsRateMean*100
dataFrame$Species <- factor(dataFrame$Species,levels = c("Human","Mouse","Cattle","Pig"))
################## runtime #######
tiff("figure/runTimeAllErr.tiff",width = 2100,height =1700,res = 300)
p<-ggplot(dataFrame,aes(x=mapper,y=cpuTimeMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Runtime (min)')+
  scale_y_continuous(limits=c(0,300),breaks=seq(0,300,50)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black"),
        axis.text.y = element_text( size = 15,colour ="black"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 15,colour ="black"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black") , legend.text=element_text(size=13,colour ="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
print(p)
dev.off()

################## rss ##############
tiff("figure/rssAllErr.tiff",width = 2100,height =1700,res = 300)
p<-ggplot(dataFrame,aes(x=mapper,y=RSSMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Memory consumption (GB)')+
  scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black"),
        axis.text.y = element_text( size = 15,colour ="black"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 15,colour ="black"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black") , legend.text=element_text(size=13,colour ="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
print(p)
dev.off()

################## uniquely mapped reads ######
tiff("figure/uniMapAllErr.tiff",width = 2100,height =1700,res = 300)
p<-ggplot(dataFrame,aes(x=mapper,y=matchedReadsRateMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Uniquely mapped reads (%)')+
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
  theme(axis.title.y = element_text( size = 15,colour ="black"),
        axis.text.y = element_text( size = 15,colour ="black"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 15,colour ="black"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black") , legend.text=element_text(size=13,colour ="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
print(p)
dev.off()
################## precisoin #####################

tiff("figure/precisionAllErr.tiff",width = 2100,height =1700,res = 300)
p<-ggplot(dataFrame,aes(x=mapper,y=macroAvgPrecisionMean*100,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Mapped precision (%)')+
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
  theme(axis.title.y = element_text( size = 15,colour ="black"),
        axis.text.y = element_text( size = 15,colour ="black"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 15,colour ="black"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=14,colour ="black") , legend.text=element_text(size=13,colour ="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
print(p)
dev.off()
################## cor #############################
SingleSpeciesRes<- function(species,species2){
  bench <- read.csv(paste("metaData\\",species,"BenchSimuDataAccuUniConbine.csv",sep = ""),header = TRUE, sep = "\t", quote = "\"", dec = ".")
  bench$cpuTime<- bench$cpusysTime+bench$cpuuserTime
  bench$errorRate<- paste(bench$errorRate,"%",sep = "", collapse = NULL)
  bench$RSS <- bench$RSS/1000000
  bench$cpuTime <- bench$cpuTime/60
  bench$matchedReadsRate <- bench$matchedReads/2000000
  
  bench_mean <- aggregate(bench[,c(6,9,15,20,21)],by=list(bench$errorRate,bench$mapper),mean)
  names(bench_mean)[1:7] <- c("errorRate","mapper","macroAvgPrecisionMean","microAvgPrecisionMean","RSSMean","cpuTimeMean","matchedReadsRateMean")
  bench_sd <- aggregate(bench[,c(6,9,15,20,21)],by=list(bench$errorRate,bench$mapper),sd)
  names(bench_sd)[3:7] <- c("macroAvgPrecisionSd","microAvgPrecisionSd","RSSSd","cpuTimeSd","matchedReadsRateSd")
  bench_mean_sd <- cbind(bench_mean,bench_sd[,3:7])
  bench_mean_sd$mapper <- as.character(bench_mean_sd$mapper)
  bench_mean_sd[bench_mean_sd$mapper=="batmeth2",2] <- "Batmeth2"
  bench_mean_sd[bench_mean_sd$mapper=="bismarkbwt2",2] <- "Bismark-bwt2-e2e"
  bench_mean_sd[bench_mean_sd$mapper=="bismarkhis2",2] <- "Bismark-his2"
  bench_mean_sd[bench_mean_sd$mapper=="bsmap",2] <- "BSMAP"
  bench_mean_sd[bench_mean_sd$mapper=="bsseeker2bt",2] <- "BSSeeker2-bwt"
  bench_mean_sd[bench_mean_sd$mapper=="bsseeker2bt2end",2] <- "BSSeeker2-bwt2-e2e"
  bench_mean_sd[bench_mean_sd$mapper=="bsseeker2bt2loc",2] <- "BSSeeker2-bwt2-local"
  bench_mean_sd[bench_mean_sd$mapper=="bsseeker2soap",2] <- "BSSeeker2-soap"
  bench_mean_sd[bench_mean_sd$mapper=="bwameth",2] <- "Bwa-meth"
  bench_mean_sd[bench_mean_sd$mapper=="walt",2] <- "Walt"
  bench_mean_sd[bench_mean_sd$errorRate=="0%",1] <- "0"
  bench_mean_sd[bench_mean_sd$errorRate=="0.5%",1] <- "0.50%"
  bench_mean_sd[bench_mean_sd$errorRate=="1%",1] <- "1.00%"
  
  names(bench_mean_sd)[1] <- "Error rate"
  bench_mean_sd$mapper <- factor(bench_mean_sd$mapper,
                                 levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt","Batmeth2",
                                            "BSSeeker2-bwt2-local","BSSeeker2-soap","Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
  bench_mean_sd$SpeMap <- paste(species2,": ",bench_mean_sd$mapper,sep = "")
  bench_mean_sd$matchedReadsRateMean <- bench_mean_sd$matchedReadsRateMean*100
  bench_mean_sd$macroAvgPrecisionMean <- bench_mean_sd$macroAvgPrecisionMean*100
  return(bench_mean_sd)
}

human<- SingleSpeciesRes("human","Human")
mouse <-SingleSpeciesRes("mouse","Mouse")
cattle <- SingleSpeciesRes("cattle","Cattle")
pig <- SingleSpeciesRes("pig","Pig")
human$species <- "Human"
mouse$species <- "Mouse"
cattle$species <- "Cattle"
pig$species <- "Pig"
all <- rbind(pig,cattle,mouse,human)
for(i in 1:length(all$`Error rate`)){
  if(all$`Error rate`[i]=="0"){all$`Error rate`[i]=0}
  if(all$`Error rate`[i]=="0.25%"){all$`Error rate`[i]=0.25}
  if(all$`Error rate`[i]=="0.50%"){all$`Error rate`[i]=0.50}
  if(all$`Error rate`[i]=="0.75%"){all$`Error rate`[i]=0.75}
  if(all$`Error rate`[i]=="1.00%"){all$`Error rate`[i]=1.00}
}

errCorClass<- function(all,line,class){
  res <- data.frame(species=NULL,Mpper=NULL,Cor=NULL,Class=NULL)
  for (spe in unique(all$species)) {
    for (map in unique(all$mapper)){
      temp <- all[which(all$species==spe & all$mapper==map),]
      co <- cor.test(as.numeric(temp$`Error rate`),temp[,line])
      tempData<- data.frame(Species=spe,
                            Mapper=map,
                            Cor=as.numeric(co$estimate),
                            Pcor=as.numeric(co$p.value))
      res <- rbind(tempData,res)
    }
  }
  names(res)[3] <- paste(class,"Cor",sep = "")
  names(res)[4] <- paste(class,"CorP",sep = "")
  return(res)
}
cpu <- errCorClass(all,6,"CpuTime")
rss <- errCorClass(all,5,"RSS")
pre <- errCorClass(all,3,"Precision")
uni <- errCorClass(all,7,"Unique")

res <- cbind(cpu,rss$RSSCor,rss$RSSCorP,pre$PrecisionCor,pre$PrecisionCorP,uni$UniqueCor,uni$UniqueCorP)
names(res)[5:10] <- c("RSScor","RSSCorP","PrecisionCor","PrecisionCorP","UniqueCor","UniqueCorP")
res <- res[order(res$Mapper),]
#write.csv(res,file = "significance\\corTest.csv",row.names = F)

# pheatmap
res$SpeMap <- paste(res$Species,res$Mapper,sep = ": ")
res <- res[c(2,12,22,32,7,17,27,37,9,19,29,39,1,11,21,31,3,13,23,33,
             4,14,24,34,10,20,30,40,6,16,26,36,8,18,28,38,5,15,25,35),]
corErr <- res[,c(3,5,9,7)]
row.names(corErr) <- res$SpeMap
names(corErr) <- c("Runtime","Memory consumption","Uniquely mapped reads","Mapped precision")
`Mapping algorithms`<- factor(c(rep("Bwa-meth",4),rep("BSMAP",4),rep("Bismark-bwt2-e2e",4),rep("Walt",4),rep("BSSeeker2-soap2",4),
                                rep("BSSeeker2-bwt2-local",4),rep("Batmeth2",4),rep("BSSeeker2-bwt",4),rep("Bismark-his2",4),rep("BSSeeker2-bwt2-e2e",4)),
                              levels = c("Bwa-meth","BSMAP","Bismark-bwt2-e2e","Walt","BSSeeker2-soap2","BSSeeker2-bwt2-local",
                                         "Batmeth2","BSSeeker2-bwt","Bismark-his2","BSSeeker2-bwt2-e2e"))
`Mapping algorithms` <- data.frame(`Mapping algorithms`)
names(`Mapping algorithms`) <- "Alignment algorithms"
rownames(`Mapping algorithms`)=rownames(corErr)
ann_colors = list(
  `Alignment algorithms`=c(`Bwa-meth`="palevioletred4",BSMAP="lightskyblue4",	`Bismark-bwt2-e2e`="darkseagreen4",
                           Walt="orchid4",`BSSeeker2-soap2`="tan4",`BSSeeker2-bwt2-local`="azure4",
                           Batmeth2="salmon3",`BSSeeker2-bwt`="hotpink3",`Bismark-his2`="indianred4",
                           `BSSeeker2-bwt2-e2e`="seagreen4")
)
bk <- c(seq(-1,-0.001,by=0.0001),seq(0,1,by=0.0001))

tiff("figure\\heatCor.tiff",width = 2000,height =1940,res = 300)
pheatmap(corErr,annotation_row =`Mapping algorithms`,annotation_colors = ann_colors,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks = bk,
         legend_breaks = seq(-1,1,1),legend_labels = c("-1", "0", "1"),
         cluster_row = FALSE, cluster_col = FALSE
         #         ,filename = "figure/heatCot.tiff"
)
dev.off()

## cor:Scatter plot
all$`Error rate` <- as.numeric(all$`Error rate`)
#runtime
for(map in unique(all$mapper)){
  temp<- all[which(all$mapper==map),]
  temp$species <- factor(temp$species,levels = c("Human","Mouse","Cattle","Pig"))
  tiff(paste0("figure/","runtime","ScatterPlot/",map,"_","runtime",".tiff"),width = 2100,height =1700,res = 300)
  p<- ggplot(temp,aes(x=`Error rate`,y=cpuTimeMean))+
    geom_point(aes(color=species))+
    labs( x="",y = '')+
    stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
    stat_cor(aes(color=species), method = "pearson",size =6)+
    theme(axis.title = element_text( size = 15,colour ="black"),
          axis.text = element_text( size = 15,colour ="black"),axis.line = element_line(color="black"),
          panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
  print(p)
  dev.off()
}

#memory
for(map in unique(all$mapper)){
  temp<- all[which(all$mapper==map),]
  temp$species <- factor(temp$species,levels = c("Human","Mouse","Cattle","Pig"))
  tiff(paste0("figure/","memory","ScatterPlot/",map,"_","memory",".tiff"),width = 2100,height =1700,res = 300)
  p<- ggplot(temp,aes(x=`Error rate`,y=RSSMean))+
    geom_point(aes(color=species))+
    labs( x="",y = '')+
    stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
    stat_cor(aes(color=species), method = "pearson",size =6)+
    theme(axis.title = element_text( size = 15,colour ="black"),
          axis.text = element_text( size = 15,colour ="black"),axis.line = element_line(color="black"),
          panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
  print(p)
  dev.off()
}

#unique
for(map in unique(all$mapper)){
  temp<- all[which(all$mapper==map),]
  temp$species <- factor(temp$species,levels = c("Human","Mouse","Cattle","Pig"))
  tiff(paste0("figure/","unique","ScatterPlot/",map,"_","unique",".tiff"),width = 2100,height =1700,res = 300)
  p<- ggplot(temp,aes(x=`Error rate`,y=matchedReadsRateMean))+
    geom_point(aes(color=species))+
    labs( x="",y = '')+
    stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
    stat_cor(aes(color=species), method = "pearson",size =6)+
    theme(axis.title = element_text( size = 15,colour ="black"),
          axis.text = element_text( size = 15,colour ="black"),axis.line = element_line(color="black"),
          panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
  print(p)
  dev.off()
}

#precision
for(map in unique(all$mapper)){
  temp<- all[which(all$mapper==map),]
  temp$species <- factor(temp$species,levels = c("Human","Mouse","Cattle","Pig"))
  tiff(paste0("figure/","precision","ScatterPlot/",map,"_","precision",".tiff"),width = 2100,height =1700,res = 300)
  p<- ggplot(temp,aes(x=`Error rate`,y=macroAvgPrecisionMean))+
    geom_point(aes(color=species))+
    labs( x="",y = '')+
    stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
    stat_cor(aes(color=species), method = "pearson",size =6)+
    theme(axis.title = element_text( size = 15,colour ="black"),
          axis.text = element_text( size = 15,colour ="black"),axis.line = element_line(color="black"),
          panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
  print(p)
  dev.off()
}


