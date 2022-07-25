setwd("..\\WGBS\\result\\bench")
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
  colnames(bench)
  bench_mean <- aggregate(bench[,c(6,7,8,15,20,21)],by=list(bench$errorRate,bench$mapper),mean)
  names(bench_mean)[1:8] <- c("errorRate","mapper","macroAvgPrecisionMean","macroAvgRecallMean","macroAvgF1ScoreMean","RSSMean","cpuTimeMean","matchedReadsRateMean")
  bench_sd <- aggregate(bench[,c(6,7,8,15,20,21)],by=list(bench$errorRate,bench$mapper),sd)
  names(bench_sd)[3:8] <- c("macroAvgPrecisionSd","macroAvgRecallSd","macroAvgF1ScoreSd","RSSSd","cpuTimeSd","matchedReadsRateSd")
  bench_mean_sd <- cbind(bench_mean,bench_sd[,3:8])
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
  bench_mean_sd[bench_mean_sd$mapper=="hisat_3n",2] <- "Hisat_3n"
  bench_mean_sd[bench_mean_sd$mapper=="hisat_3n_repeat",2] <- "Hisat_3n_repeat"
  bench_mean_sd[bench_mean_sd$mapper=="bsbolt",2] <- "BSBolt"
  bench_mean_sd[bench_mean_sd$mapper=="abismal",2] <- "Abismal"
  bench_mean_sd[bench_mean_sd$errorRate=="0%",1] <- "0"
  bench_mean_sd[bench_mean_sd$errorRate=="0.5%",1] <- "0.50%"
  bench_mean_sd[bench_mean_sd$errorRate=="1%",1] <- "1.00%"
  names(bench_mean_sd)[1] <- "Error rate"
  bench_mean_sd$mapper <- factor(bench_mean_sd$mapper,
                                 levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt","Batmeth2",
                                            "BSSeeker2-bwt2-local","BSSeeker2-soap2","Walt","Bismark-bwt2-e2e","Abismal",
                                            "BSMAP","BSBolt","Hisat_3n","Hisat_3n_repeat","Bwa-meth"))
  bench_mean_sd$SpeMap <- paste(species2,": ",bench_mean_sd$mapper,sep = "")
  bench_mean_sd$matchedReadsRateMean <- bench_mean_sd$matchedReadsRateMean*100
  bench_mean_sd$macroAvgPrecisionMean <- bench_mean_sd$macroAvgPrecisionMean*100
  bench_mean_sd$macroAvgRecallMean <- bench_mean_sd$macroAvgRecallMean*100
  return(bench_mean_sd)
}
human<- SingleSpeciesRes("human","Human")
cattle <- SingleSpeciesRes("cattle","Cattle")
pig <- SingleSpeciesRes("pig","Pig")
all <- rbind(pig,cattle,human)
all$`Error rate` <- factor(all$`Error rate`,levels=sort(unique(all$`Error rate`),decreasing = T))

################## Figure 2a(left)£ºrss ################################
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\RSS.tiff",width = 4000,height =3640,res = 600)
ggplot(all, aes(SpeMap, weight = RSSMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Memory consumption (GB)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,400),breaks=seq(0,400,80)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",3),rep("indianred4",3),
                                             rep("hotpink3",3),rep("salmon3",3),
                                             rep("azure4",3),rep("tan4",3),
                                             rep("orchid4",3),rep("darkseagreen4",3),
                                             rep("lightskyblue4",3),rep("palevioletred4",3),rep("yellow3",3),
                                             rep("burlywood4",3),rep("tomato3",3),rep("antiquewhite4",3)),size = 11,
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        legend.title=element_text(colour ="black",face = "bold") , legend.text=element_text(colour ="black",face = "bold"),
        axis.text.x = element_text(size = 11,colour ="black",face = "bold"),
        axis.title.x = element_text(face = "bold"))
dev.off()
################## Figure 2a(left)£ºruntime ##########
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\runtime.tiff",width = 4000,height =3640,res = 600)
ggplot(all, aes(SpeMap, weight = cpuTimeMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Runtime (min)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,1300),breaks=seq(0,1300,260)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",3),rep("indianred4",3),
                                             rep("hotpink3",3),rep("salmon3",3),
                                             rep("azure4",3),rep("tan4",3),
                                             rep("orchid4",3),rep("darkseagreen4",3),
                                             rep("lightskyblue4",3),rep("palevioletred4",3),rep("yellow3",3),
                                             rep("burlywood4",3),rep("tomato3",3),rep("antiquewhite4",3)),size = 11,
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        legend.title=element_text(colour ="black",face = "bold") , legend.text=element_text(colour ="black",face = "bold"),
        axis.text.x = element_text(size = 11,colour ="black",face = "bold"),
        axis.title.x = element_text(face = "bold"))
dev.off()


################## Figure 3a(left)£ºuniquely mapped reads ###############################

all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\uniqueMap.tiff",width = 4000,height =3640,res = 600)
ggplot(all, aes(SpeMap, weight = matchedReadsRateMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Uniquely mapped read (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,500),breaks=seq(0,500,100)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",3),rep("indianred4",3),
                                             rep("hotpink3",3),rep("salmon3",3),
                                             rep("azure4",3),rep("tan4",3),
                                             rep("orchid4",3),rep("darkseagreen4",3),
                                             rep("lightskyblue4",3),rep("palevioletred4",3),rep("yellow3",3),
                                             rep("burlywood4",3),rep("tomato3",3),rep("antiquewhite4",3)),size = 11,
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        legend.title=element_text(colour ="black",face = "bold") , legend.text=element_text(colour ="black",face = "bold"),
        axis.text.x = element_text(size = 11,colour ="black",face = "bold"),
        axis.title.x = element_text(face = "bold"))
dev.off()

################## Figure S5a£ºmapped precisoin #################################################
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\mapPre.tiff",width = 4000,height =3640,res = 600)
ggplot(all, aes(SpeMap, weight = macroAvgPrecisionMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Mapped precision (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,500),breaks=seq(0,500,100)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",3),rep("indianred4",3),
                                             rep("hotpink3",3),rep("salmon3",3),
                                             rep("azure4",3),rep("tan4",3),
                                             rep("orchid4",3),rep("darkseagreen4",3),
                                             rep("lightskyblue4",3),rep("palevioletred4",3),rep("yellow3",3),
                                             rep("burlywood4",3),rep("tomato3",3),rep("antiquewhite4",3)),size = 11,
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        legend.title=element_text(colour ="black",face = "bold") , legend.text=element_text(colour ="black",face = "bold"),
        axis.text.x = element_text(size = 11,colour ="black",face = "bold"),
        axis.title.x = element_text(face = "bold"))
dev.off()

################## Figure S5b£ºrecall ##################
all$SpeMap <- factor(all$SpeMap,levels = unique(all$SpeMap[order(all$mapper)]))
tiff("figure\\recall.tiff",width = 4000,height =3640,res = 600)
ggplot(all, aes(SpeMap, weight = macroAvgRecallMean, fill =`Error rate`)) +
  geom_bar(width = .9, position = 'stack') +
  labs( x="",y = 'Recall (%)') +
  guides(fill=guide_legend(reverse=TRUE))+
  coord_flip(clip = "off", expand = FALSE)+
  scale_fill_manual(values = c( "paleturquoise3","pink3","darkseagreen3","khaki3","lightskyblue3"))+
  scale_y_continuous(expand = c(0,0),limits=c(0,500),breaks=seq(0,500,100)) +
  theme(axis.text.y = element_text(colour =c(rep("seagreen4",3),rep("indianred4",3),
                                             rep("hotpink3",3),rep("salmon3",3),
                                             rep("azure4",3),rep("tan4",3),
                                             rep("orchid4",3),rep("darkseagreen4",3),
                                             rep("lightskyblue4",3),rep("palevioletred4",3),rep("yellow3",3),
                                             rep("burlywood4",3),rep("tomato3",3),rep("antiquewhite4",3)),size = 11,
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
  fileName_mean <- aggregate(fileName[,c(6,7,8,15,20,21)],by=list(fileName$mapper,fileName$species),mean)
  names(fileName_mean)[1:8] <- c("mapper","species","macroAvgPrecisionMean","macroAvgRecallMean","macroAvgF1ScoreMean","RSSMean","cpuTimeMean","matchedReadsRateMean")
  fileName_sd <- aggregate(fileName[,c(6,7,8,15,20,21)],by=list(fileName$mapper,fileName$species),sd)
  names(fileName_sd)[3:8] <-c("macroAvgPrecisionSd","macroAvgRecallSd","macroAvgF1ScoreSd","RSSSd","cpuTimeSd","matchedReadsRateSd")
  fileName_mean_sd <- cbind(fileName_mean,fileName_sd[,3:8])
  fileName_mean_sd$mapper <- as.character(fileName_mean_sd$mapper)
  oldName <- c("batmeth2","bismarkbwt2","bismarkhis2","bsmap","bsseeker2bt","bsseeker2bt2end","bsseeker2bt2loc",
               "bsseeker2soap","bwameth","walt","hisat_3n","hisat_3n_repeat","bsbolt","abismal")
  newName <- c("Batmeth2","Bismark-bwt2-e2e","Bismark-his2","BSMAP","BSSeeker2-bwt","BSSeeker2-bwt2-e2e","BSSeeker2-bwt2-local",
               "BSSeeker2-soap2","Bwa-meth","Walt","Hisat_3n","Hisat_3n_repeat","BSBolt","Abismal")
  for(i in 1:14){
    fileName_mean_sd$mapper[which(fileName_mean_sd$mapper==oldName[i])] <- c(newName[i])
  }
  fileName_mean_sd$species <- species2
  names(fileName_mean_sd)[2] <- "Species"
  return(fileName_mean_sd)
}
human <- dataRehape("human","Human")
pig <- dataRehape("pig","Pig")
cattle <- dataRehape("cattle","Cattle")
dataFrame <- rbind(human,pig,cattle)
dataFrame$mapper <- factor(dataFrame$mapper,levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt","Batmeth2",
                                                       "BSSeeker2-bwt2-local","BSSeeker2-soap2","Walt","Bismark-bwt2-e2e","Abismal",
                                                       "BSMAP","BSBolt","Hisat_3n","Hisat_3n_repeat","Bwa-meth"))
dataFrame$matchedReadsRateMean <- dataFrame$matchedReadsRateMean*100
dataFrame$Species <- factor(dataFrame$Species,levels = c("Human","Cattle","Pig"))
################## Figure 2c: runtime #######

tiff("figure/runTimeAllErr.tiff",width = 4000,height =3700,res = 600)
p<-ggplot(dataFrame,aes(x=mapper,y=cpuTimeMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Runtime (min)')+
  scale_y_continuous(limits=c(0,300),breaks=seq(0,300,50)) +
  geom_point(size=4,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 17,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 18,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 17,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=18,colour ="black",face = "bold") , 
        legend.text=element_text(size=18,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
print(p)
dev.off()

################## Figure 2d: rss ##############
tiff("figure/rssAllErr.tiff",width = 4000,height =3700,res = 600)
p<-ggplot(dataFrame,aes(x=mapper,y=RSSMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Memory consumption (GB)')+
  scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20)) +
  geom_point(size=4,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 17,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 18,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 17,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=18,colour ="black",face = "bold") , 
        legend.text=element_text(size=18,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
print(p)
dev.off()

################## Figure 3b: uniquely mapped reads ######
tiff("figure/uniMapAllErr.tiff",width = 4000,height =3700,res = 600)
p<-ggplot(dataFrame,aes(x=mapper,y=matchedReadsRateMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Uniquely mapped reads (%)')+
  geom_point(size=4,aes(shape=Species,colour=Species)) + 
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
  theme(axis.title.y = element_text( size = 17,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 18,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 17,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=18,colour ="black",face = "bold") , 
        legend.text=element_text(size=18,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
print(p)
dev.off()
################## Figure 3d: precision and recall #####################

precisonRecall<- function(data,species){
  data$mapper <- factor(data$mapper,levels =c("Bwa-meth","Hisat_3n_repeat","Hisat_3n","BSBolt","BSMAP","Abismal","Bismark-bwt2-e2e","Walt",
                                              "BSSeeker2-soap2","BSSeeker2-bwt2-local","Batmeth2","BSSeeker2-bwt","Bismark-his2","BSSeeker2-bwt2-e2e"))
  colnames(data)[1] <- "Alignment algorithms"
  tiff(paste0("figure\\",species,"_precision_recall.tiff"),width = 5550,height =3500,res = 600)
  p<- ggplot(data = data, mapping = aes(x = macroAvgRecallMean*100, y = macroAvgPrecisionMean*100,shape=`Alignment algorithms`, colour = `Alignment algorithms`)) + 
    geom_point(size = 7,stroke = 2)+
    labs( x="Recall (%)",y = 'Mapped precision (%)') +
    geom_abline(slope=1,intercept = 0,linetype="dashed", size=0.4)+
    scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,11,12,13,14,15))+
    scale_y_continuous(expand = c(0,0),limits=c(70,101),breaks=seq(70,100,10)) +
    scale_x_continuous(expand = c(0,0),limits=c(70,101),breaks=seq(70,100,10)) +
    theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.line = element_line(color="black"),
          axis.title = element_text(colour ="black",face = "bold",size=20),
          axis.text = element_text(colour ="black",face = "bold",size=20),
          legend.text=element_text(colour ="black",face = "bold",size=20),
          legend.title=element_text(colour ="black",face = "bold",size=20))  
  print(p)
  dev.off()
}
precisonRecall(human,"human")
precisonRecall(cattle,"cattle")
precisonRecall(pig,"pig")

################## Figure 3e: F1 score #########################
tiff("figure/f1scoreAllErr.tiff",width = 4000,height =3700,res = 600)
p<-ggplot(dataFrame,aes(x=mapper,y=macroAvgF1ScoreMean*100,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'F1 Score (%)')+
  scale_y_continuous(limits=c(70,100),breaks=seq(70,100,10)) +
  geom_point(size=4,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 17,colour ="black",face = "bold"),
        axis.text.y = element_text( size = 18,colour ="black",face = "bold"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 17,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        legend.title=element_text(size=18,colour ="black",face = "bold") , 
        legend.text=element_text(size=18,colour ="black",face = "bold"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
print(p)
dev.off()
################## Figure 3f: uniquely mapped reads and F1 socre: average of five error rate and four species #########################
human <- read.csv(paste("metaData\\","human","BenchSimuDataAccuUniConbine.csv",sep = ""),header = T,sep = "\t", quote = "\"", dec = ".")
cattle <- read.csv(paste("metaData\\","cattle","BenchSimuDataAccuUniConbine.csv",sep = ""),header = T,sep = "\t", quote = "\"", dec = ".")
pig <- read.csv(paste("metaData\\","pig","BenchSimuDataAccuUniConbine.csv",sep = ""),header = T,sep = "\t", quote = "\"", dec = ".")
fileName <- rbind(human,cattle,pig)
fileName$cpuTime <- fileName$cpusysTime+fileName$cpuuserTime
fileName$RSS <- fileName$RSS/1000000
fileName$cpuTime <- fileName$cpuTime/60
fileName$matchedReadsRate <- fileName$matchedReads/2000000
#  fileName<- fileName[,c(1,2,3,6,15,20,21,19)]
fileName_mean <- aggregate(fileName[,c(6,7,8,15,20,21)],by=list(fileName$mapper),mean)
names(fileName_mean)[1:7] <- c("Alignment algorithms","macroAvgPrecisionMean","macroRecallnMean","macroF1ScoreMean","RSSMean","cpuTimeMean","matchedReadsRateMean")

fileName_sd <- aggregate(fileName[,c(6,7,8,15,20,21)],by=list(fileName$mapper),sd)
names(fileName_sd)[2:7] <- c("macroAvgPrecisionSd","macroRecallSd","macroF1ScoreSd","RSSSd","cpuTimeSd","matchedReadsRateSd")
fileName_mean_sd <- cbind(fileName_mean,fileName_sd[,2:7])
fileName_mean_sd$'Alignment algorithms' <- as.character(fileName_mean_sd$'Alignment algorithms')
oldName <- c("abismal","batmeth2","bismarkbwt2","bismarkhis2","bsbolt","bsmap","bsseeker2bt","bsseeker2bt2end","bsseeker2bt2loc",
             "bsseeker2soap","bwameth","hisat_3n","hisat_3n_repeat","walt")
newName <- c("Abismal","Batmeth2","Bismark-bwt2-e2e","Bismark-his2","BSBolt","BSMAP","BSSeeker2-bwt","BSSeeker2-bwt2-e2e","BSSeeker2-bwt2-local",
             "BSSeeker2-soap2","Bwa-meth","Hisat_3n","Hisat_3n_repeat","Walt")

for(i in 1:14){
  fileName_mean_sd$'Alignment algorithms'[which(fileName_mean_sd$'Alignment algorithms'==oldName[i])] <- c(newName[i])
}

fileName_mean_sd$'Alignment algorithms' <- factor(fileName_mean_sd$'Alignment algorithms',
                                                  levels =c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt","Batmeth2",
                                                            "BSSeeker2-bwt2-local","BSSeeker2-soap2","Walt","Bismark-bwt2-e2e","Abismal",
                                                            "BSMAP","BSBolt","Hisat_3n","Hisat_3n_repeat","Bwa-meth"))
fileName_mean_sd$'Alignment algorithms' <- factor(fileName_mean_sd$'Alignment algorithms',levels =c("Bwa-meth","Hisat_3n_repeat","Hisat_3n","BSBolt","BSMAP","Abismal","Bismark-bwt2-e2e",
                                                                                                    "Walt","BSSeeker2-soap2","BSSeeker2-bwt2-local","Batmeth2","BSSeeker2-bwt",
                                                                                                    "Bismark-his2","BSSeeker2-bwt2-e2e"))
all <- fileName_mean_sd
rm(human,cattle,pig,fileName_mean,fileName_sd,fileName,fileName_mean_sd)

tiff("figure\\allSpe_unique_f1.tiff",width = 5550,height =3500,res = 600)
ggplot(data = all, mapping = aes(x = macroF1ScoreMean*100, y = matchedReadsRateMean*100,shape=`Alignment algorithms`, colour = `Alignment algorithms`)) + 
  geom_point(size = 7,stroke = 2)+
  labs( x="F1 score (%)",y = 'Uniquely mapped reads (%)') +
  geom_abline(slope=1,intercept = 0,linetype="dashed", size=0.4)+
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,11,12,13,14,15))+
  scale_y_continuous(expand = c(0,0),limits=c(70,101),breaks=seq(70,100,10)) +
  scale_x_continuous(expand = c(0,0),limits=c(70,101),breaks=seq(70,100,10)) +
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title = element_text(colour ="black",face = "bold",size=20),
        axis.text = element_text(colour ="black",face = "bold",size=20),
        legend.text=element_text(colour ="black",face = "bold",size=20),
        legend.title=element_text(colour ="black",face = "bold",size=20))  
dev.off()
################## Figure 2b: cor: heatmap ###############
SingleSpeciesRes<- function(species,species2){
  bench <- read.csv(paste("metaData\\",species,"BenchSimuDataAccuUniConbine.csv",sep = ""),header = TRUE, sep = "\t", quote = "\"", dec = ".")
  bench$cpuTime<- bench$cpusysTime+bench$cpuuserTime
  bench$errorRate<- paste(bench$errorRate,"%",sep = "", collapse = NULL)
  bench$RSS <- bench$RSS/1000000
  bench$cpuTime <- bench$cpuTime/60
  bench$matchedReadsRate <- bench$matchedReads/2000000
  colnames(bench)
  bench_mean <- aggregate(bench[,c(6,7,9,15,20,21)],by=list(bench$errorRate,bench$mapper),mean)
  names(bench_mean)[1:8] <- c("errorRate","mapper","macroAvgPrecisionMean","macroAvgRecallMean","microAvgPrecisionMean","RSSMean","cpuTimeMean","matchedReadsRateMean")
  bench_sd <- aggregate(bench[,c(6,7,9,15,20,21)],by=list(bench$errorRate,bench$mapper),sd)
  names(bench_sd)[3:8] <- c("macroAvgPrecisionSd","macroAvgRecallSd","microAvgPrecisionSd","RSSSd","cpuTimeSd","matchedReadsRateSd")
  bench_mean_sd <- cbind(bench_mean,bench_sd[,3:8])
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
  bench_mean_sd[bench_mean_sd$mapper=="hisat_3n",2] <- "Hisat_3n"
  bench_mean_sd[bench_mean_sd$mapper=="hisat_3n_repeat",2] <- "Hisat_3n_repeat"
  bench_mean_sd[bench_mean_sd$mapper=="bsbolt",2] <- "BSBolt"
  bench_mean_sd[bench_mean_sd$mapper=="abismal",2] <- "Abismal"
  bench_mean_sd[bench_mean_sd$errorRate=="0%",1] <- "0"
  bench_mean_sd[bench_mean_sd$errorRate=="0.5%",1] <- "0.50%"
  bench_mean_sd[bench_mean_sd$errorRate=="1%",1] <- "1.00%"
  
  names(bench_mean_sd)[1] <- "Error rate"
  bench_mean_sd$mapper <- factor(bench_mean_sd$mapper,
                                 levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt","Batmeth2",
                                            "BSSeeker2-bwt2-local","BSSeeker2-soap2","Walt","Bismark-bwt2-e2e","Abismal",
                                            "BSMAP","BSBolt","Hisat_3n","Hisat_3n_repeat","Bwa-meth"))
  bench_mean_sd$SpeMap <- paste(species2,": ",bench_mean_sd$mapper,sep = "")
  bench_mean_sd$matchedReadsRateMean <- bench_mean_sd$matchedReadsRateMean*100
  bench_mean_sd$macroAvgPrecisionMean <- bench_mean_sd$macroAvgPrecisionMean*100
  bench_mean_sd$macroAvgRecallMean <- bench_mean_sd$macroAvgRecallMean*100
  return(bench_mean_sd)
}

human<- SingleSpeciesRes("human","Human")
cattle <- SingleSpeciesRes("cattle","Cattle")
pig <- SingleSpeciesRes("pig","Pig")
human$species <- "Human"
cattle$species <- "Cattle"
pig$species <- "Pig"
all <- rbind(pig,cattle,human)
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
cpu <- errCorClass(all,7,"CpuTime")
rss <- errCorClass(all,6,"RSS")
pre <- errCorClass(all,3,"Precision")
recall <- errCorClass(all,4,"Recall")
uni <- errCorClass(all,8,"Unique")

res <- cbind(cpu,rss$RSSCor,rss$RSSCorP,pre$PrecisionCor,pre$PrecisionCorP,recall$RecallCor,recall$RecallCorP,uni$UniqueCor,uni$UniqueCorP)
names(res)[5:12] <- c("RSScor","RSSCorP","PrecisionCor","PrecisionCorP","RecallCor","RecallCorP","UniqueCor","UniqueCorP")

#write.csv(res,file = "significance\\corTest.csv",row.names = F)

# plot
res$SpeMap <- paste(res$Species,res$Mapper,sep = ": ")
res[4,2]

res <- res[c(as.numeric(row.names(res[which(res$Mapper=="Bwa-meth"),])),
             as.numeric(row.names(res[which(res$Mapper=="Hisat_3n_repeat"),])),
             as.numeric(row.names(res[which(res$Mapper=="Hisat_3n"),])),
             as.numeric(row.names(res[which(res$Mapper=="BSBolt"),])),
             as.numeric(row.names(res[which(res$Mapper=="BSMAP"),])),
             as.numeric(row.names(res[which(res$Mapper=="Abismal"),])),
             as.numeric(row.names(res[which(res$Mapper=="Bismark-bwt2-e2e"),])),
             as.numeric(row.names(res[which(res$Mapper=="Walt"),])),
             as.numeric(row.names(res[which(res$Mapper=="BSSeeker2-soap2"),])),
             as.numeric(row.names(res[which(res$Mapper=="BSSeeker2-bwt2-local"),])),
             as.numeric(row.names(res[which(res$Mapper=="Batmeth2"),])),
             as.numeric(row.names(res[which(res$Mapper=="BSSeeker2-bwt"),])),
             as.numeric(row.names(res[which(res$Mapper=="Bismark-his2"),])),
             as.numeric(row.names(res[which(res$Mapper=="BSSeeker2-bwt2-e2e"),]))),]

corErr <- res[,c(3,5,11,7,9)]
row.names(corErr) <- res$SpeMap
names(corErr) <- c("Runtime","Memory consumption","Uniquely mapped reads","Mapped precision","Recall")
`Mapping algorithms`<- factor(c(rep("Bwa-meth",3),rep("Hisat_3n_repeat",3),rep("Hisat_3n",3),rep("BSBolt",3),rep("BSMAP",3),rep("Abismal",3),
                                rep("Bismark-bwt2-e2e",3),rep("Walt",3),rep("BSSeeker2-soap2",3),
                                rep("BSSeeker2-bwt2-local",3),rep("Batmeth2",3),rep("BSSeeker2-bwt",3),rep("Bismark-his2",3),rep("BSSeeker2-bwt2-e2e",3)),
                              levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt","Batmeth2",
                                         "BSSeeker2-bwt2-local","BSSeeker2-soap2","Walt","Bismark-bwt2-e2e","Abismal",
                                         "BSMAP","BSBolt","Hisat_3n","Hisat_3n_repeat","Bwa-meth"))
`Mapping algorithms` <- data.frame(`Mapping algorithms`)
names(`Mapping algorithms`) <- "Alignment algorithms"
rownames(`Mapping algorithms`)=rownames(corErr)
ann_colors = list(
  `Alignment algorithms`=c(`Bwa-meth`="antiquewhite4",
                           Hisat_3n_repeat="tomato3",
                           Hisat_3n="burlywood4",
                           BSBolt="yellow3",
                           BSMAP="palevioletred4",	
                           Abismal="lightskyblue4",
                           `Bismark-bwt2-e2e`="darkseagreen4",
                           Walt="orchid4",
                           `BSSeeker2-soap2`="tan4",
                           `BSSeeker2-bwt2-local`="azure4",
                           Batmeth2="salmon3",
                           `BSSeeker2-bwt`="hotpink3",
                           `Bismark-his2`="indianred4",
                           `BSSeeker2-bwt2-e2e`="seagreen4")
)
bk <- c(seq(-1,-0.001,by=0.0001),seq(0,1,by=0.0001))

tiff("figure\\heatCor.tiff",width = 4400,height =4040,res = 600)
pheatmap(corErr,annotation_row =`Mapping algorithms`,annotation_colors = ann_colors,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks = bk,
         legend_breaks = seq(-1,1,1),legend_labels = c("-1", "0", "1"),
         cluster_row = FALSE, cluster_col = FALSE
         #         ,filename = "figure/heatCot.tiff"
)
dev.off()

################## cor: Scatter plot ######################
all$`Error rate` <- as.numeric(all$`Error rate`)
all$species <- factor(all$species,levels = c("Human","Cattle","Pig"))
all$mapper <- factor(all$mapper,
                     levels = rev(c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt","Batmeth2",
                                    "BSSeeker2-bwt2-local","BSSeeker2-soap2","Walt","Bismark-bwt2-e2e","Abismal",
                                    "BSMAP","BSBolt","Hisat_3n","Hisat_3n_repeat","Bwa-meth")))
#runtime
################## Figure S2 #######################################
tiff("figure\\scatter_runtime.tiff",width = 9100,height =5700,res = 600)
ggplot(all,aes(x=`Error rate`,y=cpuTimeMean))+
  geom_point(aes(color=species))+
  labs( x="",y = '')+
  stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
  stat_cor(aes(color=species), method = "pearson",size=5)+
  theme(axis.title = element_text( size = 15,colour ="black",face = "bold"),
        axis.text = element_text( size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        strip.text = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
  facet_wrap( .~ mapper,ncol = 5,scales= "free" )
dev.off()

################## Figure S3 ########################################
# memory consumption
tiff("figure\\scatter_memory.tiff",width = 8900,height =5700,res = 600)
ggplot(all,aes(x=`Error rate`,y=RSSMean))+
  geom_point(aes(color=species))+
  labs( x="",y = '')+
  stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
  stat_cor(aes(color=species), method = "pearson",size=5)+
  theme(axis.title = element_text( size = 15,colour ="black",face = "bold"),
        axis.text = element_text( size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        strip.text = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
  facet_wrap( .~ mapper,ncol = 5,scales= "free" )
dev.off()

################## Figure S4 #############################
# uniquely mapped reads
tiff("figure\\scatter_unique.tiff",width = 9100,height =5700,res = 600)
ggplot(all,aes(x=`Error rate`,y=matchedReadsRateMean))+
  geom_point(aes(color=species))+
  labs( x="",y = '')+
  stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
  stat_cor(aes(color=species), method = "pearson",size=5)+
  theme(axis.title = element_text( size = 15,colour ="black",face = "bold"),
        axis.text = element_text( size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        strip.text = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
  facet_wrap( .~ mapper,ncol = 5,scales= "free" )
dev.off()

################## Figure S6 ####################
# mapped precision
tiff("figure\\scatter_precision.tiff",width = 9000,height =5700,res = 600)
ggplot(all,aes(x=`Error rate`,y=macroAvgPrecisionMean))+
  geom_point(aes(color=species))+
  labs( x="",y = '')+
  stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
  stat_cor(aes(color=species), method = "pearson",size=5)+
  theme(axis.title = element_text( size = 15,colour ="black",face = "bold"),
        axis.text = element_text( size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        strip.text = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
  facet_wrap( .~ mapper,ncol = 5,scales= "free" )
dev.off()

################## Figure S7 ##################################
# recall
tiff("figure\\scatter_recall.tiff",width = 9000,height =5700,res = 600)
ggplot(all,aes(x=`Error rate`,y=macroAvgRecallMean))+
  geom_point(aes(color=species))+
  labs( x="",y = '')+
  stat_smooth(aes(color=species),method = "lm",se=FALSE, formula = 'y ~ x')+
  stat_cor(aes(color=species), method = "pearson",size=5)+
  theme(axis.title = element_text( size = 15,colour ="black",face = "bold"),
        axis.text = element_text( size = 15,colour ="black",face = "bold"),
        axis.line = element_line(color="black"),
        strip.text = element_text(colour ="black",face = "bold",size = 13),
        legend.title=element_text(colour ="black",face = "bold",size = 13),
        legend.text=element_text(colour ="black",face = "bold",size = 13),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
  facet_wrap( .~ mapper,ncol = 5,scales= "free" )
dev.off()

