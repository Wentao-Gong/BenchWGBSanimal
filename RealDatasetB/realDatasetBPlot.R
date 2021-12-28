setwd("..\\WGBS\\result\\realBench")
library(ggplot2)

SingleSpeciesRes<- function(species,species2){
  bench <- read.csv(paste(species,"BenchRealDataUniMap.csv",sep = ""),header = TRUE, sep = "\t", quote = "\"", dec = ".")
  bench$cpuTime<- bench$cpusysTime+bench$cpuuserTime
  bench$RSS <- bench$RSS/1000000
  bench$cpuTime <- bench$cpuTime/60
  
  bench_mean <- aggregate(bench[,c(7,9,13)],by=list(bench$species,bench$tool),mean)
  names(bench_mean)[1:5] <- c("Species","Mapper","matchedReadsRateMean","RSSMean","cpuTimeMean")
  
  bench_sd <- aggregate(bench[,c(7,9,13)],by=list(bench$species,bench$tool),sd)
  names(bench_sd)[3:5] <- c("matchedReadsRateSd","RSSSd","cpuTimeSd")
  bench_mean_sd <- cbind(bench_mean,bench_sd[,3:5])
  bench_mean_sd$Mapper <- as.character(bench_mean_sd$Mapper)
  bench_mean_sd[bench_mean_sd$Mapper=="batmeth2",2] <- "Batmeth2"
  bench_mean_sd[bench_mean_sd$Mapper=="bismarkbwt2",2] <- "Bismark-bwt2-e2e"
  bench_mean_sd[bench_mean_sd$Mapper=="bismarkhis2",2] <- "Bismark-his2"
  bench_mean_sd[bench_mean_sd$Mapper=="bsmap",2] <- "BSMAP"
  bench_mean_sd[bench_mean_sd$Mapper=="bsseeker2bt",2] <- "BSSeeker2-bwt"
  bench_mean_sd[bench_mean_sd$Mapper=="bsseeker2bt2end",2] <- "BSSeeker2-bwt2-e2e"
  bench_mean_sd[bench_mean_sd$Mapper=="bsseeker2bt2loc",2] <- "BSSeeker2-bwt2-local"
  bench_mean_sd[bench_mean_sd$Mapper=="bsseeker2soap",2] <- "BSSeeker2-soap"
  bench_mean_sd[bench_mean_sd$Mapper=="bwameth",2] <- "Bwa-meth"
  bench_mean_sd[bench_mean_sd$Mapper=="walt",2] <- "Walt"
  bench_mean_sd$Mapper <- factor(bench_mean_sd$Mapper,levels = c("Batmeth2","Bismark-bwt2-e2e","Bismark-his2","BSMAP","BSSeeker2-bwt","BSSeeker2-bwt2-e2e","BSSeeker2-bwt2-local","BSSeeker2-soap",
                                                                 "Bwa-meth","Walt"))
  bench_mean_sd$SpeMap <- paste(species2,": ",bench_mean_sd$Mapper,sep = "")
  bench_mean_sd$matchedReadsRateMean <- bench_mean_sd$matchedReadsRateMean*100
  bench_mean_sd$Species <- species2
  return(bench_mean_sd) 
}
human<- SingleSpeciesRes("human","Human")
mouse <-SingleSpeciesRes("mouse","Mouse")
cattle <- SingleSpeciesRes("cattle","Cattle")
pig <- SingleSpeciesRes("pig","Pig")
all <- rbind(pig,cattle,mouse,human)

################ uniquely mapped reads ######################
all$Mapper <- factor(all$Mapper,levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt", "Batmeth2","BSSeeker2-bwt2-local",
                                           "BSSeeker2-soap","Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$Species <- factor(all$Species,levels = c("Human","Mouse","Cattle","Pig"))

tiff("figure\\realUniqueMapRateV3.tiff",width = 2100,height =1700,res = 300)
ggplot(all,aes(x=Mapper,y=matchedReadsRateMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Uniquely mapped reads (%)')+
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) +
  theme(axis.title.y = element_text( size = 15,colour ="black"),
        axis.text.y = element_text( size = 15,colour ="black"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 15,colour ="black"),
        axis.line = element_line(color="black"),
#        axis.ticks = element_blank(),
        legend.title=element_text(size=14,colour ="black") , legend.text=element_text(size=13,colour ="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())

dev.off()

################ rss ##################################
all$Mapper <- factor(all$Mapper,levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt", "Batmeth2","BSSeeker2-bwt2-local",
                                           "BSSeeker2-soap","Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$Species <- factor(all$Species,levels = c("Human","Mouse","Cattle","Pig"))
tiff("figure\\realRssV3.tiff",width = 2100,height =1700,res = 300)
ggplot(all,aes(x=Mapper,y=RSSMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Maximum resident set size (GB)')+
  scale_y_continuous(limits=c(0,80),breaks=seq(0,80,20)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black"),
        axis.text.y = element_text( size = 15,colour ="black"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 15,colour ="black"),
        axis.line = element_line(color="black"),
#        axis.ticks = element_blank(),
        legend.title=element_text(size=14,colour ="black") , legend.text=element_text(size=13,colour ="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()

################ cpu time ##################################
all$Mapper <- factor(all$Mapper,levels = c("BSSeeker2-bwt2-e2e","Bismark-his2","BSSeeker2-bwt", "Batmeth2","BSSeeker2-bwt2-local",
                                           "BSSeeker2-soap","Walt","Bismark-bwt2-e2e","BSMAP","Bwa-meth"))
all$Species <- factor(all$Species,levels = c("Human","Mouse","Cattle","Pig"))

tiff("figure\\realCputimeV3.tiff",width = 2100,height =1700,res = 300)
ggplot(all,aes(x=Mapper,y=cpuTimeMean,group=Species))+
  geom_line(aes(colour=Species)) +
  labs( x="",y = 'Cpu time (min)')+
  scale_y_continuous(limits=c(0,300),breaks=seq(0,300,50)) +
  geom_point(size=3,aes(shape=Species,colour=Species)) + 
  theme(axis.title.y = element_text( size = 15,colour ="black"),
        axis.text.y = element_text( size = 15,colour ="black"),
        axis.text.x = element_text(angle = 62, hjust = 1, size = 15,colour ="black"),
        axis.line = element_line(color="black"),
#        axis.ticks = element_blank(),
        legend.title=element_text(size=14,colour ="black") , legend.text=element_text(size=13,colour ="black"),
        panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()