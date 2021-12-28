#!/bin/sh
soft=../soft
realDataDir=../data

species=human
sampleList=(SRR6818517 SRR6373926 SRR6373932)
for sample in "${sampleList[@]}"
do
	${soft}/seqtk-master/seqtk sample -s100 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_1.fastq 1000000 > ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq
	${soft}/seqtk-master/seqtk sample -s100 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_2.fastq 1000000 > ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq
done 

species=mouse
sampleList=(SRR13482506 SRR13482507 SRR13482509)
for sample in "${sampleList[@]}"
do
	${soft}/seqtk-master/seqtk sample -s100 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_1.fastq 1000000 > ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq
	${soft}/seqtk-master/seqtk sample -s100 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_2.fastq 1000000 > ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq
done

species=cattle
sampleList=(SRR7528450 SRR7528456 SRR7528458)
for sample in "${sampleList[@]}"
do
	${soft}/seqtk-master/seqtk sample -s100 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_1.fastq 1000000 > ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq
	${soft}/seqtk-master/seqtk sample -s100 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_2.fastq 1000000 > ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq
done

species=pig
sampleList=(SRR7812176 SRR7812178 SRR7812179)
for sample in "${sampleList[@]}"
do
	${soft}/seqtk-master/seqtk sample -s100 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_1.fastq 1000000 > ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq
	${soft}/seqtk-master/seqtk sample -s100 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_2.fastq 1000000 > ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq
done 



 
