#!/bin/sh
indexDir=../index
softDir=../soft
realDataDir=../data
resultDir=../result/realRes
script=../my_script/realRes
dataPath=../data

genome=susScr11
species=pig
genomeLen=2501912388

sampleList=(SRR7812176 SRR7812178 SRR7812179 SRR7812199 SRR7812200 SRR7812210)

echo "map reads to ref genome using walt"
mapper=walt
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	${softDir}/walt-master/bin/walt -i ${indexDir}/${species}/${mapper}/${genome}.dbindex \
        -t 8 \
         -sam \
         -1 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_1.fastq \
         -2 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_2.fastq \
         -o ${resultDir}/${species}/${sample}/${mapper}/${sample}.sam \
         -a \
         -u
	samtools view -b -@ 8 ${resultDir}/${species}/${sample}/${mapper}/${sample}.sam -o ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}.sam
	${softDir}/samtools-1.12/samtools sort -@ 10 ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam -o ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
	${softDir}/samtools-1.12/samtools index ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam
	${softDir}/MethylDackel extract -@ 8 ${indexDir}/${species}/${genome}.fa ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam -o ${resultDir}/${species}/${mapper}_${sample}
	python ${script}/prepareDSSinput.py -i ${resultDir}/${species}/${mapper}_${sample}_CpG.bedGraph -o ${resultDir}/${species}/${mapper}_${sample}_CpGofDSS.txt
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
done

echo "map reads to ref genome using bwameth"
mapper=bwameth
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	${softDir}/bwa-meth-master/bwameth.py \
       --reference ${indexDir}/${species}/${mapper}/${genome}.fa \
       ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_1.fastq \
       ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_2.fastq  \
       -t 8 \
       | samtools view -b - > \
       ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam
  ${softDir}/samtools-1.12/samtools sort -@ 10 ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam -o ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
	${softDir}/samtools-1.12/samtools index ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam
	${softDir}/MethylDackel extract -@ 8 ${indexDir}/${species}/${genome}.fa ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam -o ${resultDir}/${species}/${mapper}_${sample}
	python ${script}/prepareDSSinput.py -i ${resultDir}/${species}/${mapper}_${sample}_CpG.bedGraph -o ${resultDir}/${species}/${mapper}_${sample}_CpGofDSS.txt
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
done

echo "map reads to ref genome using bismarkbwt2"
mapper=bismarkbwt2
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	${softDir}/Bismark-0.22.3/bismark --path_to_bowtie2 ${softDir}/bowtie2-2.3.5.1-linux-x86_64 \
				--genome ${indexDir}/${species}/${mapper} \
        -1 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_1.fastq \
        -2 ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_2.fastq \
        --multicore 8 \
				--temp_dir ${resultDir}/${species}/${sample}/${mapper}/temp \
        -o ${resultDir}/${species}/${sample}/${mapper}
  mv ${resultDir}/${species}/${sample}/${mapper}/${sample}_clean_1_bismark_bt2_pe.bam ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam
  ${softDir}/samtools-1.12/samtools sort -@ 10 ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam -o ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
	${softDir}/samtools-1.12/samtools index ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam
	${softDir}/MethylDackel extract -@ 8 ${indexDir}/${species}/${genome}.fa ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam -o ${resultDir}/${species}/${mapper}_${sample}
	python ${script}/prepareDSSinput.py -i ${resultDir}/${species}/${mapper}_${sample}_CpG.bedGraph -o ${resultDir}/${species}/${mapper}_${sample}_CpGofDSS.txt
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
done

echo "map reads to ref genome using bsmap"
mapper=bsmap
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	${softDir}/bsmap-2.90/bsmap -d ${indexDir}/${species}/${mapper}/${genome}.fa \
       -a ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_1.fastq \
       -b ${realDataDir}/${species}/${sample}/cleandata/${sample}_clean_2.fastq \
       -p 8 \
       -o ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam
  ${softDir}/samtools-1.12/samtools sort -@ 10 ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam -o ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
	${softDir}/samtools-1.12/samtools index ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}.bam
	${softDir}/MethylDackel extract -@ 8 ${indexDir}/${species}/${genome}.fa ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam -o ${resultDir}/${species}/${mapper}_${sample}
	python ${script}/prepareDSSinput.py -i ${resultDir}/${species}/${mapper}_${sample}_CpG.bedGraph -o ${resultDir}/${species}/${mapper}_${sample}_CpGofDSS.txt
	rm ${resultDir}/${species}/${sample}/${mapper}/${sample}_sort.bam
done




