#!/bin/sh

indexDir=../index
simudataDir=../data/simudate/depth5
softDir=../soft
resultDir=../result/depth5


speciesList=(human mouse cattle pig)
genomeList=(hg38 mm10 bosTau6 susScr11)
errorRateList=(0 1)
numList=(1 2 3)


for i in $(seq 0 3)
do
	for errorRate in "${errorRateList[@}"
	do
		echo "map reads to ref genome using bwameth"
		mapper=bwameth
		for num in "${numList[@]}"
		do
		mkdir -p ${resultDir}/${speciesList[$i]}${num}/${mapper}
			${softDir}/bwa-meth-master/bwameth.py \
				--reference ${indexDir}/${speciesList[$i]}/${mapper}/${genomeList[$i]}.fa \
				${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate}Depth5Num${num}_1.fastq \
				${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate}Depth5Num${num}_2.fastq \
				-t 8 \
				> ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam
			python /share/nas2/USER_DIR/wenping/gwentao/my_script/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}${num}/${mapper}/SimuDataAccuUni${errorRate}.csv
			samtools view -b -@ 8 ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam -o ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.bam
			rm ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam
		done

		echo "map reads to ref genome using bsmap"
		mapper=bsmap
		for num in "${numList[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}${num}/${mapper}
			${softDir}/bsmap-2.90/bsmap -d ${indexDir}/${speciesList[$i]}/${mapper}/${genomeList[$i]}.fa \
				-a ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate}Depth5Num${num}_1.fastq \
				-b ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate}Depth5Num${num}_2.fastq \
				-p 8 \
				-o ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam
			python /share/nas2/USER_DIR/wenping/gwentao/my_script/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}${num}/${mapper}/SimuDataAccuUni${errorRate}.csv
			samtools view -b -@ 8 ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam -o ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.bam
			rm ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam
		done

		echo "map reads to ref genome using walt"
		mapper=walt
		for num in "${numList[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}${num}/${mapper}
			${softDir}/walt-master/bin/walt -i ${indexDir}/${speciesList[$i]}/${mapper}/${genomeList[$i]}.dbindex \
				-t 8 \
				-sam \
				-a \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate}Depth5Num${num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate}Depth5Num${num}_2.fastq \
				-o ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam
			python /share/nas2/USER_DIR/wenping/gwentao/my_script/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}${num}/${mapper}/SimuDataAccuUni${errorRate}.csv
			samtools view -b -@ 8 ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam -o ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.bam
			rm ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam
		done


		echo "map reads to ref genome using bismarkbwt2"
		mapper=bismarkbwt2
		for num in "${numList[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}${num}/${mapper}
			${softDir}/Bismark-0.22.3/bismark --path_to_bowtie2 ${softDir}/bowtie2-2.3.5.1-linux-x86_64 \
				--genome ${indexDir}/${speciesList[$i]}/${mapper} \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate}Depth5Num${num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate}Depth5Num${num}_2.fastq \
				--sam \
				--ambiguous \
				${resultDir}/${speciesList[$i]}${num}/${mapper}/temp \
				-o ${resultDir}/${speciesList[$i]}${num}/${mapper}
			mv ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}_1_bismark_bt2_pe.sam \
				${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam
			python /share/nas2/USER_DIR/wenping/gwentao/my_script/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}${num}/${mapper}/SimuDataAccuUni${errorRate}.csv
			samtools view -b -@ 3 ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam -o ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.bam
			rm ${resultDir}/${speciesList[$i]}${num}/${mapper}/simulatedErrRates${errorRate}Depth5Num${num}.sam
		done
	done
done