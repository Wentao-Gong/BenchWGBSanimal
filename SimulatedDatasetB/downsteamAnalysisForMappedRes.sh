#!/bin/sh

script=../my_script/unmapReadAnalysis
dataPath=../result/depth5
mapDataPath=../data/simudate/depth5
annoFilePath=../annotation
indexPath=../index


NumList=(1 2 3)
errorList=(0 1)
speciesList=(human cattle pig)
genomeList=(hg38 bosTau9 susScr11)


for species in "${errorList[@]}"
do
	python ${script}/CGIannoFileToBed.py -i ${annoFilePath}/${species]}/cpgIslandExtUnmasked.txt -o ${annoFilePath}/${species]}/CGIandNonCGI.bed
	cat ${annoFilePath}/${species]}/rmsk.txt |awk '{print $6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$17"\t"($8-$7)}'  > ${annoFilePath}/${species]}/rmsk.bed
	awk '{split($0,a,"chr");print a[2]}' ${annoFilePath}/${species]}/rmsk.bed > ${annoFilePath}/${species]}/rmsk.bed.new
	rm ${annoFilePath}/${species]}/rmsk.bed
	mv ${annoFilePath}/${species]}/rmsk.bed.new ${annoFilePath}/${species]}/${species]}rmsk.bed
done



for i in $(seq 0 3)
do
	for error in "${errorList[@]}"
	do
		for Num in "${NumList[@]}"
		do
			mkdir -p ${dataPath}/${speciesList[$i]}${Num}/result/readAna
			python ${script}/readID.py \
				-a ${dataPath}/${species}${Num}/bismarkbwt2/simulatedErrRates${error}Depth5Num${Num}.bam \
				-b ${dataPath}/${species}${Num}/bsmap/simulatedErrRates${error}Depth5Num${Num}.bam \
        -c ${dataPath}/${species}${Num}/bwameth/simulatedErrRates${error}Depth5Num${Num}.bam \
        -d ${dataPath}/${species}${Num}/walt/simulatedErrRates${error}Depth5Num${Num}.bam \
        -e ${dataPath}/${species}${Num}/bsbolt/simulatedErrRates${error}Depth5Num${Num}.bam \
				-an bismarkbwt2 \
				-bn bsmap \
				-cn bwameth \
				-dn walt \
				-en bsbolt \
				-f ${mapDataPath}/${speciesList[$i]}/simulatedErrRates${error}Depth5Num${Num}_1.fastq \
				-bis ${dataPath}/${speciesList[$i]}${Num}/bismarkbwt2/simulatedErrRates${error}Depth5Num${Num}_1.fastq_ambiguous_reads_1.fq.gz \
				-s ${annoFilePath}/${speciesList[$i]}/${genomeList[$i]}.chrom.sizes.txt \
				-er ${dataPath}/${speciesList[$i]}${Num}/result/abnormalReport.txt \
				-o ${dataPath}/${speciesList[$i]}${Num}/result
		done
	done
done