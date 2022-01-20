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
			python ${script}/readid2.py \
				-a ${dataPath}/${speciesList[$i]}${Num}/bismarkbwt2/simulatedErrRates${error}Depth5Num${Num}.bam \
				-b ${dataPath}/${speciesList[$i]}${Num}/bsmap/simulatedErrRates${error}Depth5Num${Num}.bam \
				-c ${dataPath}/${speciesList[$i]}${Num}/bwameth/simulatedErrRates${error}Depth5Num${Num}.bam \
				-d ${dataPath}/${speciesList[$i]}${Num}/walt/simulatedErrRates${error}Depth5Num${Num}.bam \
				-an bismarkbwt2 \
				-bn bsmap \
				-cn bwameth \
				-dn walt \
				-f ${mapDataPath}/${speciesList[$i]}/simulatedErrRates${error}Depth5Num${Num}_1.fastq \
				-bis ${dataPath}/${speciesList[$i]}${Num}/bismarkbwt2/simulatedErrRates${error}Depth5Num${Num}_1.fastq_ambiguous_reads_1.fq.gz \
				-s ${annoFilePath}/${speciesList[$i]}/${genomeList[$i]}.chrom.sizes.txt \
				-e ${dataPath}/${speciesList[$i]}${Num}/result/abnormalReport.txt \
				-o ${dataPath}/${speciesList[$i]}${Num}/result

			python ${script}/readAnnoAnalysis.py \
				-a ${dataPath}/${speciesList[$i]}${Num}/result/temp/allReadId.txt \
				-b ${dataPath}/${speciesList[$i]}${Num}/result/readId/commonNonRightMapReadId.txt \
				-c ${dataPath}/${speciesList[$i]}${Num}/result/temp/bismarkbwt2NonRightMapReadId.txt \
				-d ${dataPath}/${speciesList[$i]}${Num}/result/temp/bsmapNonRightMapReadId.txt \
				-e ${dataPath}/${speciesList[$i]}${Num}/result/temp/bwamethNonRightMapReadId.txt \
				-f ${dataPath}/${speciesList[$i]}${Num}/result/temp/waltNonRightMapReadId.txt \
				-r ${annoFilePath}/${speciesList[$i]}/${speciesList[$i]}rmskchr.bed \
				-cg ${annoFilePath}/${speciesList[$i]}/CGIandNonCGI.bed \
				-o ${dataPath}/${speciesList[$i]}${Num}/result/readAna
			rm ${dataPath}/${speciesList[$i]}${Num}/result/readAna/*bed
		done
	done
done