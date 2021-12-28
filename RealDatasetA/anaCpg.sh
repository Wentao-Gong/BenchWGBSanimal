#!/bin/sh
indexDir=../index
resultDir=../result/realRes
script=../my_script/realRes
annoFilePath=../annotation

species=human
python ${script}/CGIRegionDivision.py \
	-i ${annoFilePath}/${species}/cpgIslandExtUnmasked.txt \
	-o ${annoFilePath}/${species}/CGIRegionDivision.bed \
	-s hg38.chrom.sizes.txt
cat ${annoFilePath}/${species}/rmsk.txt |awk '{print $6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$17"\t"($8-$7)}'  > ${annoFilePath}/${species]}/${species}rmskchr.bed
sampleList=(SRR6373923 SRR6825466 SRR6825471 SRR6818517 SRR6373926 SRR6373932)
depth=10
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample}/
	python ${script}/conAndDisCpG4.py \
		-a ${resultDir}/${species}/bismarkbwt2_${sample}_CpG.bedGraph \
		-b ${resultDir}/${species}/bsmap_${sample}_CpG.bedGraph \
		-c ${resultDir}/${species}/bwameth_${sample}_CpG.bedGraph \
		-d ${resultDir}/${species}/walt_${sample}_CpG.bedGraph \
		-t ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample} \
		-o ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample}/dpeth${depth}_${sample} \
		-de ${depth} \
		-r ${annoFilePath}/${species}/${species}rmskchr.bed \
		-cg ${annoFilePath}/${species}/CGIRegionDivision.bed 
done

species=mouse
python ${script}/CGIRegionDivision.py \
	-i ${annoFilePath}/${species}/cpgIslandExtUnmasked.txt \
	-o ${annoFilePath}/${species}/CGIRegionDivision.bed \
	-s mm10.chrom.sizes.txt
cat ${annoFilePath}/${species]}/rmsk.txt |awk '{print $6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$17"\t"($8-$7)}'  > ${annoFilePath}/${species]}/${species}rmskchr.bed
sampleList=(SRR13482504 SRR13482503 SRR13482507 SRR13482509 SRR13482506 SRR13482505)
depth=10
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample}/
	python ${script}/conAndDisCpG4.py \
		-a ${resultDir}/${species}/bismarkbwt2_${sample}_CpG.bedGraph \
		-b ${resultDir}/${species}/bsmap_${sample}_CpG.bedGraph \
		-c ${resultDir}/${species}/bwameth_${sample}_CpG.bedGraph \
		-d ${resultDir}/${species}/walt_${sample}_CpG.bedGraph \
		-t ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample} \
		-o ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample}/dpeth${depth}_${sample} \
		-de ${depth} \
		-r ${annoFilePath}/${species}/${species}rmskchr.bed \
		-cg ${annoFilePath}/${species}/CGIRegionDivision.bed 
done

species=cattle
python ${script}/CGIRegionDivision.py \
	-i ${annoFilePath}/${species}/cpgIslandExtUnmasked.txt \
	-o ${annoFilePath}/${species}/CGIRegionDivision.bed \
	-s bosTau6.chrom.sizes.txt
cat ${annoFilePath}/${species]}/rmsk.txt |awk '{print $6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$17"\t"($8-$7)}'  > ${annoFilePath}/${species]}/${species}rmskchr.bed
sampleList=(SRR7528450 SRR7528456 SRR7528458 SRR7528459 SRR7528464 SRR7528465)
depth=10
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample}/
	python ${script}/conAndDisCpG4.py \
		-a ${resultDir}/${species}/bismarkbwt2_${sample}_CpG.bedGraph \
		-b ${resultDir}/${species}/bsmap_${sample}_CpG.bedGraph \
		-c ${resultDir}/${species}/bwameth_${sample}_CpG.bedGraph \
		-d ${resultDir}/${species}/walt_${sample}_CpG.bedGraph \
		-t ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample} \
		-o ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample}/dpeth${depth}_${sample} \
		-de ${depth} \
		-r ${annoFilePath}/${species}/${species}rmskchr.bed \
		-cg ${annoFilePath}/${species}/CGIRegionDivision.bed 
done

species=pig
python ${script}/CGIRegionDivision.py \
	-i ${annoFilePath}/${species}/cpgIslandExtUnmasked.txt \
	-o ${annoFilePath}/${species}/CGIRegionDivision.bed \
	-s susScr11.chrom.sizes.txt
cat ${annoFilePath}/${species]}/rmsk.txt |awk '{print $6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$17"\t"($8-$7)}'  > ${annoFilePath}/${species]}/${species}rmskchr.bed
sampleList=(SRR7812176 SRR7812178 SRR7812179 SRR7812199 SRR7812200 SRR7812210)
depth=10
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample}/
	python ${script}/conAndDisCpG4.py \
		-a ${resultDir}/${species}/bismarkbwt2_${sample}_CpG.bedGraph \
		-b ${resultDir}/${species}/bsmap_${sample}_CpG.bedGraph \
		-c ${resultDir}/${species}/bwameth_${sample}_CpG.bedGraph \
		-d ${resultDir}/${species}/walt_${sample}_CpG.bedGraph \
		-t ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample} \
		-o ${resultDir}/${species}/cpgAnaRes/depth${depth}_${sample}/dpeth${depth}_${sample} \
		-de ${depth} \
		-r ${annoFilePath}/${species}/${species}rmskchr.bed \
		-cg ${annoFilePath}/${species}/CGIRegionDivision.bed 
done


