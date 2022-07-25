#!/bin/sh

indexDir=../index
softDir=../soft
realDataDir=../data
resultDir=../result/realBench

MAPPERLIST=(bismarkbwt2 bismarkhis2 bsmap bwameth walt batmeth2 bsseeker2bt bsseeker2bt2end bsseeker2bt2loc bsseeker2soap hisat_3n hisat_3n_repeat bsbolt abismal) 
genome=bosTau9
species=cattle
readLen=150

sampleList=(SRR7528450 SRR7528456 SRR7528458)


echo "map reads to ref genome using walt"
mapper=walt
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/walt-master/bin/walt -i ${indexDir}/${species}/${mapper}/${genome}.dbindex \
			-t 1 \
			-sam \
			-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
			-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
			-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using bwameth"
mapper=bwameth
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/bwa-meth-master/bwameth.py \
			--reference ${indexDir}/${species}/${mapper}/${genome}.fa \
			${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
			${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
			-t 1 \
			> ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done


echo "map reads to ref genome using bismarkbwt2"
mapper=bismarkbwt2
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/Bismark-0.22.3/bismark --path_to_bowtie2 ${softDir}/bowtie2-2.3.5.1-linux-x86_64 \
			--genome ${indexDir}/${species}/${mapper} \
			-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
			-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
			--sam \
			--temp_dir ${resultDir}/${species}/${sample}/${mapper}/temp \
			-o ${resultDir}/${species}/${sample}/${mapper}
		mv ${resultDir}/${species}/${sample}/${mapper}/${sample}_clean_readLen150_1_bismark_bt2_pe.sam \
			${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using bsmap"
mapper=bsmap
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/bsmap-2.90/bsmap -d ${indexDir}/${species}/${mapper}/${genome}.fa \
			-a ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
			-b ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
			-p 1 \
			-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using batmeth2"
mapper=batmeth2
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/BatMeth2/bin/BatMeth2 align \
		-g ${indexDir}/${species}/${mapper}/${genome}.fa \
		-p 1 \
		-of SAM \
		-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-o realDataSample${sample}
	mv ./realDataSample${sample}.sam \
		${resultDir}/${species}/${sample}/${mapper}/
	mv realDataSample${sample}.run.log \
		${resultDir}/${species}/${sample}/${mapper}/
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done


echo "map reads to ref genome using bismarkhis2"
mapper=bismarkhis2
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/Bismark-0.22.3/bismark --hisat2 --genome ${indexDir}/${species}/${mapper} \
		--path_to_hisat2  ${softDir}/hisat2-2.1.0 \
		-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq  \
		--sam \
		--temp_dir ${resultDir}/${species}/${sample}/${mapper}/temp \
		-o ${resultDir}/${species}/${sample}/${mapper}
	mv ${resultDir}/${species}/${sample}/${mapper}/${sample}_clean_readLen150_1_bismark_hisat2_pe.sam \
		${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam		
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done


echo "map reads to ref genome using bsseeker2bt"
mapper=bsseeker2bt
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	LD_PRELOAD=/share/nas2/USER_DIR/wenping/soft/glibc-2.14/lib/libc-2.14.so \
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-align.py \
		-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-g ${genome}.fa \
		-d ${indexDir}/${species}/${mapper} \
		-f sam \
		--bt-p 1 \
		--aligner=bowtie \
		-p ${softDir}/bowtie-1.3.0-linux-x86_64 \
		--temp_dir=${resultDir}/${species}/${sample}/${mapper}/ \
		-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done


echo "map reads to ref genome using bsseeker2bt2end"
mapper=bsseeker2bt2end
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-align.py \
		-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-g ${genome}.fa \
		-d ${indexDir}/${species}/bsseeker2bt2 \
		-f sam \
		--bt2-p 1 \
		--bt2--end-to-end \
		--aligner=bowtie2 \
		-p ${softDir}/bowtie2-2.3.4.3-linux-x86_64 \
		--temp_dir=${resultDir}/${species}/${sample}/${mapper}/ \
		-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using bsseeker2bt2loc"
mapper=bsseeker2bt2loc
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-align.py \
		-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-g ${genome}.fa \
		-d ${indexDir}/${species}/bsseeker2bt2 \
		-f sam \
		--bt2-p 1 \
		--aligner=bowtie2 \
		-p ${softDir}/bowtie2-2.3.4.3-linux-x86_64 \
		--temp_dir=${resultDir}/${species}/${sample}/${mapper}/ \
		-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using bsseeker2soap"
mapper=bsseeker2soap
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-align.py \
		-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-g ${genome}.fa \
		-d ${indexDir}/${species}/${mapper} \
		-f sam \
		--soap-p 1 \
		--soap-r 1 \
		--aligner=soap \
		-p /share/nas2/genome/biosoft/soap/2.21 \
		--temp_dir=${resultDir}/${species}/${sample}/${mapper}/ \
		-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using hisat_3n"
mapper=hisat_3n
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/hisat-3n/hisat-3n -x ${indexDir}/${species}/${mapper}/${genome} \
		-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-S ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-p 1 \
		--base-change C,T
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using hisat_3n_repeat"
mapper=hisat_3n_repeat
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/hisat-3n/hisat-3n -x ${indexDir}/${species}/${mapper}/${genome} \
		-1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-S ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-p 1 \
		--base-change C,T \
		--repeat
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using bsbolt"
mapper=bsbolt
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	python -m bsbolt Align \
		-F1 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		-F2 ${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-DB ${indexDir}/${species}/${mapper}/ \
		-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t 1
	${softDir}/samtools-1.12/bin/samtools view -h -@ 3 \
		${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.bam \
		-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo "map reads to ref genome using abismal"
mapper=abismal
for sample in "${sampleList[@]}"
do
	mkdir -p ${resultDir}/${species}/${sample}/${mapper}
	echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${species}/${sample}/${mapper}/Bench.csv
	/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${species}/${sample}/${mapper}/Bench.csv -a \
	${softDir}/abismal-3.0.0/bin/abismal \
		-i ${indexDir}/${species}/${mapper}/${genome}.abismalidx \
		-o ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_1.fastq \
		${realDataDir}/${species}/seedReadLen/${sample}_clean_readLen150_2.fastq \
		-t 1
	python /share/nas2/USER_DIR/wenping/gwentao/my_script/RealDataUniMap.py \
		-i ${resultDir}/${species}/${sample}/${mapper}/realDataSample${sample}.sam \
		-t ${mapper} \
		-s ${species} \
		-d ${sample} \
		-l 0 \
		-r ${readLen} \
		-a 2000000 \
		-o ${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv
	awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
		${resultDir}/${species}/${sample}/${mapper}/Bench.csv \
		${resultDir}/${species}/${sample}/${mapper}/RealDataUniMap.csv \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv
	cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap1.csv \
		| sed -n 2p \
		> ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv
done

echo -e "tool\tspecies\tdataName\tseedLen\treadLen\tcountMatchReads\tuniMapRate\tmem\tRSS\trealTime\tcpusysTime\tcpuuserTime" > ${resultDir}/${species}/${species}BenchRealDataUniMap.csv

for sample in "${sampleList[@]}"
do
	for mapper in "${MAPPERLIST[@]}"
	do
		cat ${resultDir}/${species}/${sample}/${mapper}/BenchRealDataUniMap.csv \
		>> ${resultDir}/${species}/${species}BenchRealDataUniMap.csv
	done
done
