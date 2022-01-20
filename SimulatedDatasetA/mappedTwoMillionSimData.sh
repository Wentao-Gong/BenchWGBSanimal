#!/bin/sh

indexDir=../index
simudataDir=../data/simudate
softDir=../soft
resultDir=../result/bench
scriptDir=../my_script


speciesList=(human cattle pig)
genomeList=(hg38 bosTau9 susScr11)

for i in $(seq 0 3)
do
	MAPPERLIST=(bismarkbwt2 bismarkhis2 bsmap bwameth walt batmeth2 bsseeker2bt bsseeker2bt2end bsseeker2bt2loc bsseeker2soap) 
	ERRORRATE=(0 0.25 0.5 0.75 1)
	sampleDuplicationList=(1 2 3)

	for Num in "${sampleDuplicationList[@]}"
	do
		mkdir -p ${resultDir}/${speciesList[$i]}
		echo "map reads to ref genome using walt"
		mapper=walt
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			${softDir}/walt-master/bin/walt -i ${indexDir}/${speciesList[$i]}/${mapper}/${genomeList[$i]}.dbindex \
				-t 1 \
				-sam \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam          
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done

		echo "map reads to ref genome using batmeth2"
		mapper=batmeth2
		for errorRate in "${ERRORRATE[@]}"
		do	
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			${softDir}/BatMeth2/bin/BatMeth2 align \
				-g ${indexDir}/${speciesList[$i]}/${mapper}/${genomeList[$i]}.fa \
				-p 1 \
				-of SAM \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				-o simulatedErrRates${errorRate//./}Num${Num}
				mv ./simulatedErrRates${errorRate//./}Num${Num}.sam \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/
			mv ./simulatedErrRates${errorRate//./}Num${Num}.run.log \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/	          
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done

		echo "map reads to ref genome using bwameth"
		mapper=bwameth
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			${softDir}/bwa-meth-master/bwameth.py \
				--reference ${indexDir}/${speciesList[$i]}/${mapper}/${genomeList[$i]}.fa \
				${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				-t 1 \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done

		echo "map reads to ref genome using bismarkbwt2"
		mapper=bismarkbwt2
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			${softDir}/Bismark-0.22.3/bismark --path_to_bowtie2 ${softDir}/bowtie2-2.3.5.1-linux-x86_64 \
				--genome ${indexDir}/${speciesList[$i]}/${mapper} \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				--sam \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/temp \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			mv ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}_1_bismark_bt2_pe.sam \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done

		echo "map reads to ref genome using bismarkhis2"
		mapper=bismarkhis2
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			${softDir}/Bismark-0.22.3/bismark --hisat2 --genome ${indexDir}/${speciesList[$i]}/${mapper} \
				--path_to_hisat2  ${softDir}/hisat2-2.1.0 \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq  \
				--sam \
				--temp_dir ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/temp \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper} 
			mv ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}_1_bismark_hisat2_pe.sam \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam		
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done


		echo "map reads to ref genome using bsseeker2bt"
		mapper=bsseeker2bt
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			LD_PRELOAD=${softDir}/glibc-2.14/lib/libc-2.14.so \
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-align.py \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				-g ${genomeList[$i]}.fa \
				-d ${indexDir}/${speciesList[$i]}/${mapper} \
				-f sam \
				--bt-p 1 \
				--aligner=bowtie \
				-p ${softDir}/bowtie-1.3.0-linux-x86_64 \
				--temp_dir=${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/ \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done


		echo "map reads to ref genome using bsseeker2bt2end"
		mapper=bsseeker2bt2end
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-align.py \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				-g ${genomeList[$i]}.fa \
				-d ${indexDir}/${speciesList[$i]}/bsseeker2bt2 \
				-f sam \
				--bt2-p 1 \
				--bt2--end-to-end \
				--aligner=bowtie2 \
				-p ${softDir}/bowtie2-2.3.4.3-linux-x86_64 \
				--temp_dir=${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/ \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done

		echo "map reads to ref genome using bsseeker2bt2loc"
		mapper=bsseeker2bt2loc
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-align.py \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				-g ${genomeList[$i]}.fa \
				-d ${indexDir}/${speciesList[$i]}/bsseeker2bt2 \
				-f sam \
				--bt2-p 1 \
				--aligner=bowtie2 \
				-p ${softDir}/bowtie2-2.3.4.3-linux-x86_64 \
				--temp_dir=${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/ \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done

		echo "map reads to ref genome using bsseeker2soap"
		mapper=bsseeker2soap
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-align.py \
				-1 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-2 ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				-g ${genomeList[$i]}.fa \
				-d ${indexDir}/${speciesList[$i]}/${mapper} \
				-f sam \
				--soap-p 1 \
				--soap-r 1 \
				--aligner=soap \
				-p ${softDir}/soap/2.21 \
				--temp_dir=${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/ \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done

		echo "map reads to ref genome using bsmap"
		mapper=bsmap
		for errorRate in "${ERRORRATE[@]}"
		do
			mkdir -p ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}
			echo -e "mem\tRSS\trealTime\tcpusysTime\tcpuuserTime" >	${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv
			/usr/bin/time -f "%K\t%M\t%E\t%S\t%U" -o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv -a \
			${softDir}/bsmap-2.90/bsmap -d ${indexDir}/${speciesList[$i]}/${mapper}/${genomeList[$i]}.fa \
				-a ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_1.fastq \
				-b ${simudataDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}_2.fastq \
				-p 1 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam  	
			python ${scriptDir}/SimuDataAccuUni.py \
				-i ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/simulatedErrRates${errorRate//./}Num${Num}.sam \
				-t ${mapper} \
				-s ${speciesList[$i]} \
				-e ${errorRate} \
				-l 0 \
				-r 150 \
				-o ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv
			awk '{if(NR==FNR){a[FNR]=$0;}else{print $0 "\t" a[FNR]}}' \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/Bench.csv \
				${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/SimuDataAccuUni.csv \
				> ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapper}/BenchSimuDataAccuUni.csv
		done

		echo -e "mapper\tspecies\terrorRate\tseedLength\treadLength\tmacroAvgPrecision\tmacroAvgRecall\tmacroF1Score\tmicroAvgPrecision\tmicroAvgRecall\tmicroF1Score\tavgAccuracy\tmatchedReads\tmem\tRSS\trealTime\tcpusysTime\tcpuuserTime" > ${resultDir}/${speciesList[$i]}/${speciesList[$i]}BenchSimuDataAccuUniConbine${Num}.csv
		for errorRate in "${ERRORRATE[@]}"
		do
			for mapperlist in "${MAPPERLIST[@]}"
			do
				cat ${resultDir}/${speciesList[$i]}/simulatedErrRates${errorRate//./}Num${Num}/${mapperlist}/BenchSimuDataAccuUni.csv \
				| sed -n 2p \
				>> ${resultDir}/${speciesList[$i]}/${speciesList[$i]}BenchSimuDataAccuUniConbine${Num}.csv
			done
		done
	done
done