#!/bin/sh

indexDir=../index
simudataDir=../data/simudate/depth5
softDir=../soft

speciesList=(human cattle pig)
readNumList=(53488102 41698540 41698540)

for i in $(seq 0 3)
do
	ERRORRATE=(0 1)
	num=(1 2 3)
	depth=5
	mkdir -p ${simudataDir}/${speciesList[$i]}
	cd ${simudataDir}/${speciesList[$i]}
	for errorRate in "${ERRORRATE[@]}"
	do
		for num in "${num[@]}"
		do
			${softDir}/sherman/Sherman -l 150 -n ${readNumList[$i]} -genome_folder ${indexDir}/${speciesList[$i]} -pe --conversion_rate 100 --error_rate ${errorRate}
			mv ./simulated_1.fastq ./simulatedErrRates${errorRate//./}Depth${depth}Num${num}_1.fastq      # change file name
			mv ./simulated_2.fastq ./simulatedErrRates${errorRate//./}Depth${depth}Num${num}_2.fastq
			sed -i 's/_R1$//' ./simulatedErrRates${errorRate//./}Depth${depth}Num${num}_1.fastq           # remove _Rx from read ID
			sed -i 's/_R2$//' ./simulatedErrRates${errorRate//./}Depth${depth}Num${num}_2.fastq
		done
	done
done
