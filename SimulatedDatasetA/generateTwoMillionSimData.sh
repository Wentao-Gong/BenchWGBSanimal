#!/bin/sh

#***
#date:   9.2020
#task:   generate simulated bisulfite reads
#***

indexDir=../index
simudataDir=../data/simudate
softDir=../soft

speciesList=(human mouse cattle pig)
ERRORRATE=(0 0.25 0.5 0.75 1)
sampleDuplicationList=(1 2 3)

mkdir ${simudataDir}/${species}
cd ${simudataDir}/${species}
for species in "${speciesList[@]}"
do
	for errorRate in "${ERRORRATE[@]}"
	do
		for num in "${sampleDuplicationList[@]}"
		do
			${softDir}/sherman/Sherman -l 150 -n 1000000 -genome_folder ${indexDir}/${species} -pe --conversion_rate 100 --error_rate ${errorRate}
			mv ./simulated_1.fastq ./simulatedErrRates${errorRate//./}Num${num}_1.fastq
			mv ./simulated_2.fastq ./simulatedErrRates${errorRate//./}Num${num}_2.fastq
			sed -i 's/_R1$//' ./simulatedErrRates${errorRate//./}Num${num}_1.fastq
			sed -i 's/_R2$//' ./simulatedErrRates${errorRate//./}Num${num}_2.fastq
		done
	done
done


