#!/bin/sh
indexDir=../index
softDir=../soft

speciesList=(human cattle pig)
genomeList=(hg38 bosTau9 susScr11)

for i in $(seq 0 3)
do
	MAPPERLIST=(bismarkbwt2 bismarkhis2 bsmap bwameth walt batmeth2 bsseeker2bt bsseeker2bt2 bsseeker2soap hisat_3n hisat_3n_repeat bsbolt abismal)
	#creat folder of index
	for mapper in "${MAPPERLIST[@]}"
	do
		mkdir -p ${indexDir}/${speciesList[$i]}/${mapper}
		ln -s ${indexDir}/${speciesList[$i]}/${genomeList[$i]}.fa ${indexDir}/${speciesList[$i]}/${mapper}
	done

	echo "generate bismarkbwt2 index"
	${softDir}/Bismark-0.22.3/bismark_genome_preparation \
            --bowtie2 --path_to_aligner ${softDir}/bowtie2-2.3.5.1-linux-x86_64 \
             ${indexDir}/${speciesList[$i]}/bismarkbwt2

	echo "genereate bismarkhis2 index"
	${softDir}/Bismark-0.22.3/bismark_genome_preparation \
             --hisat2 --path_to_aligner ${softDir}/hisat2-2.1.0 \
             ${indexDir}/${speciesList[$i]}/bismarkhis2

	echo "generate bwameth index"
	${softDir}/bwa-meth-master/bwameth.py index \
             ${indexDir}/${speciesList[$i]}/bwameth/${genomeList[$i]}.fa

	echo "generate walt index"
	${softDir}/walt-master/bin/makedb \
             -c ${indexDir}/${speciesList[$i]}/walt/${genomeList[$i]}.fa \
             -o ${indexDir}/${speciesList[$i]}/walt/${genomeList[$i]}.dbindex

	echo "generate batmeth2 index"      
	${softDir}/BatMeth2/bin/BatMeth2 build_index \
            ${indexDir}/${speciesList[$i]}/batmeth2/${genomeList[$i]}.fa

	echo "generate bsseeker2bt index"
	LD_PRELOAD=${softDir}/glibc-2.14/lib/libc-2.14.so \
	python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-build.py \
      -f ${indexDir}/${speciesList[$i]}/bsseeker2bt/${genomeList[$i]}.fa \
      --aligner=bowtie \
      -p ${softDir}/bowtie-1.3.0-linux-x86_64 \
      -d ${indexDir}/${speciesList[$i]}/bsseeker2bt
      
	echo "generate bsseeker2bt2 index"
	python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-build.py \
      -f ${indexDir}/${speciesList[$i]}/bsseeker2bt2/${genomeList[$i]}.fa \
      --aligner=bowtie2 \
      -p ${softDir}/bowtie2-2.3.4.3-linux-x86_64 \
      -d ${indexDir}/${speciesList[$i]}/bsseeker2bt2
      
	echo "generate bsseeker2soap index"
	python ${softDir}/BSseeker2-BSseeker2-v2.1.8/bs_seeker2-build.py \
      -f ${indexDir}/${speciesList[$i]}/bsseeker2soap/${genomeList[$i]}.fa  \
      --aligner=soap \
      -p ${softDir}/soap/2.21 \
      -d ${indexDir}/${speciesList[$i]}/bsseeker2soap
  
  echo "generate hisat_3n index"    
  ${softDir}/hisat-3n/hisat-3n-build \
      --base-change C,T \
      ${indexDir}/${speciesList[$i]}/hista_3n/${genomeList[$i]}.fa \
      ${indexDir}/${speciesList[$i]}/hista_3n/${genomeList[$i]}
  
  echo "generate hisat_3n_repeat index"
  ${softDir}/hisat-3n/hisat-3n-build \
      --base-change T,C \
      --repeat-index ${indexDir}/${speciesList[$i]}/hista_3n_repeat/${genomeList[$i]}.fa \
      ${indexDir}/${speciesList[$i]}/hista_3n_repeat/${genomeList[$i]}
  
  echo "generate bsbolt index"
  python3 -m bsbolt \
      Index -G ${indexDir}/${speciesList[$i]}/bsbolt/${genomeList[$i]}.fa \
      -DB ${indexDir}/${speciesList[$i]}/bsbolt/
  
  echo "generate abismal index"
  ${softDir}/abismal-3.0.0/bin/abismalidx \
      ${indexDir}/${speciesList[$i]}/abismal/${genomeList[$i]}.fa \
      ${indexDir}/${speciesList[$i]}/abismal/${genomeList[$i]}.abismalidx
done