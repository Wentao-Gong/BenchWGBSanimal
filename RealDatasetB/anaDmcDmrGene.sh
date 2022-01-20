#!/bin/sh
resultDir=../meth
script=../myScript/realRes
annoFilePath=../annotation

speciesList=(human cattle pig)
for species in "${speciesList[@]}"
do
	## the analysis of DMC
	python ${script}/dmcAnalysis3.py \
		-a ${resultDir}/${species}2/dssRes2/bismarkbwt2_DML.txt \
		-b ${resultDir}/${species}2/dssRes2/bsmap_DML.txt \
		-c ${resultDir}/${species}2/dssRes2/bwameth_DML.txt \
		-d ${resultDir}/${species}2/dssRes2/walt_DML.txt \
		-o ${resultDir}/${species}2/dmcRes/ \
		-r ${annoFilePath}/${species}/${species}rmskchr.bed \
		-cg ${annoFilePath}/${species}/CGIandNonCGI.bed

	## the analysis of DMR
	python ${script}/dmrAna.py \
		-a ${resultDir}/${species}2/dssRes2/bismarkbwt2_DMR.txt \
		-b ${resultDir}/${species}2/dssRes2/bsmap_DMR.txt \
		-c ${resultDir}/${species}2/dssRes2/bwameth_DMR.txt \
		-d ${resultDir}/${species}2/dssRes2/walt_DMR.txt \
		-o ${resultDir}/${species}2/dmr/${species} \
		-r ${annoFilePath}/${species}/${species}rmskchr.bed \
		-cg ${annoFilePath}/${species}/CGIandNonCGI.bed

	## the analysis of DMR-related gene
	mapperList=(bismarkbwt2 bsmap bwameth walt)
	for mapper in "${mapperList[@]}"
	do
		python ${script}/textToBedDmr.py -i ${resultDir}/${species}2/dssRes2/${mapper}_DMR.txt -o ${resultDir}/${species}2/dssRes2/${mapper}_DMR.bed
		bedtools intersect -a ${annoFilePath}/${species}/${species}GeneCoo.bed -b ${resultDir}/${species}2/dssRes2/${mapper}_DMR.bed -wao > ${resultDir}/${species}2/gene/${mapper}_gene_DMR.txt
		python ${script}/geneRelaDmr.py -i ${resultDir}/${species}2/gene/${mapper}_gene_DMR.txt -o ${resultDir}/${species}2/gene/${mapper}_gene_DMR_Res.txt
	done
	python ${script}/geneVenn.py \
		-a ${resultDir}/${species}2/gene/bismarkbwt2_gene_DMR_Res.txt \
		-b ${resultDir}/${species}2/gene/bsmap_gene_DMR_Res.txt \
		-c ${resultDir}/${species}2/gene/bwameth_gene_DMR_Res.txt \
		-d ${resultDir}/${species}2/gene/walt_gene_DMR_Res.txt \
		-o ${resultDir}/${species}2/gene/gene_DMR_Venn_report.txt
	
	## the analysis of signaling pathway in KEGG
	python ${script}/keggVenn.py \
		-a ${resultDir}/${species}2/kegg/${species}_bismarkbwt2_enrichKEGG_depth10.txt \
		-b ${resultDir}/${species}2/kegg/${species}_bsmap_enrichKEGG_depth10.txt \
		-c ${resultDir}/${species}2/kegg/${species}_bwameth_enrichKEGG_depth10.txt \
		-d ${resultDir}/${species}2/kegg/${species}_walt_enrichKEGG_depth10.txt \
		-o ${resultDir}/${species}2/kegg/${species}_kegg_venn.txt
done

