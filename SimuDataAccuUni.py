import argparse

parser = argparse.ArgumentParser(description='Caculated the uniquely mapped reads, mapped precision, recall and F1 score in Simulated Dataset A')
parser.add_argument('-i', '--input', help='Input sam file', required=True)
parser.add_argument('-t', '--tool', help='used mapper', required=True)
parser.add_argument('-s', '--species', help='species', required=True)
parser.add_argument('-e', '--errRate', help='error rate', required=True)
parser.add_argument('-l', '--seedLen', help='seed length', required=True)
parser.add_argument('-r', '--readLen', help='read length', required=True)
parser.add_argument('-o', '--output', help='output csv file', required=True)

args = parser.parse_args()

inFile = args.input
outFile = args.output
tool = args.tool
species = args.species
errorRate = args.errRate
seedLen = args.seedLen
readLen = args.readLen
countMatchReads = 0

falseNegative = {}
falsePositive = {}
truePositive = {}
positions = set()

with open(inFile) as samFile:
    for line in samFile:
        if line.startswith('@'):
            continue
        elif line.split('\t')[0].split(':')[1].split('-')[0].split('_')[0]=='':
            continue
        else:
            bitFlag = int(line.split('\t')[1])
            predChrom = line.split('\t')[2]
            trueChrom = line.split('\t')[0].split(':')[0].split('_',1)[1]
            
			#first read forward strand (or case [67]: bsseeker2)
            if (bitFlag == 99 or bitFlag == 67):
                predPosition = int(line.split('\t')[3])
                truePosition = int(line.split('\t')[0].split(':')[1].split('-')[0].split('_')[0])
			#secound read reverse strand (or case [131]: bsseeker2)
            elif (bitFlag == 147 or bitFlag == 131):
                predPosition = int(line.split('\t')[3])
                truePosition = int(line.split('\t')[0].split(':')[1].split('-')[1].split('_')[0])-int(readLen)+1
			#first read reverse strand
            elif (bitFlag == 83 or bitFlag == 115):
                predPosition = int(line.split('\t')[3])
                truePosition = int(line.split('\t')[0].split(':')[1].split('-')[1].split('_')[0])-int(readLen)+1
			#secound read forward strand
            elif (bitFlag == 163 or bitFlag == 179):
                predPosition = int(line.split('\t')[3])
                truePosition = int(line.split('\t')[0].split(':')[1].split('-')[0].split('_')[0])  
			#read doesn't match
            else:
                continue
            
            positions.add(trueChrom + '_' + str(truePosition))
            positions.add(predChrom + '_' + str(predPosition))
            countMatchReads += 1

			#update count dictionary
            if predPosition == truePosition:
                if  predChrom + '_' + str(predPosition) in truePositive :
                    truePositive[predChrom + '_' + str(predPosition)] += 1
                else:
                    truePositive[predChrom + '_' + str(predPosition)] = 1
            else:
                if  trueChrom + '_' + str(truePosition) in falseNegative :
                    falseNegative[trueChrom + '_' + str(truePosition)] += 1
                else:
                    falseNegative[trueChrom + '_' + str(truePosition)] = 1
                
                if  predChrom + '_' + str(predPosition) in falsePositive :
                    falsePositive[predChrom + '_' + str(predPosition)] += 1
                else:
                    falsePositive[predChrom + '_' + str(predPosition)] = 1

#macro average
macroAvgPrecision = 0.0
macroAvgRecall = 0.0

#micro average
tmpSumTP = 0.0
tmpSumTPFP = 0.0
tmpSumTPFN = 0.0

avgAccuracy = 0.0
countPositions = len(positions)

for position in positions:

	FN = float(falseNegative[position]) if position in falseNegative else 0.0
	FP = float(falsePositive[position]) if position in falsePositive else 0.0
	TP = float(truePositive[position]) if position in truePositive else 0.0
	TN = countMatchReads - FN - TP - FP

	if (TP + FP != 0.0) and (TP + FN != 0.0):
		#macro average
		macroAvgPrecision = macroAvgPrecision + (TP / (TP + FP))
		macroAvgRecall = macroAvgRecall + (TP / (TP + FN))
		#micro average
		tmpSumTP = tmpSumTP + TP
		tmpSumTPFN = tmpSumTPFN + (TP + FN)
		tmpSumTPFP = tmpSumTPFP + (TP + FP)

		avgAccuracy = avgAccuracy + ((TP+TN) / countMatchReads)

#macro F1-Score
macroAvgPrecision = macroAvgPrecision / countPositions
macroAvgRecall = macroAvgRecall / countPositions
macroF1Score = 2* (macroAvgPrecision * macroAvgRecall) / (macroAvgPrecision + macroAvgRecall)

#micro F1-Score
microAvgPrecision = tmpSumTP / tmpSumTPFP
microAvgRecall = tmpSumTP / tmpSumTPFN
microF1Score = 2* (microAvgPrecision * microAvgRecall) / (microAvgPrecision + microAvgRecall)

avgAccuracy = avgAccuracy / countPositions

with open(outFile,"w") as fp:
    fp.write("mapper\tspecies\terrorRate\tseedLength\treadLength\tmacroAvgPrecision\tmacroAvgRecall\tmacroF1Score\tmicroAvgPrecision\tmicroAvgRecall\tmicroF1Score\tavgAccuracy\tmatchedReads\n")
    fp.write(tool+"\t"+species+"\t"+errorRate+"\t"+seedLen+"\t"+readLen+"\t"+str(macroAvgPrecision)+"\t"+str(macroAvgRecall)+"\t"+str(macroF1Score)+"\t"+str(microAvgPrecision)+"\t"+str(microAvgRecall)+"\t"+str(microF1Score)+"\t"+str(avgAccuracy)+"\t"+str(countMatchReads)+"\n") 

