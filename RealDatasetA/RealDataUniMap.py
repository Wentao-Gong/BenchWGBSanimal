from __future__ import division
import argparse

parser = argparse.ArgumentParser(description='The number of uniquely mapped read in Real Dataset A')
parser.add_argument('-i', '--input', help='Input sam file', required=True)
parser.add_argument('-t', '--tool', help='used mapper', required=True)
parser.add_argument('-s', '--species', help='species', required=True)
parser.add_argument('-d', '--data', help='data name', required=True)
parser.add_argument('-l', '--seedLen', help='seed length', required=True)
parser.add_argument('-r', '--readLen', help='read length', required=True)
parser.add_argument('-o', '--output', help='output csv file', required=True)
parser.add_argument('-a', '--allRead', help='num of all Reads', required=True)

args = parser.parse_args()

inFile = args.input
outFile = args.output
tool = args.tool
species = args.species
data = args.data
seedLen = args.seedLen
readLen = args.readLen
#allReads=0
allReads=args.allRead
countMatchReads = 0

with open(inFile) as samFile:
	for line in samFile:
		if line.startswith('@'):
			continue
		else:
#			allReads += 1
			bitFlag = int(line.split('\t')[1])          
			if (bitFlag == 99 or bitFlag == 67):
				countMatchReads += 1
			elif (bitFlag == 147 or bitFlag == 131):
				countMatchReads += 1
			elif (bitFlag == 83 or bitFlag == 115):
				countMatchReads += 1
			elif (bitFlag == 163 or bitFlag == 179):
				countMatchReads += 1
			else:
				continue

uniMapRate=countMatchReads/int(allReads)

with open(outFile,'w') as f:
	f.write("tool\tspecies\tdataName\tseedLen\treadLen\tcountMatchReads\tuniMapRate\n")
	f.write(tool+"\t"+species+"\t"+data+"\t"+str(seedLen)+"\t"+str(readLen)+"\t"+str(countMatchReads)+"\t"+str(uniMapRate)+"\n")