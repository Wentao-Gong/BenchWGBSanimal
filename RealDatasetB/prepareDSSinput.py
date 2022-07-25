import argparse
parser = argparse.ArgumentParser(description='Prepare the required input files for DSS')
parser.add_argument('-i', '--input', help='Prepare the result file obtained by MethylDackel', required=True)
parser.add_argument('-o', '--output', help='The required input files for DSS', required=True)
args = parser.parse_args()

inputFileName=args.input
outputFileName=args.output

def inputFileOfDss(inputFileName,outputFileName):
	inputFile=open(inputFileName,"r")
	outputFile=open(outputFileName,"w")
	outputFile.write("chr"+"\t"+"pos"+"\t"+"N"+"\t"+"X"+"\n")
	for line in inputFile:
		if line.startswith("track"):
			continue
		else:
			line=line.split()
			chrom=line[0]
			pos=line[2]
			N=str(int(line[4])+int(line[5]))
			X=line[4]
			outputFile.write(chrom+"\t"+pos+"\t"+N+"\t"+X+"\n")
			
inputFileOfDss(inputFileName,outputFileName)
