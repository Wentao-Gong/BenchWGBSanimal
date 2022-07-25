import argparse

parser = argparse.ArgumentParser(description='Convert the format of  CGI annnotation file to bed file format')
parser.add_argument('-i', '--input', help='Input a CGI annotataion file,end by .txt', required=True)
parser.add_argument('-o', '--output', help='The output file name,end by .bed', required=True)
args = parser.parse_args()

inFileName = args.input
outFileName = args.output

def toBed(inFileName,outFileName):
	inFile=open(inFileName,"r")
	outFile=open(outFileName,"w")
	for line in inFile:
		line=line.split()
		chrom=line[1]
		start=int(line[2])
		end=line[3]
		outFile.write(chrom+"\t"+str(start)+"\t"+end+"\n")

toBed(inFileName,outFileName)
