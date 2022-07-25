
import argparse
parser = argparse.ArgumentParser(description='The statistics of DMR-related genes')
parser.add_argument('-i', '--input', help='The file contain the overlap between DMR and Gene', required=True)
parser.add_argument('-o', '--output', help='The file contain the gene that was related to DMR', required=True)

args = parser.parse_args()

inputFileName=args.input
outputFileName=args.output

def geneRelaDmr(inputFileName,outputFileName):
    inputFile=open(inputFileName,"r")
    outputFile=open(outputFileName,"w")
    resDict = {}
    for line in inputFile:
        line=line.split()
        if line[7]!="-1":
            if line[3] in resDict:
                continue
            else:
                lineChr=line[0]
                lineStart=line[1]
                lineEnd=line[2]
                lineId=line[3]
                lineName=line[4]
                lineName2=line[5]
                resDict[line[3]]=lineName2+"\n"
                resDict[line[3]]=lineChr+"\t"+lineStart+"\t"+lineEnd+"\t"+lineId+"\t"+lineName+"\t"+lineName2+"\n"
    for key in resDict.keys():
        outputFile.write(resDict[key])

geneRelaDmr(inputFileName,outputFileName)

