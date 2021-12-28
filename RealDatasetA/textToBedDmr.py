# -*- coding: utf-8 -*-
"""
Created on Thu May 20 21:55:09 2021



@author: 滔
"""

import argparse
parser = argparse.ArgumentParser(description='Convert DMR txt format file to bed file format')
parser.add_argument('-i', '--input', help='The file of DMR obtained by DSS', required=True)
parser.add_argument('-o', '--output', help='The bed file of DMR', required=True)
args = parser.parse_args()

inputFileName=args.input
outputFileName=args.output

def inputFileOfDss(inputFileName,outputFileName):
    inputFile=open(inputFileName,"r")
    outputFile=open(outputFileName,"w")
    for line in inputFile:
        if "start" in line:
            continue
        else:
            line=line.split()
            if "e" in line[1]:
                line[1]=str(int(float(line[1])))
            elif "e" in line[2]:
                line[2]=str(int(float(line[2])))
            a="\t".join(line)
            outputFile.write(a+"\n")
		
inputFileOfDss(inputFileName,outputFileName)
