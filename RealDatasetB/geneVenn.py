# -*- coding: utf-8 -*-


import argparse
import os
parser = argparse.ArgumentParser(description='The venn results of DMR-related gene')
parser.add_argument('-a', '--afile', help='a file', required=True)
parser.add_argument('-b', '--bfile', help='b file', required=True)
parser.add_argument('-c', '--cfile', help='c file', required=True)
parser.add_argument('-d', '--dfile', help='d file', required=True)
parser.add_argument('-o', '--output', help='The file name of result', required=True)

args = parser.parse_args()
aFile = args.afile
bFile = args.bfile
cFile = args.cfile
dFile = args.dfile
outputFile =args.output

def geneVenn(aFile,bFile,cFile,dFile,outputFile):
    output=open(outputFile,"w")
    aMapperName=os.path.basename(aFile).split('_')[0]
    bMapperName=os.path.basename(bFile).split('_')[0]
    cMapperName=os.path.basename(cFile).split('_')[0]
    dMapperName=os.path.basename(dFile).split('_')[0]
    aSet=set()
    bSet=set()
    cSet=set()
    dSet=set()
    
    with open(aFile,"r") as f:
        for line in f:
            line=line.split()
            aSet.add(line[0])
    with open(bFile,"r") as f:
        for line in f:
            line=line.split()
            bSet.add(line[0])
    with open(cFile,"r") as f:
        for line in f:
            line=line.split()
            cSet.add(line[0])
    with open(dFile,"r") as f:
        for line in f:
            line=line.split()
            dSet.add(line[0])
    
    aNum=len(aSet)
    bNum=len(bSet)
    cNum=len(cSet)
    dNum=len(dSet)
    abNum = len(aSet & bSet)
    acNum = len(aSet & cSet)
    adNum = len(aSet & dSet)
    bcNum = len(bSet & cSet)
    bdNum = len(bSet & dSet)
    cdNum = len(cSet & dSet)
    abcNum = len(aSet & bSet & cSet)
    abdNum = len(aSet & bSet & dSet)
    acdNum = len(aSet & cSet & dSet)
    bcdNum = len(bSet & cSet & dSet)
    abcdNum = len(aSet & bSet & cSet & dSet)
    allNum = len(aSet | bSet | cSet | dSet)
    output.write("aFileName"+"\t"+aMapperName+"\n"+\
                    "bFileName"+"\t"+bMapperName+"\n"+\
                    "cFileName"+"\t"+cMapperName+"\n"+\
                    "dFileName"+"\t"+dMapperName+"\n"+\
                    "aNum"+"\t"+str(aNum)+"\n"+\
                    "bNum"+"\t"+str(bNum)+"\n"+\
                    "cNum"+"\t"+str(cNum)+"\n"+\
                    "dNum"+"\t"+str(dNum)+"\n"+\
                    "abNum"+"\t"+str(abNum)+"\n"+\
                    "acNum"+"\t"+str(acNum)+"\n"+\
                    "adNum"+"\t"+str(adNum)+"\n"+\
                    "bcNum"+"\t"+str(bcNum)+"\n"+\
                    "bdNum"+"\t"+str(bdNum)+"\n"+\
                    "cdNum"+"\t"+str(cdNum)+"\n"+\
                    "abcNum"+"\t"+str(abcNum)+"\n"+\
                    "abdNum"+"\t"+str(abdNum)+"\n"+\
                    "acdNum"+"\t"+str(acdNum)+"\n"+\
                    "bcdNum"+"\t"+str(bcdNum)+"\n"+\
                    "abcdNum"+"\t"+str(abcdNum)+"\n"+\
                    "allNum"+"\t"+str(allNum)+"\n")
    output.close()

geneVenn(aFile,bFile,cFile,dFile,outputFile)
    
    