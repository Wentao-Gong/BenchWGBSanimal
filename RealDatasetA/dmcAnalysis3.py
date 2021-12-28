# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 10:33:41 2021
处理的数据主要是基于四个比对软件的比对结果得到的DMC
一致性DMC：四个软件DMC结果的交集
非一致性DMC：四个软件整体全部的DMC-四个软件DMC结果的交集


@author: 滔
"""

import argparse
import os 
parser = argparse.ArgumentParser(description='Base quality analysis')
parser.add_argument('-a', '--afile', help='a file (bismarkbwt2)', required=True)
parser.add_argument('-b', '--bfile', help='b file (bsmap)', required=True)
parser.add_argument('-c', '--cfile', help='c file (bwameth)', required=True)
parser.add_argument('-d', '--dfile', help='d file (walt)', required=True)
parser.add_argument('-o', '--outpath', help='Enter the path to store the results', required=True)
parser.add_argument('-r', '--repanno', help='Input annotation file of repeat sequence', required=True)
parser.add_argument('-cg', '--cgianno', help='Input annotation file of CGI related area', required=True)
args = parser.parse_args()

aFile=args.afile
bFile=args.bfile
cFile=args.cfile
dFile=args.dfile
outPath=args.outpath
repAnnoFile = args.repanno
cgiAnnoFile = args.cgianno

### Statistical results of each data set

def resOut(dataSet,outputFileName):
    outputFile=open(outputFileName,"w")
    for line in dataSet:
        line=line.split("+")
        chrom=line[0]
        start=str(int(line[1])-1)
        end=str(int(line[1])+2)
        outputFile.write(chrom+"\t"+start+"\t"+end+"\n")

def intersecionOfFourDataSet(aFile,bFile,cFile,dFile,outPath):
    aSet=set()
    bSet=set()
    cSet=set()
    dSet=set()
    

    aFileName=os.path.basename(aFile).split('.')[0]
    bFileName=os.path.basename(bFile).split('.')[0]
    cFileName=os.path.basename(cFile).split('.')[0]
    dFileName=os.path.basename(dFile).split('.')[0]
#    fileName="/"+aFileName.split("_")[1]+"_"+aFileName.split("_")[2]
    fileName=aFileName.split("_")[1]
    
    comDmcFileName=outPath+fileName+"commonDmc.bed"
    allNonComDmcFileName=outPath+fileName+"allNonCommonDmc.bed"
    aDiffDmcFileName=outPath+"/"+aFileName+"_DiffDmc.bed"
    bDiffDmcFileName=outPath+"/"+bFileName+"_DiffDmc.bed"
    cDiffDmcFileName=outPath+"/"+cFileName+"_DiffDmc.bed"
    dDiffDmcFileName=outPath+"/"+dFileName+"_DiffDmc.bed"
    reportFileName=outPath+fileName+"report.txt"
    
    reportFile=open(reportFileName,"w")
    
    with open(aFile,"r") as f:
        for line in f:
            line=line.split()
            if "pos" not in line:
                aSet.add(line[0]+"+"+line[1])
    with open(bFile,"r") as f:
        for line in f:
            line=line.split()
            if "pos" not in line:
                bSet.add(line[0]+"+"+line[1])
    with open(cFile,"r") as f:
        for line in f:
            line=line.split()
            if "pos" not in line:
                cSet.add(line[0]+"+"+line[1])
    with open(dFile,"r") as f:
        for line in f:
            line=line.split()
            if "pos" not in line:
                dSet.add(line[0]+"+"+line[1])
    
    interSet = aSet & bSet & cSet & dSet
    allSet = aSet | bSet | cSet | dSet
    allDiffSet = allSet - interSet
    aDiffSet = aSet - interSet
    bDiffSet = bSet - interSet
    cDiffSet = cSet - interSet
    dDiffSet = dSet - interSet
    
    resOut(interSet,comDmcFileName)
    resOut(allDiffSet,allNonComDmcFileName)
    resOut(aDiffSet,aDiffDmcFileName)
    resOut(bDiffSet,bDiffDmcFileName)
    resOut(cDiffSet,cDiffDmcFileName)
    resOut(dDiffSet,dDiffDmcFileName)
    
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
    abcdNum = len(interSet)
    allNum = len(allSet)
    allDiffNum = len(allDiffSet)
    aDiffNum = len(aDiffSet)
    bDiffNum = len(bDiffSet)
    cDiffNum = len(cDiffSet)
    dDiffNum = len(dDiffSet)
    abcDiffNum = len(aDiffSet & bDiffSet & cDiffSet)
    abdDiffNum = len(aDiffSet & bDiffSet & dDiffSet)
    acdDiffNum = len(aDiffSet & cDiffSet & dDiffSet)
    bcdDiffNum = len(bDiffSet & cDiffSet & dDiffSet)

    reportFile.write("aFileName"+"\t"+aFileName+"\n"+\
                    "bFileName"+"\t"+bFileName+"\n"+\
                    "cFileName"+"\t"+cFileName+"\n"+\
                    "dFileName"+"\t"+dFileName+"\n"+\
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
                    "allNum"+"\t"+str(allNum)+"\n"+\
                    "allDiffNum"+"\t"+str(allDiffNum)+"\n"+\
                    "aDiffNum"+"\t"+str(aDiffNum)+"\n"+\
                    "bDiffNum"+"\t"+str(bDiffNum)+"\n"+\
                    "cDiffNum"+"\t"+str(cDiffNum)+"\n"+\
                    "dDiffNum"+"\t"+str(dDiffNum)+"\n"+
                    "abcDiffNum"+"\t"+str(abcDiffNum)+"\n"+\
                    "abdDiffNum"+"\t"+str(abdDiffNum)+"\n"+\
                    "acdDiffNum"+"\t"+str(acdDiffNum)+"\n"+\
                    "bcdDiffNum"+"\t"+str(bcdDiffNum)+"\n")
    resNumDict = {}
    resNumDict["abcdNum"]=abcdNum
    resNumDict["allDiffNum"]=allDiffNum
    resNumDict["aDiffNum"]=aDiffNum
    resNumDict["bDiffNum"]=bDiffNum
    resNumDict["cDiffNum"]=cDiffNum
    resNumDict["dDiffNum"]=dDiffNum
    return(resNumDict)

resNumDict=intersecionOfFourDataSet(aFile,bFile,cFile,dFile,outPath)

### region annotation
aFileName=os.path.basename(aFile).split('.')[0]
bFileName=os.path.basename(bFile).split('.')[0]
cFileName=os.path.basename(cFile).split('.')[0]
dFileName=os.path.basename(dFile).split('.')[0]
#fileName="/"+aFileName.split("_")[1]+"_"+aFileName.split("_")[2] 
fileName=aFileName.split("_")[1]

comDmcFileName=outPath+fileName+"commonDmc.bed"
allNonComDmcFileName=outPath+fileName+"allNonCommonDmc.bed"
aDiffDmcFileName=outPath+"/"+aFileName+"_DiffDmc.bed"
bDiffDmcFileName=outPath+"/"+bFileName+"_DiffDmc.bed"
cDiffDmcFileName=outPath+"/"+cFileName+"_DiffDmc.bed"
dDiffDmcFileName=outPath+"/"+dFileName+"_DiffDmc.bed"

def singleResStaRep(inFileName,allCpGNum,class1,class2):
    dictName={}
    dictName[class1]=0
    inFile=open(inFileName,"r")
    for line in inFile:
        resClass=line.split()[5]
        if resClass == "-1":
            dictName[class1] += 1
    dictName[class2] = allCpGNum -dictName[class1]
    return(dictName)
    
def allResStaRep(comDmcAnnoResCGI,allNonComDmAnnoResCGI,aDiffDmcAnnoResCGI,bDiffDmcAnnoResCGI,cDiffDmcAnnoResCGI,dDiffDmcAnnoResCGI,resNumDict,outFileName,class1,class2):
    outFile=open(outFileName,"w")
    
    comDict=singleResStaRep(comDmcAnnoResCGI,resNumDict["abcdNum"],class1,class2)
    allNonComDict=singleResStaRep(allNonComDmAnnoResCGI,resNumDict["allDiffNum"],class1,class2)
    aDiffDict=singleResStaRep(aDiffDmcAnnoResCGI,resNumDict["aDiffNum"],class1,class2)
    bDiffDict=singleResStaRep(bDiffDmcAnnoResCGI,resNumDict["bDiffNum"],class1,class2)
    cDiffDict=singleResStaRep(cDiffDmcAnnoResCGI,resNumDict["cDiffNum"],class1,class2)
    dDiffDict=singleResStaRep(dDiffDmcAnnoResCGI,resNumDict["dDiffNum"],class1,class2)
    
    comDmcAnnoResCGIFileName=os.path.basename(comDmcAnnoResCGI).split(".")[0]
    allNonComDmAnnoResCGIFileName=os.path.basename(allNonComDmAnnoResCGI).split(".")[0]
    aDiffDmcAnnoResCGIFileName=os.path.basename(aDiffDmcAnnoResCGI).split(".")[0]
    bDiffDmcAnnoResCGIFileName=os.path.basename(bDiffDmcAnnoResCGI).split(".")[0]
    cDiffDmcAnnoResCGIFileName=os.path.basename(cDiffDmcAnnoResCGI).split(".")[0]
    dDiffDmcAnnoResCGIFileName=os.path.basename(dDiffDmcAnnoResCGI).split(".")[0]
    
    outFile.write("resClass"+"\t"+comDmcAnnoResCGIFileName+"\t"+allNonComDmAnnoResCGIFileName+"\t"+\
                  aDiffDmcAnnoResCGIFileName+"\t"+bDiffDmcAnnoResCGIFileName+"\t"+\
                      cDiffDmcAnnoResCGIFileName+"\t"+dDiffDmcAnnoResCGIFileName+"\n")
    classSet=set(comDict.keys())
    for resClass in classSet:
        com=comDict[resClass]
        allNonCom=allNonComDict[resClass]
        aDiff=aDiffDict[resClass]
        bDiff=bDiffDict[resClass]
        cDiff=cDiffDict[resClass]
        dDiff=dDiffDict[resClass]
        outFile.write(resClass+"\t"+str(com)+"\t"+str(allNonCom)+"\t"+str(aDiff)+"\t"+str(bDiff)+"\t"+str(cDiff)+"\t"+str(dDiff)+"\n")
### CGI related region annotation

comDmcAnnoResCGI=outPath+fileName+"comDmcAnnoResCGI.bed"
allNonComDmAnnoResCGI=outPath+fileName+"allNonComDmcAnnoResCGI.bed"
aDiffDmcAnnoResCGI=outPath+"/"+aFileName+"_DiffDmcAnnoResCGI.bed"
bDiffDmcAnnoResCGI=outPath+"/"+bFileName+"_DiffDmcAnnoResCGI.bed"
cDiffDmcAnnoResCGI=outPath+"/"+cFileName+"_DiffDmcAnnoResCGI.bed"
dDiffDmcAnnoResCGI=outPath+"/"+dFileName+"_DiffDmcAnnoResCGI.bed"

os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(comDmcFileName,cgiAnnoFile,comDmcAnnoResCGI))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(allNonComDmcFileName,cgiAnnoFile,allNonComDmAnnoResCGI))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(aDiffDmcFileName,cgiAnnoFile,aDiffDmcAnnoResCGI))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(bDiffDmcFileName,cgiAnnoFile,bDiffDmcAnnoResCGI))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(cDiffDmcFileName,cgiAnnoFile,cDiffDmcAnnoResCGI))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(dDiffDmcFileName,cgiAnnoFile,dDiffDmcAnnoResCGI))

dmcAnnoResStaCGI=outPath+fileName+"dmcAnnoResStaCGI.txt"
allResStaRep(comDmcAnnoResCGI,allNonComDmAnnoResCGI,aDiffDmcAnnoResCGI,bDiffDmcAnnoResCGI,cDiffDmcAnnoResCGI,dDiffDmcAnnoResCGI,resNumDict,dmcAnnoResStaCGI,"nonCgi","cgi")

os.system("rm -r %s" % comDmcAnnoResCGI)
os.system("rm -r %s" % allNonComDmAnnoResCGI)
os.system("rm -r %s" % aDiffDmcAnnoResCGI)
os.system("rm -r %s" % bDiffDmcAnnoResCGI)
os.system("rm -r %s" % cDiffDmcAnnoResCGI)
os.system("rm -r %s" % dDiffDmcAnnoResCGI)

### repeat region annotation
comDmcAnnoResRep=outPath+fileName+"comDmcAnnoResRep.bed"
allNonComDmAnnoResRep=outPath+fileName+"allNonComDmcAnnoResRep.bed"
aDiffDmcAnnoResRep=outPath+"/"+aFileName+"_DiffDmcAnnoResRep.bed"
bDiffDmcAnnoResRep=outPath+"/"+bFileName+"_DiffDmcAnnoResRep.bed"
cDiffDmcAnnoResRep=outPath+"/"+cFileName+"_DiffDmcAnnoResRep.bed"
dDiffDmcAnnoResRep=outPath+"/"+dFileName+"_DiffDmcAnnoResRep.bed"

os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(comDmcFileName,repAnnoFile,comDmcAnnoResRep))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(allNonComDmcFileName,repAnnoFile,allNonComDmAnnoResRep))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(aDiffDmcFileName,repAnnoFile,aDiffDmcAnnoResRep))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(bDiffDmcFileName,repAnnoFile,bDiffDmcAnnoResRep))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(cDiffDmcFileName,repAnnoFile,cDiffDmcAnnoResRep))
os.system('/project/pub/software/bedtools2-2.25.0/bin/bedtools intersect -a %s -b %s -wao > %s' %(dDiffDmcFileName,repAnnoFile,dDiffDmcAnnoResRep))

dmcAnnoResStaRep=outPath+fileName+"dmcAnnoResStaRep.txt"
#cgiAndGeneAnnoResSta(comDmcAnnoResRep,allNonComDmAnnoResRep,aDiffDmcAnnoResRep,bDiffDmcAnnoResRep,cDiffDmcAnnoResRep,dDiffDmcAnnoResRep,dmcAnnoResStaRep,8)

allResStaRep(comDmcAnnoResRep,allNonComDmAnnoResRep,aDiffDmcAnnoResRep,bDiffDmcAnnoResRep,cDiffDmcAnnoResRep,dDiffDmcAnnoResRep,resNumDict,dmcAnnoResStaRep,"nonRepeat","repeat")


os.system("rm -r %s" % comDmcAnnoResRep)
os.system("rm -r %s" % allNonComDmAnnoResRep)
os.system("rm -r %s" % aDiffDmcAnnoResRep)
os.system("rm -r %s" % bDiffDmcAnnoResRep)
os.system("rm -r %s" % cDiffDmcAnnoResRep)
os.system("rm -r %s" % dDiffDmcAnnoResRep)

