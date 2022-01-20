

from __future__ import division
import os
import argparse
parser = argparse.ArgumentParser(description='The analysis of Concordant and Discordant CpG sites')
parser.add_argument('-a', '--afile', help='a file', required=True)
parser.add_argument('-b', '--bfile', help='b file', required=True)
parser.add_argument('-c', '--cfile', help='c file', required=True)
parser.add_argument('-d', '--dfile', help='d file', required=True)
parser.add_argument('-t', '--temppath', help='temp path', required=True)
parser.add_argument('-o', '--output', help='Result file output prefix', required=True)
parser.add_argument('-de', '--cpgdepth', help='Minimum depth of cpg site', required=True)
parser.add_argument('-r', '--repanno', help='Input annotation file of repeat sequence', required=True)
parser.add_argument('-cg', '--cgianno', help='Input annotation file of CGI related area', required=True)

args = parser.parse_args()
aFile = args.afile
bFile = args.bfile
cFile = args.cfile
dFile = args.dfile
tempPath = args.temppath
outPut =args.output
cpgDepth=args.cpgdepth
repAnnoFile = args.repanno
cgiAnnoFile = args.cgianno


## According to the chromosome number, cut large files into small files
fileList=[]

def cutFilesByChrom(infile,fileList,tempPath,cpgDepth):
    with open(infile,"r") as f:
        for line in f:
            if "track" not in line:                
                if int(line.split()[4])+int(line.split()[5]) >= int(cpgDepth):
                    chrom = line.split()[0]                
                    if chrom[3:].isdigit():
                        if chrom in fileList:
                            filename=''.join((tempPath,"/",chrom,".txt"))
                            with open(filename,"a") as df:
                                df.write(line)
                        else:
                            fileList.append(chrom)
                            filename=''.join((tempPath,"/",chrom,".txt"))
                            with open(filename,"w") as df:
                                df .write(line)
                    else:
                        if "other" in fileList:
                            filename=''.join((tempPath,"/other.txt"))
                            with open(filename,"a") as df:
                                df.write(line)
                        else:
                            fileList.append("other")
                            filename=''.join((tempPath,"/other.txt"))
                            with open(filename,"w") as df:
                                df.write(line)
    return(fileList)

atempfilename=''.join((tempPath,"/atemp"))
if not os.path.exists(atempfilename):
	os.mkdir(atempfilename)
fileList=cutFilesByChrom(aFile,fileList,atempfilename,cpgDepth)

btempfilename=''.join((tempPath,"/btemp"))
if not os.path.exists(btempfilename):
	os.mkdir(btempfilename)
fileList=cutFilesByChrom(bFile,fileList,btempfilename,cpgDepth)

ctempfilename=''.join((tempPath,"/ctemp"))
if not os.path.exists(ctempfilename):
	os.mkdir(ctempfilename)
fileList=cutFilesByChrom(cFile,fileList,ctempfilename,cpgDepth)

dtempfilename=''.join((tempPath,"/dtemp"))
if not os.path.exists(dtempfilename):
	os.mkdir(dtempfilename)
fileList=cutFilesByChrom(dFile,fileList,dtempfilename,cpgDepth)

## Screening consistent CpG sites and statistical results
def resOut(fileList,atempfilename,btempfilename,ctempfilename,dtempfilename,outPut,aFile,bFile,cFile,dFile):
    aNum=0
    bNum=0
    cNum=0
    dNum=0
    abNum=0
    acNum=0
    adNum=0
    bcNum=0
    bdNum=0
    cdNum=0
    abcNum=0
    abdNum=0
    acdNum=0
    bcdNum=0
    abcdNum=0
    consisNum=0
    nonConsisNum=0
    allNonConsisNum=0
    allCpgNum=0             ## The total number of CpG sites in all software (duplicated)
    
    reportFileName=outPut+"_report.txt"
    consisCpgFileName=outPut+"_fourDetec_consisCpg.bed"
    nonConsisCpgFileName=outPut+"_fourDetec_nonConsisCpg.bed"
    consisCpgMethFileName=outPut+"_fourDetect_consisCpgMeth.bed"
    nonConsisCpgMethFileName=outPut+"_fourDetec_nonConsisCpgMeth.bed"
    diffNonConsisCpgFileName=outPut+"_fourDiff_all.bed"
    
    reportFile=open(reportFileName,"w")
    consisCpgFile=open(consisCpgFileName,"w")
    nonConsisCpgFile=open(nonConsisCpgFileName,"w")
    consisCpgMethFile=open(consisCpgMethFileName,"w")
    nonConsisCpgMethFile=open(nonConsisCpgMethFileName,"w")
    diffNonConsisCpgFile=open(diffNonConsisCpgFileName,"w")
    
    aMapperName=os.path.basename(aFile).split('_')[0]
    bMapperName=os.path.basename(bFile).split('_')[0]
    cMapperName=os.path.basename(cFile).split('_')[0]
    dMapperName=os.path.basename(dFile).split('_')[0]
    
    aNonConsisCpgFileName=outPut+"_"+aMapperName+"_diffDetecCpg.bed"
    bNonConsisCpgFileName=outPut+"_"+bMapperName+"_diffDetecCpg.bed"
    cNonConsisCpgFileName=outPut+"_"+cMapperName+"_diffDetecCpg.bed"
    dNonConsisCpgFileName=outPut+"_"+dMapperName+"_diffDetecCpg.bed"
    
    aNonConsisCpgFile=open(aNonConsisCpgFileName,"w")
    bNonConsisCpgFile=open(bNonConsisCpgFileName,"w")
    cNonConsisCpgFile=open(cNonConsisCpgFileName,"w")
    dNonConsisCpgFile=open(dNonConsisCpgFileName,"w")
    
    for fileid in fileList:
        aDict={}
        bDict={}
        cDict={}
        dDict={}
        allDict={}
        afilename=''.join((atempfilename,"/",fileid,".txt"))
        if os.path.exists(afilename):
             with open(afilename,'r') as f:
                for line in f:
                    aNum += 1
                    line=line.split()
                    aDict[line[0]+"_"+line[1]+"_"+line[2]]=line
                    allDict[line[0]+"_"+line[1]+"_"+line[2]]=line
        bfilename=''.join((btempfilename,"/",fileid,".txt"))
        if os.path.exists(bfilename):
             with open(bfilename,'r') as f:
                for line in f:
                    bNum += 1
                    line=line.split()
                    bDict[line[0]+"_"+line[1]+"_"+line[2]]=line
                    allDict[line[0]+"_"+line[1]+"_"+line[2]]=line
        cfilename=''.join((ctempfilename,"/",fileid,".txt"))
        if os.path.exists(cfilename):
             with open(cfilename,'r') as f:
                for line in f:
                    cNum += 1
                    line=line.split()
                    cDict[line[0]+"_"+line[1]+"_"+line[2]]=line
                    allDict[line[0]+"_"+line[1]+"_"+line[2]]=line
        dfilename=''.join((dtempfilename,"/",fileid,".txt"))
        if os.path.exists(dfilename):
             with open(dfilename,'r') as f:
                for line in f:
                    dNum += 1
                    line=line.split()
                    dDict[line[0]+"_"+line[1]+"_"+line[2]]=line
                    allDict[line[0]+"_"+line[1]+"_"+line[2]]=line
        def twoDict(aDict,bDict):
            conNum=0
            aSet=set(aDict.keys())
            bSet=set(bDict.keys())
            inter=aSet & bSet
            for i in inter:
                aMeth=int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5]))
                bMeth=int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5]))
                diffMeth=abs(aMeth-bMeth)
                if diffMeth <= 0.05:
                    conNum += 1
            return(conNum)
        def threeDict(aDict,bDict,cDict):
            conNum=0
            aSet=set(aDict.keys())
            bSet=set(bDict.keys())
            cSet=set(cDict.keys())
            inter=aSet & bSet & cSet
            for i in inter:
                aMeth=int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5]))
                bMeth=int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5]))
                cMeth=int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5]))
                diffMeth=max(aMeth,bMeth,cMeth)-min(aMeth,bMeth,cMeth)
                if diffMeth <= 0.05:
                    conNum += 1
            return(conNum)
        def fourDict(aDict,bDict,cDict,dDict):
            conNum=0
            aSet=set(aDict.keys())
            bSet=set(bDict.keys())
            cSet=set(cDict.keys())
            dSet=set(dDict.keys())
            inter=aSet & bSet & cSet & dSet
            for i in inter:
                aMeth=int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5]))
                bMeth=int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5]))
                cMeth=int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5]))
                dMeth=int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5]))
                diffMeth=max(aMeth,bMeth,cMeth,dMeth)-min(aMeth,bMeth,cMeth,dMeth)
                if diffMeth <= 0.05:
                    conNum += 1
            return(conNum)
        
        aSet=set(aDict.keys())
        bSet=set(bDict.keys())
        cSet=set(cDict.keys())
        dSet=set(dDict.keys())
        interSet=aSet & bSet & cSet & dSet
        allCpgSet=aSet | bSet | cSet | dSet
        ## Intersection of CpG site results of various software
        abNum += twoDict(aDict, bDict)
        acNum += twoDict(aDict, cDict)
        adNum += twoDict(aDict, dDict)
        bcNum += twoDict(bDict, cDict)
        bdNum += twoDict(bDict, dDict)
        cdNum += twoDict(cDict, dDict)
        abcNum += threeDict(aDict,bDict,cDict)
        abdNum += threeDict(aDict,bDict,dDict)
        acdNum += threeDict(aDict,cDict,dDict)
        bcdNum += threeDict(bDict,cDict,dDict)
        abcdNum += fourDict(aDict,bDict,cDict,dDict)
        allCpgNum += len(allCpgSet)
        
        ## Distinguish between consistent CpG sites and non-identical CpG sites      
        a_interCpgSet = aSet - interSet
        b_interCpgSet = bSet - interSet
        c_interCpgSet = cSet - interSet
        d_interCpgSet = dSet - interSet
        all_interCpgSet = allCpgSet - interSet

        for i in a_interCpgSet:
            aNonConsisCpgFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+str(int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5])))+"\n")
        for i in b_interCpgSet:
            bNonConsisCpgFile.write(bDict[i][0]+"\t"+bDict[i][1]+"\t"+bDict[i][2]+"\t"+str(int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5])))+"\n")  
        for i in c_interCpgSet:
            cNonConsisCpgFile.write(cDict[i][0]+"\t"+cDict[i][1]+"\t"+cDict[i][2]+"\t"+str(int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5])))+"\n")  
        for i in d_interCpgSet:
            dNonConsisCpgFile.write(dDict[i][0]+"\t"+dDict[i][1]+"\t"+dDict[i][2]+"\t"+str(int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5])))+"\n")  
        for i in all_interCpgSet:
           diffNonConsisCpgFile.write(allDict[i][0]+"\t"+allDict[i][1]+"\t"+allDict[i][2]+"\t"+"0"+"\n")
            
        for i in interSet:
            aMethLevel=int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5]))
            bMethLevel=int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5]))
            cMethLevel=int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5]))
            dMethLevel=int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5]))
            diffMethLevel=max(aMethLevel,bMethLevel,cMethLevel,dMethLevel)-min(aMethLevel,bMethLevel,cMethLevel,dMethLevel)
            if diffMethLevel <= 0.05:
                consisNum += 1
                consisCpgFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+"0"+"\n")
                consisCpgMethFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+str(int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5])))+"\n") 
                consisCpgMethFile.write(bDict[i][0]+"\t"+bDict[i][1]+"\t"+bDict[i][2]+"\t"+str(int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5])))+"\n") 
                consisCpgMethFile.write(cDict[i][0]+"\t"+cDict[i][1]+"\t"+cDict[i][2]+"\t"+str(int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5])))+"\n") 
                consisCpgMethFile.write(dDict[i][0]+"\t"+dDict[i][1]+"\t"+dDict[i][2]+"\t"+str(int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5])))+"\n") 
            else:
                nonConsisNum += 1
                nonConsisCpgFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+"0"+"\n")        
                nonConsisCpgMethFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+str(int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5])))+"\n") 
                nonConsisCpgMethFile.write(bDict[i][0]+"\t"+bDict[i][1]+"\t"+bDict[i][2]+"\t"+str(int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5])))+"\n") 
                nonConsisCpgMethFile.write(cDict[i][0]+"\t"+cDict[i][1]+"\t"+cDict[i][2]+"\t"+str(int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5])))+"\n") 
                nonConsisCpgMethFile.write(dDict[i][0]+"\t"+dDict[i][1]+"\t"+dDict[i][2]+"\t"+str(int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5])))+"\n") 
    
    allNonConsisNum = allCpgNum - consisNum         
    reportFile.write("aMapperName"+"\t"+aMapperName+"\n"+
                     "bMapperName"+"\t"+bMapperName+"\n"+
                     "cMapperName"+"\t"+cMapperName+"\n"+
                     "dMapperName"+"\t"+dMapperName+"\n"+
                     "aNum"+"\t"+str(aNum)+"\n"+
                     "bNum"+"\t"+str(bNum)+"\n"+
                     "cNum"+"\t"+str(cNum)+"\n"+
                     "dNum"+"\t"+str(dNum)+"\n"+
                     "abNum"+"\t"+str(abNum)+"\n"+
                     "acNum"+"\t"+str(acNum)+"\n"+
                     "adNum"+"\t"+str(adNum)+"\n"+
                     "bcNum"+"\t"+str(bcNum)+"\n"+
                     "bdNum"+"\t"+str(bdNum)+"\n"+
                     "cdNum"+"\t"+str(cdNum)+"\n"+
                     "abcNum"+"\t"+str(abcNum)+"\n"+
                     "abdNum"+"\t"+str(abdNum)+"\n"+
                     "acdNum"+"\t"+str(acdNum)+"\n"+
                     "bcdNum"+"\t"+str(bcdNum)+"\n"+
                     "abcdNum"+"\t"+str(abcdNum)+"\n"+
                     "consisNum"+"\t"+str(consisNum)+"\n"+
                     "nonConsisNum"+"\t"+str(nonConsisNum)+"\n"+
                     "allNonConsisNum"+"\t"+str(allNonConsisNum)+"\n"+
                     "allCpgNum"+"\t"+str(allCpgNum)+"\n")
    os.system("rm -r %s" % atempfilename)
    os.system("rm -r %s" % btempfilename)
    os.system("rm -r %s" % ctempfilename)
    os.system("rm -r %s" % dtempfilename)

resOut(fileList,atempfilename,btempfilename,ctempfilename,dtempfilename,outPut,aFile,bFile,cFile,dFile)

           
### Regional annotation and result statistics of Cpg sites
aMapperName=os.path.basename(aFile).split('_')[0]
bMapperName=os.path.basename(bFile).split('_')[0]
cMapperName=os.path.basename(cFile).split('_')[0]
dMapperName=os.path.basename(dFile).split('_')[0]

consisCpgFileName=outPut+"_fourDetec_consisCpg.bed"
nonConsisCpgFileName=outPut+"_fourDetec_nonConsisCpg.bed"
diffNonConsisCpgFileName=outPut+"_fourDiff_all.bed"
aNonConsisCpgFileName=outPut+"_"+aMapperName+"_diffDetecCpg.bed"
bNonConsisCpgFileName=outPut+"_"+bMapperName+"_diffDetecCpg.bed"
cNonConsisCpgFileName=outPut+"_"+cMapperName+"_diffDetecCpg.bed"
dNonConsisCpgFileName=outPut+"_"+dMapperName+"_diffDetecCpg.bed"


def resSta(inFileName,dictName,lineNum):
    dictName={}
    inFile=open(inFileName,"r")
    for line in inFile:
        resClass=line.split()[lineNum]
        if resClass in dictName:
            dictName[resClass] += 1
        else:
            dictName[resClass] = 1
    return(dictName) 

def GgiAndGeneAnnoResSta(consisCpgAnnoResBed,nonConsisCpgAnnoResBed,diffNonConsisCpgAnnoResBed,aNonConsisCpgAnnoResBed,bNonConsisCpgAnnoResBed,cNonConsisCpgAnnoResBed,dNonConsisCpgAnnoResBed,outFileName,lineNum):
    outFile=open(outFileName,"w")
    
    consisCpgDict={}
    nonConsisCpgDict={}
    diffNonConsisCpgDict={}
    aNonConsisCpgDict={}
    bNonConsisCpgDict={}
    cNonConsisCpgDict={}
    dNonConsisCpgDict={}
    
    consisCpgDict=resSta(consisCpgAnnoResBed,consisCpgDict,lineNum)
    nonConsisCpgDict=resSta(nonConsisCpgAnnoResBed,nonConsisCpgDict,lineNum)
    diffNonConsisCpgDict=resSta(diffNonConsisCpgAnnoResBed,diffNonConsisCpgDict,lineNum)
    aNonConsisCpgDict=resSta(aNonConsisCpgAnnoResBed,aNonConsisCpgDict,lineNum)
    bNonConsisCpgDict=resSta(bNonConsisCpgAnnoResBed,bNonConsisCpgDict,lineNum)
    cNonConsisCpgDict=resSta(cNonConsisCpgAnnoResBed,cNonConsisCpgDict,lineNum)
    dNonConsisCpgDict=resSta(dNonConsisCpgAnnoResBed,dNonConsisCpgDict,lineNum)
    
    allClassSet=set(consisCpgDict.keys()) | set(nonConsisCpgDict.keys()) | set(diffNonConsisCpgDict.keys()) | set(aNonConsisCpgDict.keys()) | set(bNonConsisCpgDict.keys()) | set(cNonConsisCpgDict.keys()) | set(dNonConsisCpgDict.keys())
    
    consisCpgAnnoResBedFileName=os.path.basename(consisCpgAnnoResBed).split(".")[0]
    nonConsisCpgAnnoResBedFileName=os.path.basename(nonConsisCpgAnnoResBed).split(".")[0]
    diffNonConsisCpgAnnoResBedFileName=os.path.basename(diffNonConsisCpgAnnoResBed).split(".")[0]
    aNonConsisCpgAnnoResBedFileName=os.path.basename(aNonConsisCpgAnnoResBed).split(".")[0]
    bNonConsisCpgAnnoResBedFileName=os.path.basename(bNonConsisCpgAnnoResBed).split(".")[0]
    cNonConsisCpgAnnoResBedFileName=os.path.basename(cNonConsisCpgAnnoResBed).split(".")[0]
    dNonConsisCpgAnnoResBedFileName=os.path.basename(dNonConsisCpgAnnoResBed).split(".")[0]
    
    outFile.write("resClass"+"\t"+consisCpgAnnoResBedFileName+"\t"+nonConsisCpgAnnoResBedFileName+"\t"+\
                  diffNonConsisCpgAnnoResBedFileName+"\t"+aNonConsisCpgAnnoResBedFileName+"\t"+\
                     bNonConsisCpgAnnoResBedFileName+"\t"+cNonConsisCpgAnnoResBedFileName+"\t"+\
                         dNonConsisCpgAnnoResBedFileName+"\n")
    for resClass in allClassSet:
        if resClass in consisCpgDict:
            consisCpg = consisCpgDict[resClass]
        else:
            consisCpg = 0
        if resClass in nonConsisCpgDict:
            nonConsisCpg = nonConsisCpgDict[resClass]
        else:
            nonConsisCpg = 0
        if resClass in diffNonConsisCpgDict:
            diffNonConsisCpg = diffNonConsisCpgDict[resClass]
        else:
            diffNonConsisCpg = 0
        if resClass in aNonConsisCpgDict:
            aNonConsisCpg = aNonConsisCpgDict[resClass]
        else:
            aNonConsisCpg = 0
        if resClass in bNonConsisCpgDict:
            bNonConsisCpg = bNonConsisCpgDict[resClass]
        else:
            bNonConsisCpg = 0
        if resClass in cNonConsisCpgDict:
            cNonConsisCpg = cNonConsisCpgDict[resClass]
        else:
            cNonConsisCpg = 0
        if resClass in dNonConsisCpgDict:
            dNonConsisCpg = dNonConsisCpgDict[resClass]
        else:
            dNonConsisCpg = 0 
        outFile.write(resClass+"\t"+str(consisCpg)+"\t"+str(nonConsisCpg)+"\t"+str(diffNonConsisCpg)+"\t"+
                      str(aNonConsisCpg)+"\t"+str(bNonConsisCpg)+"\t"+str(cNonConsisCpg)+"\t"+str(dNonConsisCpg)+"\n")


## Cpg realated region annotation results of cpg locus
consisCpgAnnoResBedCGI=outPut+"_fourDetec_consisCpgAnnoResCGI.bed"
nonConsisCpgAnnoResBedCGI=outPut+"_fourDetec_nonConsisCpgAnnoResCGI.bed"
diffNonConsisCpgAnnoResBedCGI=outPut+"_fourDiff_allAnnoResCGI.bed"
aNonConsisCpgAnnoResBedCGI=outPut+"_"+aMapperName+"_diffDetecCpgAnnoResCGI.bed"
bNonConsisCpgAnnoResBedCGI=outPut+"_"+bMapperName+"_diffDetecCpgAnnoResCGI.bed"
cNonConsisCpgAnnoResBedCGI=outPut+"_"+cMapperName+"_diffDetecCpgAnnoResCGI.bed"
dNonConsisCpgAnnoResBedCGI=outPut+"_"+dMapperName+"_diffDetecCpgAnnoResCGI.bed"

os.system('bedtools intersect -a %s -b %s -wao > %s' %(consisCpgFileName,cgiAnnoFile,consisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(nonConsisCpgFileName,cgiAnnoFile,nonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(diffNonConsisCpgFileName,cgiAnnoFile,diffNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(aNonConsisCpgFileName,cgiAnnoFile,aNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(bNonConsisCpgFileName,cgiAnnoFile,bNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(cNonConsisCpgFileName,cgiAnnoFile,cNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(dNonConsisCpgFileName,cgiAnnoFile,dNonConsisCpgAnnoResBedCGI))

resStaOfCgiAnno=outPut+"resStaOfCgiAnno.txt"
GgiAndGeneAnnoResSta(consisCpgAnnoResBedCGI,nonConsisCpgAnnoResBedCGI,diffNonConsisCpgAnnoResBedCGI,aNonConsisCpgAnnoResBedCGI,bNonConsisCpgAnnoResBedCGI,cNonConsisCpgAnnoResBedCGI,dNonConsisCpgAnnoResBedCGI,resStaOfCgiAnno,7)

## repeat region annotation results of cpg locus
consisCpgAnnoResBedRep=outPut+"_fourDetec_consisCpgAnnoResRep.bed"
nonConsisCpgAnnoResBedRep=outPut+"_fourDetec_nonConsisCpgAnnoResRep.bed"
diffNonConsisCpgAnnoResBedRep=outPut+"_fourDiff_allAnnoResRep.bed"
aNonConsisCpgAnnoResBedRep=outPut+"_"+aMapperName+"_diffDetecCpgAnnoResRep.bed"
bNonConsisCpgAnnoResBedRep=outPut+"_"+bMapperName+"_diffDetecCpgAnnoResRep.bed"
cNonConsisCpgAnnoResBedRep=outPut+"_"+cMapperName+"_diffDetecCpgAnnoResRep.bed"
dNonConsisCpgAnnoResBedRep=outPut+"_"+dMapperName+"_diffDetecCpgAnnoResRep.bed"

os.system('bedtools intersect -a %s -b %s -wao > %s' %(consisCpgFileName,repAnnoFile,consisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(nonConsisCpgFileName,repAnnoFile,nonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(diffNonConsisCpgFileName,repAnnoFile,diffNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(aNonConsisCpgFileName,repAnnoFile,aNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(bNonConsisCpgFileName,repAnnoFile,bNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(cNonConsisCpgFileName,repAnnoFile,cNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(dNonConsisCpgFileName,repAnnoFile,dNonConsisCpgAnnoResBedRep))

resStaOfRepAnno=outPut+"resStaOfRepAnno.txt"
GgiAndGeneAnnoResSta(consisCpgAnnoResBedRep,nonConsisCpgAnnoResBedRep,diffNonConsisCpgAnnoResBedRep,aNonConsisCpgAnnoResBedRep,bNonConsisCpgAnnoResBedRep,cNonConsisCpgAnnoResBedRep,dNonConsisCpgAnnoResBedRep,resStaOfRepAnno,9)


## Distribution of methylation levels of consistent CpG sites and non-identical CpG sites
consisCpgMethFileName=outPut+"_fourDetect_consisCpgMeth.bed"
nonConsisCpgMethFileName=outPut+"_fourDetec_nonConsisCpgMeth.bed"
resStaOfRepAnno=outPut+"cpgMehtLevelDistri.txt"

def methLevelSta(consisCpgMethFileName,nonConsisCpgMethFileName,resStaOfRepAnnoFileName):
    consisCpgMethFile=open(consisCpgMethFileName,"r")
    nonConsisCpgMethFile=open(nonConsisCpgMethFileName,"r")
    resStaOfRepAnnoFile=open(resStaOfRepAnnoFileName,"w")
    
    consisHighMethNum = 0
    consisMediumMethNum = 0
    consisLowMehtNum = 0
    
    nonConsisHighMethNum = 0
    nonConsisMediumMethNum = 0
    nonConsisLowMehtNum = 0    
    
    for line in consisCpgMethFile:
        methLevel=float(line.split()[3])
        if methLevel >= 2/3:
            consisHighMethNum += 1
        elif methLevel < 1/3:
            consisLowMehtNum += 1
        else:
            consisMediumMethNum += 1
    
    for line in nonConsisCpgMethFile:
        methLevel = float(line.split()[3])
        if methLevel >= 2/3:
            nonConsisHighMethNum += 1
        elif methLevel < 1/3:
            nonConsisLowMehtNum += 1
        else:
            nonConsisMediumMethNum += 1
            
    resStaOfRepAnnoFile.write("class"+"\t"+"consis"+"\t"+"nonConsis"+"\n")
    resStaOfRepAnnoFile.write("highMeth"+"\t"+str(consisHighMethNum)+"\t"+str(nonConsisHighMethNum)+"\n")
    resStaOfRepAnnoFile.write("mediumMeth"+"\t"+str(consisMediumMethNum)+"\t"+str(nonConsisMediumMethNum)+"\n") 
    resStaOfRepAnnoFile.write("lowMeth"+"\t"+str(consisLowMehtNum)+"\t"+str(nonConsisLowMehtNum)+"\n") 

methLevelSta(consisCpgMethFileName,nonConsisCpgMethFileName,resStaOfRepAnno)    
