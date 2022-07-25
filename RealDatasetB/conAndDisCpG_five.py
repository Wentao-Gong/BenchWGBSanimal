# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:04:57 2022

一致性CpG位点：是指同一CpG位点被两个或多个软件检测（即测序深度>5或者>10）,且其甲基化水平在两个或多个软件中变化范围小于0.05.
非一致性CpG位点：是指一致性CpG位点外的CpG位点

本脚本运行的主要流程：
1、将大文件切割成小文件
2、筛选一致性的CpG位点和统计结果
3、基因注释及结果统计
4、CGI注释及结果统计
5、重复序列区域注释及结果统计
6、甲基化水平分布的统计

@author: 滔
"""
from __future__ import division
import os
import argparse
import threading
parser = argparse.ArgumentParser(description='Extract the result of the difference set of two files:set(a)-set(b)')
parser.add_argument('-a', '--afile', help='a file', required=True)
parser.add_argument('-b', '--bfile', help='b file', required=True)
parser.add_argument('-c', '--cfile', help='c file', required=True)
parser.add_argument('-d', '--dfile', help='d file', required=True)
parser.add_argument('-e', '--efile', help='e file', required=True)
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
eFile = args.efile
tempPath = args.temppath
outPut =args.output
cpgDepth=args.cpgdepth
repAnnoFile = args.repanno
cgiAnnoFile = args.cgianno

## multiple thread
class MyThread(threading.Thread):
    def __init__(self,func,args=()):
        super(MyThread, self).__init__()
        self.func=func
        self.args=args
    
    def run(self):
        self.result=self.func(*self.args)
        
    def get_result(self):
        threading.Thread.join(self)
        try:
            return self.result
        except Exception:
            return None

## According to the chromosome number, cut large files into small files

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
btempfilename=''.join((tempPath,"/btemp"))
if not os.path.exists(btempfilename):
	os.mkdir(btempfilename)
ctempfilename=''.join((tempPath,"/ctemp"))
if not os.path.exists(ctempfilename):
	os.mkdir(ctempfilename)
dtempfilename=''.join((tempPath,"/dtemp"))
if not os.path.exists(dtempfilename):
	os.mkdir(dtempfilename)
etempfilename=''.join((tempPath,"/etemp"))
if not os.path.exists(etempfilename):
	os.mkdir(etempfilename)

threads=[]
AfileList=[]
BfileList=[]
CfileList=[]
DfileList=[]
EfileList=[]
AcutFilesByChrom=MyThread(cutFilesByChrom,args=(aFile,AfileList,atempfilename,cpgDepth,))
BcutFilesByChrom=MyThread(cutFilesByChrom,args=(bFile,BfileList,btempfilename,cpgDepth,))
CcutFilesByChrom=MyThread(cutFilesByChrom,args=(cFile,CfileList,ctempfilename,cpgDepth,))
DcutFilesByChrom=MyThread(cutFilesByChrom,args=(dFile,DfileList,dtempfilename,cpgDepth,))
EcutFilesByChrom=MyThread(cutFilesByChrom,args=(eFile,EfileList,etempfilename,cpgDepth,))
threads.append(AcutFilesByChrom)
threads.append(BcutFilesByChrom)
threads.append(CcutFilesByChrom)
threads.append(DcutFilesByChrom)
threads.append(EcutFilesByChrom)

for thr in threads:
    thr.start()
for thr in threads:
    if thr.is_alive():
        thr.join()
AfileList=AcutFilesByChrom.get_result()
BfileList=BcutFilesByChrom.get_result()
CfileList=CcutFilesByChrom.get_result()
DfileList=DcutFilesByChrom.get_result()
EfileList=EcutFilesByChrom.get_result()

fileList=list(set(AfileList+BfileList+CfileList+DfileList+EfileList))


## Screening consistent CpG sites and statistical results
def resOut(fileList,atempfilename,btempfilename,ctempfilename,dtempfilename,etempfilename,outPut,aFile,bFile,cFile,dFile,eFile):
    aNum=0
    bNum=0
    cNum=0
    dNum=0
    eNum=0
    abNum=0
    acNum=0
    adNum=0
    aeNum=0
    bcNum=0
    bdNum=0
    beNum=0
    cdNum=0
    ceNum=0
    deNum=0
    abcNum=0
    abdNum=0
    abeNum=0
    acdNum=0
    aceNum=0
    adeNum=0
    bcdNum=0
    bceNum=0
    bdeNum=0
    cdeNum=0
    abcdNum=0
    abceNum=0
    abdeNum=0
    acdeNum=0
    bcdeNum=0
    abcdeNum=0
    consisNum=0
    nonConsisNum=0
    allNonConsisNum=0
    allCpgNum=0             ## The total number of CpG sites in all software (duplicated)
    
    reportFileName=outPut+"_report.txt"
    consisCpgFileName=outPut+"_fiveDetec_consisCpg.bed"
    nonConsisCpgFileName=outPut+"_fiveDetec_nonConsisCpg.bed"
    consisCpgMethFileName=outPut+"_fiveDetect_consisCpgMeth.bed"
    nonConsisCpgMethFileName=outPut+"_fiveDetec_nonConsisCpgMeth.bed"
    diffNonConsisCpgFileName=outPut+"_fiveDiff_all.bed"
    
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
    eMapperName=os.path.basename(eFile).split('_')[0]
    
    aNonConsisCpgFileName=outPut+"_"+aMapperName+"_diffDetecCpg.bed"
    bNonConsisCpgFileName=outPut+"_"+bMapperName+"_diffDetecCpg.bed"
    cNonConsisCpgFileName=outPut+"_"+cMapperName+"_diffDetecCpg.bed"
    dNonConsisCpgFileName=outPut+"_"+dMapperName+"_diffDetecCpg.bed"
    eNonConsisCpgFileName=outPut+"_"+eMapperName+"_diffDetecCpg.bed"
    
    aNonConsisCpgFile=open(aNonConsisCpgFileName,"w")
    bNonConsisCpgFile=open(bNonConsisCpgFileName,"w")
    cNonConsisCpgFile=open(cNonConsisCpgFileName,"w")
    dNonConsisCpgFile=open(dNonConsisCpgFileName,"w")
    eNonConsisCpgFile=open(eNonConsisCpgFileName,"w")
    
    for fileid in fileList:
        aDict={}
        bDict={}
        cDict={}
        dDict={}
        eDict={}
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
        efilename=''.join((etempfilename,"/",fileid,".txt"))
        if os.path.exists(efilename):
             with open(efilename,'r') as f:
                for line in f:
                    eNum += 1
                    line=line.split()
                    eDict[line[0]+"_"+line[1]+"_"+line[2]]=line
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
        def fiveDict(aDict,bDict,cDict,dDict,eDict):
            conNum=0
            aSet=set(aDict.keys())
            bSet=set(bDict.keys())
            cSet=set(cDict.keys())
            dSet=set(dDict.keys())
            eSet=set(eDict.keys())
            inter=aSet & bSet & cSet & dSet & eSet
            for i in inter:
                aMeth=int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5]))
                bMeth=int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5]))
                cMeth=int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5]))
                dMeth=int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5]))
                eMeth=int(eDict[i][4])/(int(eDict[i][4])+int(eDict[i][5]))
                diffMeth=max(aMeth,bMeth,cMeth,dMeth,eMeth)-min(aMeth,bMeth,cMeth,dMeth,eMeth)
                if diffMeth <= 0.05:
                    conNum += 1
            return(conNum)
        
        
        aSet=set(aDict.keys())
        bSet=set(bDict.keys())
        cSet=set(cDict.keys())
        dSet=set(dDict.keys())
        eSet=set(eDict.keys())
        interSet=aSet & bSet & cSet & dSet & eSet
        allCpgSet=aSet | bSet | cSet | dSet | eSet
        ## Intersection of CpG site results of various software
        abNum += twoDict(aDict, bDict)
        acNum += twoDict(aDict, cDict)
        adNum += twoDict(aDict, dDict)
        aeNum += twoDict(aDict, eDict)
        bcNum += twoDict(bDict, cDict)
        bdNum += twoDict(bDict, dDict)
        beNum += twoDict(bDict, eDict)
        cdNum += twoDict(cDict, dDict)
        ceNum += twoDict(cDict, eDict)
        deNum += twoDict(dDict, eDict)
        abcNum += threeDict(aDict,bDict,cDict)
        abdNum += threeDict(aDict,bDict,dDict)
        abeNum += threeDict(aDict,bDict,eDict)        
        acdNum += threeDict(aDict,cDict,dDict)
        aceNum += threeDict(aDict,cDict,eDict)
        adeNum += threeDict(aDict,dDict,eDict)
        bcdNum += threeDict(bDict,cDict,dDict)
        bceNum += threeDict(bDict,cDict,eDict)
        bdeNum += threeDict(bDict,dDict,eDict)
        cdeNum += threeDict(cDict,dDict,eDict)
        abcdNum += fourDict(aDict,bDict,cDict,dDict)
        abceNum += fourDict(aDict,bDict,cDict,eDict)
        abdeNum += fourDict(aDict,bDict,dDict,eDict)
        acdeNum += fourDict(aDict,cDict,dDict,eDict)
        bcdeNum += fourDict(bDict,cDict,dDict,eDict)
        abcdeNum += fiveDict(aDict,bDict,cDict,dDict,eDict)
        allCpgNum += len(allCpgSet)
        
        ## Distinguish between consistent CpG sites and non-identical CpG sites      
        a_interCpgSet = aSet - interSet
        b_interCpgSet = bSet - interSet
        c_interCpgSet = cSet - interSet
        d_interCpgSet = dSet - interSet
        e_interCpgSet = eSet - interSet
        all_interCpgSet = allCpgSet - interSet

        for i in a_interCpgSet:
            aNonConsisCpgFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+str(int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5])))+"\n")
        for i in b_interCpgSet:
            bNonConsisCpgFile.write(bDict[i][0]+"\t"+bDict[i][1]+"\t"+bDict[i][2]+"\t"+str(int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5])))+"\n")  
        for i in c_interCpgSet:
            cNonConsisCpgFile.write(cDict[i][0]+"\t"+cDict[i][1]+"\t"+cDict[i][2]+"\t"+str(int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5])))+"\n")  
        for i in d_interCpgSet:
            dNonConsisCpgFile.write(dDict[i][0]+"\t"+dDict[i][1]+"\t"+dDict[i][2]+"\t"+str(int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5])))+"\n")  
        for i in e_interCpgSet:
            eNonConsisCpgFile.write(eDict[i][0]+"\t"+eDict[i][1]+"\t"+eDict[i][2]+"\t"+str(int(eDict[i][4])/(int(eDict[i][4])+int(eDict[i][5])))+"\n")  
        for i in all_interCpgSet:
           diffNonConsisCpgFile.write(allDict[i][0]+"\t"+allDict[i][1]+"\t"+allDict[i][2]+"\t"+"0"+"\n")
            
        for i in interSet:
            aMethLevel=int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5]))
            bMethLevel=int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5]))
            cMethLevel=int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5]))
            dMethLevel=int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5]))
            eMethLevel=int(eDict[i][4])/(int(eDict[i][4])+int(eDict[i][5]))
            diffMethLevel=max(aMethLevel,bMethLevel,cMethLevel,dMethLevel,eMethLevel)-min(aMethLevel,bMethLevel,cMethLevel,dMethLevel,eMethLevel)
            if diffMethLevel <= 0.05:
                consisNum += 1
                consisCpgFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+"0"+"\n")
                consisCpgMethFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+str(int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5])))+"\n") 
                consisCpgMethFile.write(bDict[i][0]+"\t"+bDict[i][1]+"\t"+bDict[i][2]+"\t"+str(int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5])))+"\n") 
                consisCpgMethFile.write(cDict[i][0]+"\t"+cDict[i][1]+"\t"+cDict[i][2]+"\t"+str(int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5])))+"\n") 
                consisCpgMethFile.write(dDict[i][0]+"\t"+dDict[i][1]+"\t"+dDict[i][2]+"\t"+str(int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5])))+"\n")
                consisCpgMethFile.write(eDict[i][0]+"\t"+eDict[i][1]+"\t"+eDict[i][2]+"\t"+str(int(eDict[i][4])/(int(eDict[i][4])+int(eDict[i][5])))+"\n") 
            else:
                nonConsisNum += 1
                nonConsisCpgFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+"0"+"\n")        
                nonConsisCpgMethFile.write(aDict[i][0]+"\t"+aDict[i][1]+"\t"+aDict[i][2]+"\t"+str(int(aDict[i][4])/(int(aDict[i][4])+int(aDict[i][5])))+"\n") 
                nonConsisCpgMethFile.write(bDict[i][0]+"\t"+bDict[i][1]+"\t"+bDict[i][2]+"\t"+str(int(bDict[i][4])/(int(bDict[i][4])+int(bDict[i][5])))+"\n") 
                nonConsisCpgMethFile.write(cDict[i][0]+"\t"+cDict[i][1]+"\t"+cDict[i][2]+"\t"+str(int(cDict[i][4])/(int(cDict[i][4])+int(cDict[i][5])))+"\n") 
                nonConsisCpgMethFile.write(dDict[i][0]+"\t"+dDict[i][1]+"\t"+dDict[i][2]+"\t"+str(int(dDict[i][4])/(int(dDict[i][4])+int(dDict[i][5])))+"\n") 
                nonConsisCpgMethFile.write(eDict[i][0]+"\t"+eDict[i][1]+"\t"+eDict[i][2]+"\t"+str(int(eDict[i][4])/(int(eDict[i][4])+int(eDict[i][5])))+"\n") 
    
    allNonConsisNum = allCpgNum - consisNum         
    reportFile.write("aMapperName"+"\t"+aMapperName+"\n"+
                     "bMapperName"+"\t"+bMapperName+"\n"+
                     "cMapperName"+"\t"+cMapperName+"\n"+
                     "dMapperName"+"\t"+dMapperName+"\n"+
                     "eMapperName"+"\t"+eMapperName+"\n"+
                     "aNum"+"\t"+str(aNum)+"\n"+
                     "bNum"+"\t"+str(bNum)+"\n"+
                     "cNum"+"\t"+str(cNum)+"\n"+
                     "dNum"+"\t"+str(dNum)+"\n"+
                     "eNum"+"\t"+str(eNum)+"\n"+
                     "abNum"+"\t"+str(abNum)+"\n"+
                     "acNum"+"\t"+str(acNum)+"\n"+
                     "adNum"+"\t"+str(adNum)+"\n"+
                     "aeNum"+"\t"+str(aeNum)+"\n"+
                     "bcNum"+"\t"+str(bcNum)+"\n"+
                     "bdNum"+"\t"+str(bdNum)+"\n"+
                     "beNum"+"\t"+str(beNum)+"\n"+
                     "cdNum"+"\t"+str(cdNum)+"\n"+
                     "ceNum"+"\t"+str(ceNum)+"\n"+
                     "deNum"+"\t"+str(deNum)+"\n"+
                     "abcNum"+"\t"+str(abcNum)+"\n"+
                     "abdNum"+"\t"+str(abdNum)+"\n"+
                     "abeNum"+"\t"+str(abeNum)+"\n"+
                     "acdNum"+"\t"+str(acdNum)+"\n"+
                     "aceNum"+"\t"+str(aceNum)+"\n"+
                     "adeNum"+"\t"+str(adeNum)+"\n"+                     
                     "bcdNum"+"\t"+str(bcdNum)+"\n"+
                     "bceNum"+"\t"+str(bceNum)+"\n"+
                     "bdeNum"+"\t"+str(bdeNum)+"\n"+
                     "cdeNum"+"\t"+str(cdeNum)+"\n"+
                     "abcdNum"+"\t"+str(abcdNum)+"\n"+
                     "abceNum"+"\t"+str(abceNum)+"\n"+
                     "abdeNum"+"\t"+str(abdeNum)+"\n"+
                     "acdeNum"+"\t"+str(acdeNum)+"\n"+
                     "bcdeNum"+"\t"+str(bcdeNum)+"\n"+
                     "abcdeNum"+"\t"+str(abcdeNum)+"\n"+
                     "consisNum"+"\t"+str(consisNum)+"\n"+
                     "nonConsisNum"+"\t"+str(nonConsisNum)+"\n"+
                     "allNonConsisNum"+"\t"+str(allNonConsisNum)+"\n"+
                     "allCpgNum"+"\t"+str(allCpgNum)+"\n")
    os.system("rm -r %s" % atempfilename)
    os.system("rm -r %s" % btempfilename)
    os.system("rm -r %s" % ctempfilename)
    os.system("rm -r %s" % dtempfilename)
    os.system("rm -r %s" % etempfilename)

resOut(fileList,atempfilename,btempfilename,ctempfilename,dtempfilename,etempfilename,outPut,aFile,bFile,cFile,dFile,eFile)

           
### Regional annotation and result statistics of Cpg sites
aMapperName=os.path.basename(aFile).split('_')[0]
bMapperName=os.path.basename(bFile).split('_')[0]
cMapperName=os.path.basename(cFile).split('_')[0]
dMapperName=os.path.basename(dFile).split('_')[0]
eMapperName=os.path.basename(eFile).split('_')[0]

consisCpgFileName=outPut+"_fiveDetec_consisCpg.bed"
nonConsisCpgFileName=outPut+"_fiveDetec_nonConsisCpg.bed"
diffNonConsisCpgFileName=outPut+"_fiveDiff_all.bed"
aNonConsisCpgFileName=outPut+"_"+aMapperName+"_diffDetecCpg.bed"
bNonConsisCpgFileName=outPut+"_"+bMapperName+"_diffDetecCpg.bed"
cNonConsisCpgFileName=outPut+"_"+cMapperName+"_diffDetecCpg.bed"
dNonConsisCpgFileName=outPut+"_"+dMapperName+"_diffDetecCpg.bed"
eNonConsisCpgFileName=outPut+"_"+eMapperName+"_diffDetecCpg.bed"

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

def GgiAndGeneAnnoResSta(consisCpgAnnoResBed,nonConsisCpgAnnoResBed,diffNonConsisCpgAnnoResBed,aNonConsisCpgAnnoResBed,bNonConsisCpgAnnoResBed,cNonConsisCpgAnnoResBed,dNonConsisCpgAnnoResBed,eNonConsisCpgAnnoResBed,outFileName,lineNum):
    outFile=open(outFileName,"w")
    
    consisCpgDict={}
    nonConsisCpgDict={}
    diffNonConsisCpgDict={}
    aNonConsisCpgDict={}
    bNonConsisCpgDict={}
    cNonConsisCpgDict={}
    dNonConsisCpgDict={}
    eNonConsisCpgDict={}
    
    consisCpgDict=resSta(consisCpgAnnoResBed,consisCpgDict,lineNum)
    nonConsisCpgDict=resSta(nonConsisCpgAnnoResBed,nonConsisCpgDict,lineNum)
    diffNonConsisCpgDict=resSta(diffNonConsisCpgAnnoResBed,diffNonConsisCpgDict,lineNum)
    aNonConsisCpgDict=resSta(aNonConsisCpgAnnoResBed,aNonConsisCpgDict,lineNum)
    bNonConsisCpgDict=resSta(bNonConsisCpgAnnoResBed,bNonConsisCpgDict,lineNum)
    cNonConsisCpgDict=resSta(cNonConsisCpgAnnoResBed,cNonConsisCpgDict,lineNum)
    dNonConsisCpgDict=resSta(dNonConsisCpgAnnoResBed,dNonConsisCpgDict,lineNum)
    eNonConsisCpgDict=resSta(eNonConsisCpgAnnoResBed,eNonConsisCpgDict,lineNum)
    
    allClassSet=set(consisCpgDict.keys()) | set(nonConsisCpgDict.keys()) | set(diffNonConsisCpgDict.keys()) | set(aNonConsisCpgDict.keys()) | set(bNonConsisCpgDict.keys()) | set(cNonConsisCpgDict.keys()) | set(dNonConsisCpgDict.keys()) | set(eNonConsisCpgDict.keys())
    
    consisCpgAnnoResBedFileName=os.path.basename(consisCpgAnnoResBed).split(".")[0]
    nonConsisCpgAnnoResBedFileName=os.path.basename(nonConsisCpgAnnoResBed).split(".")[0]
    diffNonConsisCpgAnnoResBedFileName=os.path.basename(diffNonConsisCpgAnnoResBed).split(".")[0]
    aNonConsisCpgAnnoResBedFileName=os.path.basename(aNonConsisCpgAnnoResBed).split(".")[0]
    bNonConsisCpgAnnoResBedFileName=os.path.basename(bNonConsisCpgAnnoResBed).split(".")[0]
    cNonConsisCpgAnnoResBedFileName=os.path.basename(cNonConsisCpgAnnoResBed).split(".")[0]
    dNonConsisCpgAnnoResBedFileName=os.path.basename(dNonConsisCpgAnnoResBed).split(".")[0]
    eNonConsisCpgAnnoResBedFileName=os.path.basename(eNonConsisCpgAnnoResBed).split(".")[0]
    
    outFile.write("resClass"+"\t"+consisCpgAnnoResBedFileName+"\t"+nonConsisCpgAnnoResBedFileName+"\t"+\
                  diffNonConsisCpgAnnoResBedFileName+"\t"+aNonConsisCpgAnnoResBedFileName+"\t"+\
                     bNonConsisCpgAnnoResBedFileName+"\t"+cNonConsisCpgAnnoResBedFileName+"\t"+\
                         dNonConsisCpgAnnoResBedFileName+"\t"+eNonConsisCpgAnnoResBedFileName+"\n")
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
        if resClass in eNonConsisCpgDict:
            eNonConsisCpg = eNonConsisCpgDict[resClass]
        else:
            eNonConsisCpg = 0
        outFile.write(resClass+"\t"+str(consisCpg)+"\t"+str(nonConsisCpg)+"\t"+str(diffNonConsisCpg)+"\t"+
                      str(aNonConsisCpg)+"\t"+str(bNonConsisCpg)+"\t"+str(cNonConsisCpg)+"\t"+str(dNonConsisCpg)+"\t"+
                      str(eNonConsisCpg)+"\n")


## Cpg realated region annotation results of cpg locus
consisCpgAnnoResBedCGI=outPut+"_fiveDetec_consisCpgAnnoResCGI.bed"
nonConsisCpgAnnoResBedCGI=outPut+"_fiveDetec_nonConsisCpgAnnoResCGI.bed"
diffNonConsisCpgAnnoResBedCGI=outPut+"_fiveDiff_allAnnoResCGI.bed"
aNonConsisCpgAnnoResBedCGI=outPut+"_"+aMapperName+"_diffDetecCpgAnnoResCGI.bed"
bNonConsisCpgAnnoResBedCGI=outPut+"_"+bMapperName+"_diffDetecCpgAnnoResCGI.bed"
cNonConsisCpgAnnoResBedCGI=outPut+"_"+cMapperName+"_diffDetecCpgAnnoResCGI.bed"
dNonConsisCpgAnnoResBedCGI=outPut+"_"+dMapperName+"_diffDetecCpgAnnoResCGI.bed"
eNonConsisCpgAnnoResBedCGI=outPut+"_"+eMapperName+"_diffDetecCpgAnnoResCGI.bed"

os.system('bedtools intersect -a %s -b %s -wao > %s' %(consisCpgFileName,cgiAnnoFile,consisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(nonConsisCpgFileName,cgiAnnoFile,nonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(diffNonConsisCpgFileName,cgiAnnoFile,diffNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(aNonConsisCpgFileName,cgiAnnoFile,aNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(bNonConsisCpgFileName,cgiAnnoFile,bNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(cNonConsisCpgFileName,cgiAnnoFile,cNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(dNonConsisCpgFileName,cgiAnnoFile,dNonConsisCpgAnnoResBedCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(eNonConsisCpgFileName,cgiAnnoFile,eNonConsisCpgAnnoResBedCGI))

resStaOfCgiAnno=outPut+"resStaOfCgiAnno.txt"
GgiAndGeneAnnoResSta(consisCpgAnnoResBedCGI,nonConsisCpgAnnoResBedCGI,diffNonConsisCpgAnnoResBedCGI,aNonConsisCpgAnnoResBedCGI,bNonConsisCpgAnnoResBedCGI,cNonConsisCpgAnnoResBedCGI,dNonConsisCpgAnnoResBedCGI,eNonConsisCpgAnnoResBedCGI,resStaOfCgiAnno,7)

## repeat region annotation results of cpg locus
consisCpgAnnoResBedRep=outPut+"_fiveDetec_consisCpgAnnoResRep.bed"
nonConsisCpgAnnoResBedRep=outPut+"_fiveDetec_nonConsisCpgAnnoResRep.bed"
diffNonConsisCpgAnnoResBedRep=outPut+"_fiveDiff_allAnnoResRep.bed"
aNonConsisCpgAnnoResBedRep=outPut+"_"+aMapperName+"_diffDetecCpgAnnoResRep.bed"
bNonConsisCpgAnnoResBedRep=outPut+"_"+bMapperName+"_diffDetecCpgAnnoResRep.bed"
cNonConsisCpgAnnoResBedRep=outPut+"_"+cMapperName+"_diffDetecCpgAnnoResRep.bed"
dNonConsisCpgAnnoResBedRep=outPut+"_"+dMapperName+"_diffDetecCpgAnnoResRep.bed"
eNonConsisCpgAnnoResBedRep=outPut+"_"+eMapperName+"_diffDetecCpgAnnoResRep.bed"

os.system('bedtools intersect -a %s -b %s -wao > %s' %(consisCpgFileName,repAnnoFile,consisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(nonConsisCpgFileName,repAnnoFile,nonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(diffNonConsisCpgFileName,repAnnoFile,diffNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(aNonConsisCpgFileName,repAnnoFile,aNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(bNonConsisCpgFileName,repAnnoFile,bNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(cNonConsisCpgFileName,repAnnoFile,cNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(dNonConsisCpgFileName,repAnnoFile,dNonConsisCpgAnnoResBedRep))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(eNonConsisCpgFileName,repAnnoFile,eNonConsisCpgAnnoResBedRep))

resStaOfRepAnno=outPut+"resStaOfRepAnno.txt"
GgiAndGeneAnnoResSta(consisCpgAnnoResBedRep,nonConsisCpgAnnoResBedRep,diffNonConsisCpgAnnoResBedRep,aNonConsisCpgAnnoResBedRep,bNonConsisCpgAnnoResBedRep,cNonConsisCpgAnnoResBedRep,dNonConsisCpgAnnoResBedRep,eNonConsisCpgAnnoResBedRep,resStaOfRepAnno,9)


## Distribution of methylation levels of consistent CpG sites and non-identical CpG sites
consisCpgMethFileName=outPut+"_fiveDetect_consisCpgMeth.bed"
nonConsisCpgMethFileName=outPut+"_fiveDetec_nonConsisCpgMeth.bed"
resStaOfMethFileName=outPut+"cpgMehtLevelDistri.txt"
aMapperName=os.path.basename(aFile).split('_')[0]
bMapperName=os.path.basename(bFile).split('_')[0]
cMapperName=os.path.basename(cFile).split('_')[0]
dMapperName=os.path.basename(dFile).split('_')[0]
eMapperName=os.path.basename(eFile).split('_')[0]

resStaOfMeth=open(resStaOfMethFileName,"w")
resStaOfMeth.write("class"+"\t"+"highMeth"+"\t"+"mediumMeth"+"\t"+"lowMeth"+"\n")

def methClass(fileName,resStaOfMeth,className):
    file=open(fileName,"r")
    highMethNum=0
    mediumMethNum=0
    lowMethNum=0
    for line in file:
        methLevel=float(line.split()[3])
        if methLevel >= 2/3:
            highMethNum += 1
        elif methLevel < 1/3:
            lowMethNum += 1
        else:
            mediumMethNum += 1
    resStaOfMeth.write(className+"\t"+str(highMethNum)+"\t"+str(mediumMethNum)+"\t"+str(lowMethNum)+"\n")

methClass(consisCpgMethFileName,resStaOfMeth,"consis")
methClass(nonConsisCpgMethFileName,resStaOfMeth,"nonConsis")
methClass(aNonConsisCpgFileName,resStaOfMeth,aMapperName)
methClass(bNonConsisCpgFileName,resStaOfMeth,bMapperName)
methClass(cNonConsisCpgFileName,resStaOfMeth,cMapperName)
methClass(dNonConsisCpgFileName,resStaOfMeth,dMapperName)
methClass(eNonConsisCpgFileName,resStaOfMeth,eMapperName)

