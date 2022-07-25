# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:31:35 2022

@author: 滔
"""

import threading
import argparse
import os
from subprocess import Popen,PIPE
from operator import itemgetter
import gzip
parser = argparse.ArgumentParser(description='Analysis of incorrect and unique comparison reads')
parser.add_argument('-a', '--afile', help='Input a file,end by .sam or .bam', required=True)
parser.add_argument('-b', '--bfile', help='Input b file,end by .sam or .bam', required=True)
parser.add_argument('-c', '--cfile', help='Input c file,end by .sam or .bam', required=True)
parser.add_argument('-d', '--dfile', help='Input d file,end by .sam or .bam', required=True)
parser.add_argument('-e', '--efile', help='Input e file,end by .sam or .bam', required=True)
parser.add_argument('-an', '--afileName', help='Input a name of afile result,such as bismarkbwt2', required=True)
parser.add_argument('-bn', '--bfileName', help='Input a name of bfile result,such as bsmap', required=True)
parser.add_argument('-cn', '--cfileName', help='Input a name of cfile result,such as bwameth', required=True)
parser.add_argument('-dn', '--dfileName', help='Input a name of dfile result,such as walt', required=True)
parser.add_argument('-en', '--efileName', help='Input a name of dfile result,such as bsbolt', required=True)
parser.add_argument('-bis', '--bismarkbwt2', help='Input a file of bismarkbwt2 ambiguouse,end by fq.gz', required=True)
parser.add_argument('-f', '--fastq', help='Input a file of fastq,end by .fq or .fastq', required=True)
parser.add_argument('-s', '--chromSize', help='Input a file of chrom size,end by .txt', required=True)
parser.add_argument('-er', '--errorReport', help='Statistical report of abnormal data.end by .txt', required=True)
parser.add_argument('-o', '--output', help='Folder to store output results,such as /path', required=True)
args = parser.parse_args()

## Parameter naming
aFile = args.afile
bFile = args.bfile
cFile = args.cfile
dFile = args.dfile
eFile = args.efile
aResName= args.afileName
bResName= args.bfileName
cResName= args.cfileName
dResName= args.dfileName
eResName= args.efileName
bismarkbwt2Mul=args.bismarkbwt2
fqFile = args.fastq
chromSizeFile = args.chromSize
errorReportFile = args.errorReport
outputPath = args.output

## Create a folder to store the results
tempFolder=outputPath+"/temp"
os.system("mkdir -p %s" %tempFolder)

readIdFolder=outputPath+"/readId"
os.system("mkdir -p %s" %readIdFolder)

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

## Identify the input file format(sam or bam)
def samOrBam(infile_name):
	if os.path.splitext(infile_name)[1] == '.sam':
		infile=open(infile_name,'r')
		return(infile)
	elif os.path.splitext(infile_name)[1] == '.bam':
		comBam="samtools view -h %s" % infile_name
		p=Popen(comBam,stdout=PIPE,shell=True,universal_newlines=True)
		infile=p.stdout
		return(infile)	

def wcL(file_name):
	comBam="wc -l %s" % file_name
	p=Popen(comBam,stdout=PIPE,shell=True,universal_newlines=True)
	infile=p.stdout
	num=int(infile.readline().split()[0])
	return(num)

## Get all readID
def creatChromSizeDict(chromSize):
  chromSizeDict={}
  with open(chromSize,"r") as f:
    for line in f:
      chromSizeDict[line.split()[0]] = int(line.split()[1])
  return(chromSizeDict)
  
def extractReadID(inFileName,outFileName,errorReportFile,chromSizeFileName):
    fragLenLess150Num = 0
    startCoorNegativeNum = 0
    fragEndMoreChromEndNum = 0
    chromSizeDict=creatChromSizeDict(chromSizeFileName)
    if inFileName[-3:]==".gz":
        inFile=gzip.open(inFileName,"r")
    else:
        inFile=open(inFileName,"r")
    outFile=open(outFileName,"w")
    errorReport =open(errorReportFile,"w")
    for line in inFile:
        if inFileName[-3:]==".gz":
            line=line.decode()
        if line.startswith('@'):
            line=line.strip('@\n')
            chrom = line.split(":")[0].split("_",1)[1]
            chromEnd = chromSizeDict[chrom]
            if line.split(":")[1].split("-")[0] == "":
                startCoorNegativeNum += 1
                continue
            elif int(line.split(":")[1].split("-")[-1]) > chromEnd:
                fragEndMoreChromEndNum += 1
                continue
            else:
                r1ReadID=line+"_R1\n"
                r2ReadID=line+"_R2\n"
                outFile.write(r1ReadID+r2ReadID)
        else:
            continue
    errorReport.write("startCoorNegativeNum"+"\t"+str(startCoorNegativeNum)+"\n"+
										"fragLenLess150Num"+"\t"+str(fragLenLess150Num)+"\n"+
										"fragEndMoreChromEndNum"+"\t"+str(fragEndMoreChromEndNum)+"\n")
	
			
allReadIdFile=tempFolder+"/allReadId.txt"
extractReadID(fqFile,allReadIdFile,errorReportFile,chromSizeFile)


## Get the readID correctly mapped by each software
def readIdRightError(inFile,rightMap,errorMap,chromSizeFileName):
    samFile=samOrBam(inFile)
    rightFile=open(rightMap,"w")
    errorFile=open(errorMap,"w")
    chromSizeDict=creatChromSizeDict(chromSizeFileName)
    for line in samFile:
        if line.startswith('@'):
            continue
        else:
            bitFlag = int(line.split('\t')[1])
            predChrom = line.split('\t')[2]
            trueChrom = line.split('\t')[0].split(':')[0].split('_',1)[1]
            readId=line.split('\t')[0]
            chromEnd = chromSizeDict[trueChrom]
            if readId.split(":")[1].split("-")[0] == "":
                continue
            elif int(readId.split(":")[1].split("-")[-1]) > chromEnd:
                continue
            else:    
			  #first read forward strand (or case [67]: bsseeker2)
                if (bitFlag == 99 or bitFlag == 67):
                    predPosition = int(line.split('\t')[3])
                    if line.split('\t')[0].split(':')[1].split('-')[0].split('_')[0]!='':
                        truePosition = int(line.split('\t')[0].split(':')[1].split('-')[0].split('_')[0])
                    if (predChrom == trueChrom and predPosition==truePosition):
                        rightFile.write(readId+"_R1"+"\n")
                    elif (predChrom != trueChrom or predPosition!=truePosition):
                        errorFile.write(readId+"_R1"+"\n")
                        
			  #secound read reverse strand (or case [131]: bsseeker2)
                elif (bitFlag == 147 or bitFlag == 131):
                    predPosition = int(line.split('\t')[3])
                    if line.split('\t')[0].split(':')[1].split('-')[1].split('_')[0]!='':
                        truePosition = int(line.split('\t')[0].split(':')[1].split('-')[1].split('_')[0])-150+1
                    if (predChrom == trueChrom and predPosition==truePosition):
                        rightFile.write(readId+"_R2"+"\n")
                    elif (predChrom != trueChrom or predPosition!=truePosition):
                        errorFile.write(readId+"_R2"+"\n")
                        
			     #first read reverse strand
                elif (bitFlag == 83 or bitFlag == 115):
                    predPosition = int(line.split('\t')[3])
                    if line.split('\t')[0].split(':')[1].split('-')[1].split('_')[0]!='':
                        truePosition = int(line.split('\t')[0].split(':')[1].split('-')[1].split('_')[0])-150+1
                    if (predChrom == trueChrom and predPosition==truePosition):
                        rightFile.write(readId+"_R1"+"\n")
                    elif (predChrom != trueChrom or predPosition!=truePosition):
                        errorFile.write(readId+"_R1"+"\n")
			    #secound read forward strand
                elif (bitFlag == 163 or bitFlag == 179):
                    predPosition = int(line.split('\t')[3])
                    if line.split('\t')[0].split(':')[1].split('-')[0].split('_')[0]!='':
                        truePosition = int(line.split('\t')[0].split(':')[1].split('-')[0].split('_')[0]) 
                    if (predChrom == trueChrom and predPosition==truePosition):
                        rightFile.write(readId+"_R2"+"\n")
                    elif (predChrom != trueChrom or predPosition!=truePosition):
                        errorFile.write(readId+"_R2"+"\n")
			#read doesn't match
                else:
                    continue
    

aRightMapReadId=tempFolder+"/"+aResName+"RightMapReadId.txt"
bRightMapReadId=tempFolder+"/"+bResName+"RightMapReadId.txt"
cRightMapReadId=tempFolder+"/"+cResName+"RightMapReadId.txt"
dRightMapReadId=tempFolder+"/"+dResName+"RightMapReadId.txt"
eRightMapReadId=tempFolder+"/"+eResName+"RightMapReadId.txt"

aErrorMapReadId=tempFolder+"/"+aResName+"ErrorMapReadId.txt"
bErrorMapReadId=tempFolder+"/"+bResName+"ErrorMapReadId.txt"
cErrorMapReadId=tempFolder+"/"+cResName+"ErrorMapReadId.txt"
dErrorMapReadId=tempFolder+"/"+dResName+"ErrorMapReadId.txt"
eErrorMapReadId=tempFolder+"/"+eResName+"ErrorMapReadId.txt"

threads=[]
AreadIdRightError=MyThread(readIdRightError,args=(aFile,aRightMapReadId,aErrorMapReadId,chromSizeFile,))
BreadIdRightError=MyThread(readIdRightError,args=(bFile,bRightMapReadId,bErrorMapReadId,chromSizeFile,))
CreadIdRightError=MyThread(readIdRightError,args=(cFile,cRightMapReadId,cErrorMapReadId,chromSizeFile,))
DreadIdRightError=MyThread(readIdRightError,args=(dFile,dRightMapReadId,dErrorMapReadId,chromSizeFile,))
EreadIdRightError=MyThread(readIdRightError,args=(eFile,eRightMapReadId,eErrorMapReadId,chromSizeFile,))
threads.append(AreadIdRightError)
threads.append(BreadIdRightError)
threads.append(CreadIdRightError)
threads.append(DreadIdRightError)
threads.append(EreadIdRightError)
for thr in threads:
    thr.start()
for thr in threads:
    if thr.is_alive():
        thr.join()

aErrorMapNum=wcL(aErrorMapReadId)
bErrorMapNum=wcL(bErrorMapReadId)
cErrorMapNum=wcL(cErrorMapReadId)
dErrorMapNum=wcL(dErrorMapReadId)
eErrorMapNum=wcL(eErrorMapReadId)


## Incorrect unique mapped readID
def cutFiles(infile,size,fileList,tempPath):
	with open(infile,"r") as f:
		for line in f:
			readnum=int(line.split("_")[0])
			fileid=readnum//size+1
			if fileid in fileList:
				filename=''.join((tempPath,"temp",str(fileid),".txt"))
				with open(filename,"a") as df:
					df.write(line)
			else:
				fileList.append(fileid)
				filename=''.join((tempPath,"temp",str(fileid),".txt"))   
				with open(filename,"w") as df:
					df.write(line)
	return(fileList)
	
def diff_seta_setb(setaFileName,setbFileName,tempPath,outFileName):
  size=1000000
  fileList=[]
  infilelist=["a","b"]
  forrunnum=0
  NonRightMapNum=0
  atempfilename=''.join((tempPath,"/atemp/"))
  if not os.path.exists(atempfilename):
     os.mkdir(atempfilename)
  fileList=cutFiles(setaFileName,size,fileList,atempfilename)
	
  btempfilename=''.join((tempPath,"/btemp/"))
  if not os.path.exists(btempfilename):
    os.mkdir(btempfilename)
  fileList=cutFiles(setbFileName,size,fileList,btempfilename)
	
  for fileid in fileList:
    forrunnum += 1
    Aset=set()
    Bset=set()
    for infile in infilelist:
      if infile=="a":
        filename=''.join((atempfilename,"temp",str(fileid),".txt"))
        if os.path.exists(filename):
          with open(filename,"r") as f:
            for line in f:
              Aset.add(line)
      if infile=="b":
        filename=''.join((btempfilename,"temp",str(fileid),".txt"))
        if os.path.exists(filename):
          with open(filename,"r") as f:
            for line in f:
              Bset.add(line)
    if forrunnum==1:
       A_Bset=Aset.difference(Bset)
       with open(outFileName,"w") as f:
         for line in A_Bset:
           NonRightMapNum += 1
           f.write(line)
    else:
       A_Bset=Aset.difference(Bset)
       with open(outFileName,"a") as f:
         for line in A_Bset:
           NonRightMapNum += 1
           f.write(line)
    os.system("rm -r %s" % atempfilename)
    os.system("rm -r %s" % btempfilename)
    return(NonRightMapNum)	
	
aNonRightMapReadId=tempFolder+"/"+aResName+"NonRightMapReadId.txt"
bNonRightMapReadId=tempFolder+"/"+bResName+"NonRightMapReadId.txt"
cNonRightMapReadId=tempFolder+"/"+cResName+"NonRightMapReadId.txt"
dNonRightMapReadId=tempFolder+"/"+dResName+"NonRightMapReadId.txt"
eNonRightMapReadId=tempFolder+"/"+eResName+"NonRightMapReadId.txt"

os.system("awk 'NR==FNR{ a[$1]=$1 } NR>FNR{ if(a[$1] == \"\"){ print $1}}' %s %s > %s" %(aRightMapReadId,allReadIdFile,aNonRightMapReadId))
os.system("awk 'NR==FNR{ a[$1]=$1 } NR>FNR{ if(a[$1] == \"\"){ print $1}}' %s %s > %s" %(bRightMapReadId,allReadIdFile,bNonRightMapReadId))
os.system("awk 'NR==FNR{ a[$1]=$1 } NR>FNR{ if(a[$1] == \"\"){ print $1}}' %s %s > %s" %(cRightMapReadId,allReadIdFile,cNonRightMapReadId))
os.system("awk 'NR==FNR{ a[$1]=$1 } NR>FNR{ if(a[$1] == \"\"){ print $1}}' %s %s > %s" %(dRightMapReadId,allReadIdFile,dNonRightMapReadId))
os.system("awk 'NR==FNR{ a[$1]=$1 } NR>FNR{ if(a[$1] == \"\"){ print $1}}' %s %s > %s" %(eRightMapReadId,allReadIdFile,eNonRightMapReadId))


aNonRightMapNum=wcL(aNonRightMapReadId)
bNonRightMapNum=wcL(bNonRightMapReadId)
cNonRightMapNum=wcL(cNonRightMapReadId)
dNonRightMapNum=wcL(dNonRightMapReadId)
eNonRightMapNum=wcL(eNonRightMapReadId)


## Incorrect unique mapped readID of common and difference
def intersecStaResOut(afile,bfile,cfile,dfile,efile,aResName,bResName,cResName,dResName,eResName,tempPath,reportResFileName):
 size=1000000
 Anum=0
 Bnum=0
 Cnum=0
 Dnum=0
 Enum=0
 ABnum=0
 ACnum=0
 ADnum=0
 AEnum=0
 BCnum=0
 BDnum=0
 BEnum=0
 CDnum=0
 CEnum=0
 DEnum=0
 ABCnum=0
 ABDnum=0
 ABEnum=0
 ACDnum=0
 ACEnum=0
 ADEnum=0
 BCDnum=0
 BCEnum=0
 BDEnum=0
 CDEnum=0
 ABCDnum=0
 ABCEnum=0
 ABDEnum=0
 ACDEnum=0
 BCDEnum=0
 ABCDEnum=0
 forrunnum=0
 filelist=[]
 mapperlist=["a","b","c","d","e"]
 reportResFile = open(reportResFileName,"w")
	
 atempfilename=''.join((tempPath,"/atemp/"))
 if not os.path.exists(atempfilename):
   os.mkdir(atempfilename)
 filelist=cutFiles(afile,size,filelist,atempfilename)

 btempfilename=''.join((tempPath,"/btemp/"))
 if not os.path.exists(btempfilename):
   os.mkdir(btempfilename)
 filelist=cutFiles(bfile,size,filelist,btempfilename)	
					
 ctempfilename=''.join((tempPath,"/ctemp/"))
 if not os.path.exists(ctempfilename):
   os.mkdir(ctempfilename)
 filelist=cutFiles(cfile,size,filelist,ctempfilename)	

 dtempfilename=''.join((tempPath,"/dtemp/"))
 if not os.path.exists(dtempfilename):
   os.mkdir(dtempfilename)
 filelist=cutFiles(dfile,size,filelist,dtempfilename)
	
 etempfilename=''.join((tempPath,"/etemp/"))
 if not os.path.exists(etempfilename):
   os.mkdir(etempfilename)
 filelist=cutFiles(efile,size,filelist,etempfilename)
    
 for fileid in filelist:
   forrunnum += 1
   Aset=set()
   Bset=set()
   Cset=set()
   Dset=set()
   Eset=set()
   for mapper in mapperlist:
     if mapper == "a":
       filename=''.join((atempfilename,"temp",str(fileid),".txt"))
       if os.path.exists(filename):
         with open(filename,"r") as f:
           for line in f:
             Anum += 1
             Aset.add(line)
     if mapper == "b":
       filename=''.join((btempfilename,"temp",str(fileid),".txt"))
       if os.path.exists(filename):
         with open(filename,"r") as f:
           for line in f:
             Bnum += 1
             Bset.add(line)
     if mapper =="c":
       filename=''.join((ctempfilename,"temp",str(fileid),".txt"))
       if os.path.exists(filename):
         with open(filename,"r") as f:
           for line in f:
             Cnum += 1
             Cset.add(line)
     if mapper == "d":
       filename=''.join((dtempfilename,"temp",str(fileid),".txt"))
       if os.path.exists(filename):
          with open(filename,"r") as f:
            for line in f:
              Dnum += 1
              Dset.add(line)
     if mapper == "e":
       filename=''.join((etempfilename,"temp",str(fileid),".txt"))
       if os.path.exists(filename):
          with open(filename,"r") as f:
            for line in f:
              Enum += 1
              Eset.add(line)

   ABset=Aset.intersection(Bset)
   ABnum += len(ABset)
			
   ACset=Aset.intersection(Cset)
   ACnum += len(ACset)
			
   ADset=Aset.intersection(Dset)
   ADnum += len(ADset)
       
   AEset=Aset.intersection(Eset)
   AEnum += len(AEset)
       
   BCset=Bset.intersection(Cset)
   BCnum += len(BCset)
		
   BDset=Bset.intersection(Dset)
   BDnum += len(BDset)
       
   BEset=Bset.intersection(Eset)
   BEnum += len(BEset)
              	
   CDset=Cset.intersection(Dset)
   CDnum += len(CDset)
       
   CEset=Cset.intersection(Eset)
   CEnum += len(CEset)
       
   DEset=Dset.intersection(Eset)
   DEnum += len(DEset)
       
   ABCset=Aset.intersection(BCset)
   ABCnum += len(ABCset)
       
   ABDset=Aset.intersection(BDset)
   ABDnum += len(ABDset)
       
   ABEset=Aset.intersection(BEset)
   ABEnum += len(ABEset)
		
   ACDset=Aset.intersection(CDset)
   ACDnum += len(ACDset)
       
   ACEset=Aset.intersection(CEset)
   ACEnum += len(ACEset)
       
   ADEset=Aset.intersection(DEset)
   ADEnum += len(ADEset)
        
   BCDset=Bset.intersection(CDset)
   BCDnum += len(BCDset)
       
   BCEset=Bset.intersection(CEset)
   BCEnum += len(BCEset)
       
   BDEset=Bset.intersection(DEset)
   BDEnum += len(BDEset)
       
   CDEset=Cset.intersection(DEset)
   CDEnum += len(CDEset)

   ABCDset=Aset.intersection(BCDset)
   ABCDnum += len(ABCDset)
       
   ABCEset=Aset.intersection(BCEset)
   ABCEnum += len(ABCEset)
       
   ABDEset=Aset.intersection(BDEset)
   ABDEnum += len(ABDEset)
       
   ACDEset=Aset.intersection(CDEset)
   ACDEnum += len(ACDEset)
       
   BCDEset=Bset.intersection(CDEset)
   BCDEnum += len(BCDEset)
       
   ABCDEset=Aset.intersection(BCDEset)
   ABCDEnum += len(ABCDEset)
 reportResFile.write("A\t"+aResName+"\n"+
                           "B\t"+bResName+"\n"+
                           "C\t"+cResName+"\n"+
                           "D\t"+dResName+"\n"+
                           "E\t"+eResName+"\n"+
                           "Anum\t"+str(Anum)+"\n"+
                           "Bnum\t"+str(Bnum)+"\n"+
                           "Cnum\t"+str(Cnum)+"\n"+
                           "Dnum\t"+str(Dnum)+"\n"+
                           "Enum\t"+str(Enum)+"\n"+
                           "ABnum\t"+str(ABnum)+"\n"+
                           "ACnum\t"+str(ACnum)+"\n"+
                           "ADnum\t"+str(ADnum)+"\n"+
                           "AEnum\t"+str(AEnum)+"\n"+
                           "BCnum\t"+str(BCnum)+"\n"+
                           "BDnum\t"+str(BDnum)+"\n"+
                           "BEnum\t"+str(BEnum)+"\n"+
                           "CDnum\t"+str(CDnum)+"\n"+
                           "CEnum\t"+str(CEnum)+"\n"+
                           "DEnum\t"+str(DEnum)+"\n"+
                           "ABCnum\t"+str(ABCnum)+"\n"+
                           "ABDnum\t"+str(ABDnum)+"\n"+
                           "ABEnum\t"+str(ABEnum)+"\n"+
                           "ACDnum\t"+str(ACDnum)+"\n"+
                           "ACEnum\t"+str(ACEnum)+"\n"+
                           "ADEnum\t"+str(ADEnum)+"\n"+
                           "BCDnum\t"+str(BCDnum)+"\n"+
                           "BCEnum\t"+str(BCEnum)+"\n"+
                           "BDEnum\t"+str(BDEnum)+"\n"+
                           "CDEnum\t"+str(CDEnum)+"\n"+
                           "ABCDnum\t"+str(ABCDnum)+"\n"+
                           "ABCEnum\t"+str(ABCEnum)+"\n"+
                           "ABDEnum\t"+str(ABDEnum)+"\n"+
                           "ACDEnum\t"+str(ACDEnum)+"\n"+
                           "BCDEnum\t"+str(BCDEnum)+"\n"+
                           "ABCDEnum\t"+str(ABCDEnum)+"\n")                           
 os.system("rm -r %s" %atempfilename)
 os.system("rm -r %s" %btempfilename)
 os.system("rm -r %s" %ctempfilename)
 os.system("rm -r %s" %dtempfilename)
 os.system("rm -r %s" %etempfilename)


nonRightMapReadIdReport=readIdFolder+"/nonRightMapReadIdReport.txt"
intersecStaResOut(afile=aNonRightMapReadId,\
                  bfile=bNonRightMapReadId,\
                  cfile=cNonRightMapReadId,\
                  dfile=dNonRightMapReadId,\
                  efile=eNonRightMapReadId,\
                  aResName=aResName,\
                  bResName=bResName,\
                  cResName=cResName,\
                  dResName=dResName,\
                  eResName=eResName,\
				  tempPath=tempFolder,\
				  reportResFileName=nonRightMapReadIdReport)
									
## mapping class of nonRight mapped reads
def read1And2Match(readIdDict,outFileName,num,multMapNum):
	output=0
	if num == 1:
		outFile=open(outFileName,'w')
	else:
		outFile=open(outFileName,'a')
	if readIdDict['339'] != [] and readIdDict['419'] != []:
		readIdDict['339'].sort(key=itemgetter(2,3,7))
		readIdDict['419'].sort(key=itemgetter(2,7,3))
		iRead1,iRead2=0,0
		while iRead1 < len(readIdDict['339']) and iRead2 < len(readIdDict['419']):
			r1Line,r2Line=readIdDict['339'][iRead1],readIdDict['419'][iRead2]
			if r1Line[2] == r2Line[2]:
				if r1Line[3] == r2Line[7]:
					if r1Line[7] == r2Line[3]:
						iRead1 += 1
						iRead2 += 1
						multMapNum += 2
						outFile.write(r1Line[0]+'_R1\n'+r2Line[0]+'_R2\n')
						output += 1
						break
					elif r1Line[7] > r2Line[3]:
						iRead2 += 1
					elif r1Line[7] < r2Line[3]:
						iRead1 += 1
				elif r1Line[3] > r2Line[7]:
					iRead2 += 1
				elif r1Line[3] < r2Line[7]:
					iRead1 += 1
			elif r1Line[2] > r2Line[2]:
				iRead2 += 1
			elif r1Line[2] < r2Line[2]:
				iRead1 += 1	
	elif readIdDict['355'] != [] and readIdDict['403'] != [] and output == 0:
		readIdDict['355'].sort(key=itemgetter(2,3,7))
		readIdDict['403'].sort(key=itemgetter(2,7,3))
		iRead1,iRead2=0,0
		while iRead1 < len(readIdDict['355']) and iRead2 < len(readIdDict['403']):
			r1Line,r2Line=readIdDict['355'][iRead1],readIdDict['403'][iRead2]
			if r1Line[2] == r2Line[2]:
				if r1Line[3] == r2Line[7]:
					if r1Line[7] == r2Line[3]:
						iRead1 += 1
						iRead2 += 1
						multMapNum += 2
						outFile.write(r1Line[0]+'_R1\n'+r2Line[0]+'_R2\n')
						break
					elif r1Line[7] > r2Line[3]:
						iRead2 += 1
					elif r1Line[7] < r2Line[3]:
						iRead1 += 1
				elif r1Line[3] > r2Line[7]:    
					iRead2 += 1
				elif r1Line[3] < r2Line[7]:
					iRead1 += 1
			elif r1Line[2] > r2Line[2]:
				iRead2 += 1
			elif r1Line[2] < r2Line[2]:
				iRead1 += 1
	return(multMapNum)

def multMapReadIdExtract(bamOrSamFile,outFileName,chromSizeFileName):
    com="samtools view %s|awk 'END{print NR}'" % bamOrSamFile
    p=Popen(com,stdout=PIPE,shell=True,universal_newlines=True)
    maxLineNum=int(p.stdout.read().strip())
    samFile=samOrBam(bamOrSamFile)
#    samFile=open(bamOrSamFile,"r")
    readId=''
    lineNum=0
    num=0
    multMapNum=0
    readIdDict={'339':[],'355':[],'419':[],'403':[]}
    chromSizeDict=creatChromSizeDict(chromSizeFileName)
    for line in samFile:
        if line.startswith('@'):
            continue
        else:
            lineNum += 1
            line=line.split('\t')
            
            chromEnd=chromSizeDict[line[0].split(":")[0].split("_",1)[1]]
            
            if line[0].split(":")[1].split("-")[0] == "":
                continue
            elif int(line[0].split(":")[1].split("-")[-1]) > chromEnd:
                continue
            else:
                lineFlags=line[1]
                if lineNum == maxLineNum:
                    if lineFlags in readIdDict.keys():
                        readIdDict[lineFlags].append(line)
                    multMapNum=read1And2Match(readIdDict,outFileName,num,multMapNum)
                else:
                    if lineFlags in readIdDict.keys():
                        lineReadId = line[0]
                        if readId=='':
                            readId=lineReadId
                            readIdDict[lineFlags].append(line)
                        else:
                            if lineReadId != readId:
                                num += 1
                                multMapNum=read1And2Match(readIdDict,outFileName,num,multMapNum)
                                readId=lineReadId
                                readIdDict={'339':[],'355':[],'419':[],'403':[]}
                                readIdDict[lineFlags].append(line)
                            elif lineReadId == readId:
                                readIdDict[lineFlags].append(line)
    return(multMapNum)

aMultMapReadId=tempFolder+"/"+aResName+"multMapReadId.txt"
bMultMapReadId=tempFolder+"/"+bResName+"multMapReadId.txt"
cMultMapReadId=tempFolder+"/"+cResName+"multMapReadId.txt"
dMultMapReadId=tempFolder+"/"+dResName+"multMapReadId.txt"
eMultMapReadId=tempFolder+"/"+eResName+"multMapReadId.txt"

def extractReadID2(inFileName,outFileName,chromSizeFileName):
    chromSizeDict=creatChromSizeDict(chromSizeFileName)
    multMapNum=0
    if inFileName[-3:]==".gz":
        inFile=gzip.open(inFileName,"r")
    else:
        inFile=open(inFileName,"r")
    outFile=open(outFileName,"w")
    for line in inFile:
        if inFileName[-3:]==".gz":
            line=line.decode()
        if line.startswith('@'):
            line=line.strip('@\n')
            chrom = line.split(":")[0].split("_",1)[1]
            chromEnd = chromSizeDict[chrom]
            if line.split(":")[1].split("-")[0] == "":
                continue
            elif int(line.split(":")[1].split("-")[-1]) > chromEnd:
                continue
            else:
                multMapNum += 2
                r1ReadID=line+"_R1\n"
                r2ReadID=line+"_R2\n"
                outFile.write(r1ReadID+r2ReadID)
        else:
            continue
    return(multMapNum)

amultMapNum=extractReadID2(bismarkbwt2Mul,aMultMapReadId,chromSizeFile)

threads=[]
BmultMapReadIdExtract=MyThread(multMapReadIdExtract,args=(bFile,bMultMapReadId,chromSizeFile,))
CmultMapReadIdExtract=MyThread(multMapReadIdExtract,args=(cFile,cMultMapReadId,chromSizeFile,))
DmultMapReadIdExtract=MyThread(multMapReadIdExtract,args=(dFile,dMultMapReadId,chromSizeFile,))
EmultMapReadIdExtract=MyThread(multMapReadIdExtract,args=(eFile,eMultMapReadId,chromSizeFile,))

threads.append(BmultMapReadIdExtract)
threads.append(CmultMapReadIdExtract)
threads.append(DmultMapReadIdExtract)
threads.append(EmultMapReadIdExtract)
for thr in threads:
    thr.start()
for thr in threads:
    if thr.is_alive():
        thr.join()
bmultMapNum=BmultMapReadIdExtract.get_result()
cmultMapNum=CmultMapReadIdExtract.get_result()
dmultMapNum=DmultMapReadIdExtract.get_result()
emultMapNum=EmultMapReadIdExtract.get_result()

aNonMapNum=aNonRightMapNum-aErrorMapNum-amultMapNum
bNonMapNum=bNonRightMapNum-bErrorMapNum-bmultMapNum
cNonMapNum=cNonRightMapNum-cErrorMapNum-cmultMapNum
dNonMapNum=dNonRightMapNum-dErrorMapNum-dmultMapNum
eNonMapNum=eNonRightMapNum-eErrorMapNum-emultMapNum

nonRightMapReadMappingClassReport=readIdFolder+"/"+"nonRightMapReadMappingClassReport.txt"

with open(nonRightMapReadMappingClassReport,"w") as df:
 df.write(aResName+"\t"+"aNonRightMapNum"+"\t"+str(aNonRightMapNum)+"\n"+
          aResName+"\t"+"aErrorMapNum"+"\t"+str(aErrorMapNum)+"\n"+
          aResName+"\t"+"aMultMapNum"+"\t"+str(amultMapNum)+"\n"+
          aResName+"\t"+"aNonMapNum"+"\t"+str(aNonMapNum)+"\n"+
          bResName+"\t"+"bNonRightMapNum"+"\t"+str(bNonRightMapNum)+"\n"+
          bResName+"\t"+"bErrorMapNum"+"\t"+str(bErrorMapNum)+"\n"+
          bResName+"\t"+"bMultMapNum"+"\t"+str(bmultMapNum)+"\n"+
          bResName+"\t"+"bNonMapNum"+"\t"+str(bNonMapNum)+"\n"+
          cResName+"\t"+"cNonRightMapNum"+"\t"+str(cNonRightMapNum)+"\n"+
          cResName+"\t"+"cErrorMapNum"+"\t"+str(cErrorMapNum)+"\n"+
          cResName+"\t"+"cMultMapNum"+"\t"+str(cmultMapNum)+"\n"+
          cResName+"\t"+"cNonMapNum"+"\t"+str(cNonMapNum)+"\n"+
          dResName+"\t"+"dNonRightMapNum"+"\t"+str(dNonRightMapNum)+"\n"+
          dResName+"\t"+"dErrorMapNum"+"\t"+str(dErrorMapNum)+"\n"+
          dResName+"\t"+"dMultMapNum"+"\t"+str(dmultMapNum)+"\n"+
          dResName+"\t"+"dNonMapNum"+"\t"+str(dNonMapNum)+"\n"+
          eResName+"\t"+"eNonRightMapNum"+"\t"+str(eNonRightMapNum)+"\n"+
          eResName+"\t"+"eErrorMapNum"+"\t"+str(eErrorMapNum)+"\n"+
          eResName+"\t"+"eMultMapNum"+"\t"+str(emultMapNum)+"\n"+
          eResName+"\t"+"eNonMapNum"+"\t"+str(eNonMapNum)+"\n")
          
	
