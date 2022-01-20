# -*- coding: utf-8 -*-

import argparse
import os
from subprocess import Popen,PIPE
from operator import itemgetter
import gzip
parser = argparse.ArgumentParser(description='The analysis of unsatisfactory aligned reads')
parser.add_argument('-a', '--afile', help='Input a file,end by .sam or .bam', required=True)
parser.add_argument('-b', '--bfile', help='Input b file,end by .sam or .bam', required=True)
parser.add_argument('-c', '--cfile', help='Input c file,end by .sam or .bam', required=True)
parser.add_argument('-d', '--dfile', help='Input d file,end by .sam or .bam', required=True)
parser.add_argument('-an', '--afileName', help='Input a name of afile result,such as bismarkbwt2', required=True)
parser.add_argument('-bn', '--bfileName', help='Input a name of bfile result,such as bsmap', required=True)
parser.add_argument('-cn', '--cfileName', help='Input a name of cfile result,such as bwameth', required=True)
parser.add_argument('-dn', '--dfileName', help='Input a name of dfile result,such as walt', required=True)
parser.add_argument('-bis', '--bismarkbwt2', help='Input a file of bismarkbwt2 ambiguouse,end by fq.gz', required=True)
parser.add_argument('-f', '--fastq', help='Input a file of fastq,end by .fq or .fastq', required=True)
parser.add_argument('-s', '--chromSize', help='Input a file of chrom size,end by .txt', required=True)
parser.add_argument('-e', '--errorReport', help='Statistical report of abnormal data.end by .txt', required=True)
parser.add_argument('-o', '--output', help='Folder to store output results,such as /path', required=True)
args = parser.parse_args()

## Parameter naming
aFile = args.afile
bFile = args.bfile
cFile = args.cfile
dFile = args.dfile
aResName= args.afileName
bResName= args.bfileName
cResName= args.cfileName
dResName= args.dfileName
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
aErrorMapReadId=tempFolder+"/"+aResName+"ErrorMapReadId.txt"
readIdRightError(aFile,aRightMapReadId,aErrorMapReadId,chromSizeFile)

bRightMapReadId=tempFolder+"/"+bResName+"RightMapReadId.txt"
bErrorMapReadId=tempFolder+"/"+bResName+"ErrorMapReadId.txt"
readIdRightError(bFile,bRightMapReadId,bErrorMapReadId,chromSizeFile)

cRightMapReadId=tempFolder+"/"+cResName+"RightMapReadId.txt"
cErrorMapReadId=tempFolder+"/"+cResName+"ErrorMapReadId.txt"
readIdRightError(cFile,cRightMapReadId,cErrorMapReadId,chromSizeFile)

dRightMapReadId=tempFolder+"/"+dResName+"RightMapReadId.txt"
dErrorMapReadId=tempFolder+"/"+dResName+"ErrorMapReadId.txt"
readIdRightError(dFile,dRightMapReadId,dErrorMapReadId,chromSizeFile)

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
					f.write(line)
		else:
			A_Bset=Aset.difference(Bset)
			with open(outFileName,"a") as f:
				for line in A_Bset:
					f.write(line)
	os.system("rm -r %s" % atempfilename)
	os.system("rm -r %s" % btempfilename)	
	
aNonRightMapReadId=tempFolder+"/"+aResName+"NonRightMapReadId.txt"
bNonRightMapReadId=tempFolder+"/"+bResName+"NonRightMapReadId.txt"
cNonRightMapReadId=tempFolder+"/"+cResName+"NonRightMapReadId.txt"
dNonRightMapReadId=tempFolder+"/"+dResName+"NonRightMapReadId.txt"

diff_seta_setb(allReadIdFile,aRightMapReadId,tempFolder,aNonRightMapReadId)	
diff_seta_setb(allReadIdFile,bRightMapReadId,tempFolder,bNonRightMapReadId)
diff_seta_setb(allReadIdFile,cRightMapReadId,tempFolder,cNonRightMapReadId)
diff_seta_setb(allReadIdFile,dRightMapReadId,tempFolder,dNonRightMapReadId)

## Incorrect unique mapped readID of common and difference
def intersecStaResOut(afile,bfile,cfile,dfile,tempPath,aDiffResFileName,bDiffResFileName,cDiffResFileName,dDiffResFileName,commonResFileName,reportResFileName):
	size=1000000
	Anum=0
	Bnum=0
	Cnum=0
	Dnum=0
	ABnum=0
	ACnum=0
	ADnum=0
	BCnum=0
	BDnum=0
	CDnum=0
	ABCnum=0
	ABDnum=0
	ACDnum=0
	BCDnum=0
	ABCDnum=0
	A_ABCDnum=0
	B_ABCDnum=0
	C_ABCDnum=0
	D_ABCDnum=0
	forrunnum=0
	filelist=[]
	mapperlist=["a","b","c","d"]
	aDiffResFile = open(aDiffResFileName,"w")
	bDiffResFile = open(bDiffResFileName,"w")
	cDiffResFile = open(cDiffResFileName,"w")
	dDiffResFile = open(dDiffResFileName,"w")
	commonResFile = open(commonResFileName,"w")
	reportResFile = open(reportResFileName,"w")
	
	atempfilename=''.join((tempPath,"/atemp/"))
	if not os.path.exists(atempfilename):
		os.mkdir(atempfilename)
	cutFiles(afile,size,filelist,atempfilename)

	btempfilename=''.join((tempPath,"/btemp/"))
	if not os.path.exists(btempfilename):
		os.mkdir(btempfilename)
	cutFiles(bfile,size,filelist,btempfilename)	
					
	ctempfilename=''.join((tempPath,"/ctemp/"))
	if not os.path.exists(ctempfilename):
		os.mkdir(ctempfilename)
	cutFiles(cfile,size,filelist,ctempfilename)	

	dtempfilename=''.join((tempPath,"/dtemp/"))
	if not os.path.exists(dtempfilename):
		os.mkdir(dtempfilename)
	cutFiles(dfile,size,filelist,dtempfilename)
	
	for fileid in filelist:
		forrunnum += 1
		Aset=set()
		Bset=set()
		Cset=set()
		Dset=set()
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
		ABset=Aset.intersection(Bset)
		ABnum += len(ABset)
			
		ACset=Aset.intersection(Cset)
		ACnum += len(ACset)
			
		ADset=Aset.intersection(Dset)
		ADnum += len(ADset)
			
		BCset=Bset.intersection(Cset)
		BCnum += len(BCset)
		
		BDset=Bset.intersection(Dset)
		BDnum += len(BDset)
		
		CDset=Cset.intersection(Dset)
		CDnum += len(CDset)
		
		ABCset=Aset.intersection(BCset)
		ABCnum += len(ABCset)
		
		ABDset=Aset.intersection(BDset)
		ABDnum += len(ABDset)
		
		ACDset=Aset.intersection(CDset)
		ACDnum += len(ACDset)
		
		BCDset=Bset.intersection(CDset)
		BCDnum += len(BCDset)
		
		ABCDset = Aset.intersection(BCDset)
		for line in ABCDset:
			ABCDnum += 1
			commonResFile.write(line)
		
		A_ABCDset = Aset.difference(ABCDset)
		for line in A_ABCDset:
			A_ABCDnum += 1
			aDiffResFile.write(line)
		
		B_ABCDset = Bset.difference(ABCDset)
		for line in B_ABCDset:
			B_ABCDnum += 1
			bDiffResFile.write(line)
		
		C_ABCDset = Cset.difference(ABCDset)
		for line in C_ABCDset:
			C_ABCDnum += 1
			cDiffResFile.write(line)
		
		D_ABCDset = Dset.difference(ABCDset)
		for line in D_ABCDset:
			D_ABCDnum += 1
			dDiffResFile.write(line)
	
	unionABCDnum=Anum+Bnum-ABnum+Cnum-ACnum-BCnum+ABCnum+Dnum-ADnum-BDnum+ABDnum-CDnum+ACDnum+BCDnum-ABCDnum
	reportResFile.write("A\t"+aDiffResFileName+"\n"+
											"B\t"+bDiffResFileName+"\n"+
											"C\t"+cDiffResFileName+"\n"+
											"D\t"+dDiffResFileName+"\n"+
											"Anum\t"+str(Anum)+"\n"+
											"Bnum\t"+str(Bnum)+"\n"+
											"Cnum\t"+str(Cnum)+"\n"+
											"Dnum\t"+str(Dnum)+"\n"+
											"ABnum\t"+str(ABnum)+"\n"+
											"ACnum\t"+str(ACnum)+"\n"+
											"ADnum\t"+str(ADnum)+"\n"+
											"BCnum\t"+str(BCnum)+"\n"+
											"BDnum\t"+str(BDnum)+"\n"+
											"CDnum\t"+str(CDnum)+"\n"+
											"ABCnum\t"+str(ABCnum)+"\n"+
											"ABDnum\t"+str(ABDnum)+"\n"+
											"ACDnum\t"+str(ACDnum)+"\n"+
											"BCDnum\t"+str(BCDnum)+"\n"+
											"ABCDnum\t"+str(ABCDnum)+"\n"+
											"A_ABCDnum\t"+str(A_ABCDnum)+"\n"+
											"B_ABCDnum\t"+str(B_ABCDnum)+"\n"+
											"C_ABCDnum\t"+str(C_ABCDnum)+"\n"+
											"D_ABCDnum\t"+str(D_ABCDnum)+"\n"+
											"unionSetOfANBCD\t"+str(unionABCDnum)+"\n")
	os.system("rm -r %s" %atempfilename)
	os.system("rm -r %s" %btempfilename)
	os.system("rm -r %s" %ctempfilename)
	os.system("rm -r %s" %dtempfilename)
	return(ABCDnum)

aDiffNonRightMapReadId=readIdFolder+"/"+aResName+"DiffNonRightMapReadId.txt"
bDiffNonRightMapReadId=readIdFolder+"/"+bResName+"DiffNonRightMapReadId.txt"
cDiffNonRightMapReadId=readIdFolder+"/"+cResName+"DiffNonRightMapReadId.txt"
dDiffNonRightMapReadId=readIdFolder+"/"+dResName+"DiffNonRightMapReadId.txt" 
commonNonRightMapReadId=readIdFolder+"/commonNonRightMapReadId.txt"
nonRightMapReadIdReport=readIdFolder+"/nonRightMapReadIdReport.txt"
commonNonRightMapReadIdNum=intersecStaResOut(afile=aNonRightMapReadId,\
									bfile=bNonRightMapReadId,\
									cfile=cNonRightMapReadId,\
									dfile=dNonRightMapReadId,\
									tempPath=tempFolder,\
									aDiffResFileName=aDiffNonRightMapReadId,\
									bDiffResFileName=bDiffNonRightMapReadId,\
									cDiffResFileName=cDiffNonRightMapReadId,\
									dDiffResFileName=dDiffNonRightMapReadId,\
									commonResFileName=commonNonRightMapReadId,\
									reportResFileName=nonRightMapReadIdReport)
									
## mapping class of nonRight mapped reads
def read1And2Match(readIdDict,outFileName,num):
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

def multMapReadIdExtract(bamOrSamFile,outFileName,chromSizeFileName):
    com="samtools view %s|awk 'END{print NR}'" % bamOrSamFile
    p=Popen(com,stdout=PIPE,shell=True,universal_newlines=True)
    maxLineNum=int(p.stdout.read().strip())
    samFile=samOrBam(bamOrSamFile)
    readId=''
    lineNum=0
    num=0
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
                    read1And2Match(readIdDict,outFileName,num)
                else:
                    if lineFlags in readIdDict.keys():
                        lineReadId = line[0]
                        if readId=='':
                            readId=lineReadId
                            readIdDict[lineFlags].append(line)
                        else:
                            if lineReadId != readId:
                                num += 1
                                read1And2Match(readIdDict,outFileName,num)
                                readId=lineReadId
                                readIdDict={'339':[],'355':[],'419':[],'403':[]}
                                readIdDict[lineFlags].append(line)
                            elif lineReadId == readId:
                                readIdDict[lineFlags].append(line)
aMultMapReadId=tempFolder+"/"+aResName+"multMapReadId.txt"
bMultMapReadId=tempFolder+"/"+bResName+"multMapReadId.txt"
cMultMapReadId=tempFolder+"/"+cResName+"multMapReadId.txt"
dMultMapReadId=tempFolder+"/"+dResName+"multMapReadId.txt"

def extractReadID2(inFileName,outFileName,chromSizeFileName):
    chromSizeDict=creatChromSizeDict(chromSizeFileName)
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
                r1ReadID=line+"_R1\n"
                r2ReadID=line+"_R2\n"
                outFile.write(r1ReadID+r2ReadID)
        else:
            continue

extractReadID2(bismarkbwt2Mul,aMultMapReadId,chromSizeFile)
multMapReadIdExtract(bFile,bMultMapReadId,chromSizeFile)
multMapReadIdExtract(cFile,cMultMapReadId,chromSizeFile)
multMapReadIdExtract(dFile,dMultMapReadId,chromSizeFile)

def nonRightMapReadMappingClass(errorMapReadIdFileName,multMapReadIdFileName,diffNonRightMapReadIdFileName,commonNonRightMapReadIdNum,outFileName,tempPath,resName):
	size=1000000
	infilelist=["error","mult","diff"]
	fileList=[]
	errorNum=0
	multNum=0
	diffNum=0
	diffErrorNum=0
	diffMultNum=0
	
	errorTempFileName=''.join((tempPath,"/errorTemp/"))
	if not os.path.exists(errorTempFileName):
		os.mkdir(errorTempFileName)
	fileList=cutFiles(errorMapReadIdFileName,size,fileList,errorTempFileName)
	
	if os.stat(multMapReadIdFileName).st_size > 0:
		multTempFileName=''.join((tempPath,"/multTemp/"))
		if not os.path.exists(multTempFileName):
			os.mkdir(multTempFileName)
		fileList=cutFiles(multMapReadIdFileName,size,fileList,multTempFileName)
	
	diffTempFileName=''.join((tempPath,"/diffTemp/"))
	if not os.path.exists(diffTempFileName):
		os.mkdir(diffTempFileName)
	fileList=cutFiles(diffNonRightMapReadIdFileName,size,fileList,diffTempFileName)
	
	for fileid in fileList:
		errorSet=set()
		multSet=set()
		diffSet=set()
		for infile in infilelist:
			if infile=="error":
				filename=''.join((errorTempFileName,"temp",str(fileid),".txt"))
				if os.path.exists(filename):
					with open(filename,"r") as f:
						for line in f:
							errorNum += 1
							errorSet.add(line)
			if os.stat(multMapReadIdFileName).st_size > 0:
				if infile=="mult":
					filename=''.join((multTempFileName,"temp",str(fileid),".txt"))
					if os.path.exists(filename):
						with open(filename,"r") as f:
							for line in f:
								multNum += 1
								multSet.add(line)
			if infile=="diff":
				filename=''.join((diffTempFileName,"temp",str(fileid),".txt"))
				if os.path.exists(filename):
					with open(filename,"r") as f:
						for line in f:
							diffNum += 1
							diffSet.add(line)
		
		diffErrorNum += len(diffSet.intersection(errorSet))
		if os.path.exists(multMapReadIdFileName):
			diffMultNum += len(diffSet.intersection(multSet))
		else:
			diffMultNum = 0
	if os.stat(multMapReadIdFileName).st_size > 0:
		diffNonNum = diffNum-diffErrorNum - diffMultNum
		commonErrorNum = errorNum - diffErrorNum
		commonMultNum = multNum -diffMultNum
		commonNonNum = commonNonRightMapReadIdNum - commonErrorNum - commonMultNum
		nonNum = diffNonNum + commonNonNum
	else:
		diffNonNum = diffNum-diffErrorNum - diffMultNum
		commonErrorNum = errorNum - diffErrorNum
		commonMultNum=0
		multNum=0
		commonNonNum = commonNonRightMapReadIdNum - commonErrorNum - commonMultNum
		nonNum = diffNonNum + commonNonNum
	with open(outFileName,"w") as df:
		df.write(resName+"\t"+"allErrorMap"+"\t"+str(errorNum)+"\n"+
						resName+"\t"+"allMultMap"+"\t"+str(multNum)+"\n"+
						resName+"\t"+"allNonMap"+"\t"+str(nonNum)+"\n"+
						resName+"\t"+"commonErrorMap"+"\t"+str(commonErrorNum)+"\n"+
						resName+"\t"+"commonMultMap"+"\t"+str(commonMultNum)+"\n"+
						resName+"\t"+"commonNonMap"+"\t"+str(commonNonNum)+"\n"+
						resName+"\t"+"diffErrorMap"+"\t"+str(diffErrorNum)+"\n"+
						resName+"\t"+"diffMultMap"+"\t"+str(diffMultNum)+"\n"+
						resName+"\t"+"diffNonMap"+"\t"+str(diffNonNum)+"\n")
	os.system("rm -r %s" % errorTempFileName)
	if os.stat(multMapReadIdFileName).st_size > 0:
		os.system("rm -r %s" % multTempFileName)
	os.system("rm -r %s" % diffTempFileName)

aNonRightMapReadMappingClassReport=readIdFolder+"/"+aResName+"NonRightMapReadMappingClassReport.txt"
bNonRightMapReadMappingClassReport=readIdFolder+"/"+bResName+"NonRightMapReadMappingClassReport.txt"
cNonRightMapReadMappingClassReport=readIdFolder+"/"+cResName+"NonRightMapReadMappingClassReport.txt"
dNonRightMapReadMappingClassReport=readIdFolder+"/"+dResName+"NonRightMapReadMappingClassReport.txt"

#commonNonRightMapReadIdNum=5332503

nonRightMapReadMappingClass(errorMapReadIdFileName=aErrorMapReadId,\
														multMapReadIdFileName=aMultMapReadId,\
														diffNonRightMapReadIdFileName=aDiffNonRightMapReadId,\
														commonNonRightMapReadIdNum=commonNonRightMapReadIdNum,\
														outFileName=aNonRightMapReadMappingClassReport,\
														tempPath=tempFolder,\
														resName=aResName)
nonRightMapReadMappingClass(errorMapReadIdFileName=bErrorMapReadId,\
														multMapReadIdFileName=bMultMapReadId,\
														diffNonRightMapReadIdFileName=bDiffNonRightMapReadId,\
														commonNonRightMapReadIdNum=commonNonRightMapReadIdNum,\
														outFileName=bNonRightMapReadMappingClassReport,\
														tempPath=tempFolder,\
														resName=bResName)
nonRightMapReadMappingClass(errorMapReadIdFileName=cErrorMapReadId,\
														multMapReadIdFileName=cMultMapReadId,\
														diffNonRightMapReadIdFileName=cDiffNonRightMapReadId,\
														commonNonRightMapReadIdNum=commonNonRightMapReadIdNum,\
														outFileName=cNonRightMapReadMappingClassReport,\
														tempPath=tempFolder,\
														resName=cResName)
nonRightMapReadMappingClass(errorMapReadIdFileName=dErrorMapReadId,\
														multMapReadIdFileName=dMultMapReadId,\
														diffNonRightMapReadIdFileName=dDiffNonRightMapReadId,\
														commonNonRightMapReadIdNum=commonNonRightMapReadIdNum,\
														outFileName=dNonRightMapReadMappingClassReport,\
														tempPath=tempFolder,\
														resName=dResName)
