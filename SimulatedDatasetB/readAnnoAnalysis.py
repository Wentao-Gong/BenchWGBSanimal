import os
import argparse
from subprocess import Popen,PIPE

parser = argparse.ArgumentParser(description='The influence of repetitive sequence and CGI on the unsatisfactory aligned reads')
parser.add_argument('-a', '--afile', help='Input a file,end by .txt', required=True)
parser.add_argument('-b', '--bfile', help='Input b file,end by .txt', required=True)
parser.add_argument('-c', '--cfile', help='Input c file,end by .txt', required=True)
parser.add_argument('-d', '--dfile', help='Input d file,end by .txt', required=True)
parser.add_argument('-e', '--efile', help='Input e file,end by .txt', required=True)
parser.add_argument('-f', '--ffile', help='Input f file,end by .txt', required=True)
parser.add_argument('-r', '--repanno', help='Input annotation file of repeat sequence', required=True)
parser.add_argument('-cg', '--cgianno', help='Input annotation file of CGI region', required=True)
parser.add_argument('-o', '--outpath', help='Folder to store output results,such as /path', required=True)
args = parser.parse_args()

aFile = args.afile
bFile = args.bfile
cFile = args.cfile
dFile = args.dfile
eFile = args.efile
fFile = args.ffile
repAnnoFile = args.repanno
cgiAnnoFile = args.cgianno
outPath = args.outpath

# Determine whether the path exists, if it does not exist, build it
if not os.path.exists(outPath):
	os.system("mkdir -p %s" %outPath)

# Convert the readid file to the bed file format
def readidToBed(inFileName,outFileName):
    inFile=open(inFileName,"r")
    outFile=open(outFileName,"w")
    for line in inFile:
        chrom = line.split(':')[0].split('_',1)[1]
        if line.split('\t')[0].split(':')[1].split('-')[0] == "":  ## start is negative
            continue
        else:
            if line.split('_')[-1]=="R1\n":
                start = str(int(line.split('\t')[0].split(':')[1].split('-')[0])-1)
                end = str(int(line.split('\t')[0].split(':')[1].split('-')[0])+149)
            elif line.split('_')[-1]=="R2\n":
                end=line.split('\t')[0].split(':')[1].split('-')[1].split('_')[0]
                start=str(int(end)-150)
        outFile.write(chrom+"\t"+start+"\t"+end+"\t"+line)
            


aFileBed=outPath+"/"+os.path.basename(aFile)[:-4]+".bed"
bFileBed=outPath+"/"+os.path.basename(bFile)[:-4]+".bed"
cFileBed=outPath+"/"+os.path.basename(cFile)[:-4]+".bed"
dFileBed=outPath+"/"+os.path.basename(dFile)[:-4]+".bed"
eFileBed=outPath+"/"+os.path.basename(eFile)[:-4]+".bed"
fFileBed=outPath+"/"+os.path.basename(fFile)[:-4]+".bed"

readidToBed(aFile,aFileBed)
readidToBed(bFile,bFileBed)
readidToBed(cFile,cFileBed)
readidToBed(dFile,dFileBed)
readidToBed(eFile,eFileBed)
readidToBed(fFile,fFileBed)


def countLineNum(fileName):
    com="wc -l %s" % fileName
    p=Popen(com,stdout=PIPE,shell=True,universal_newlines=True)
    line_num=int(p.stdout.read().split()[0])
    return(line_num)
    
aFileBedLineNum=countLineNum(aFileBed)
bFileBedLineNum=countLineNum(bFileBed)
cFileBedLineNum=countLineNum(cFileBed)
dFileBedLineNum=countLineNum(dFileBed)
eFileBedLineNum=countLineNum(eFileBed)
fFileBedLineNum=countLineNum(fFileBed)

### repeat annotation analysis
# File comment
aAnnoResultFileName=outPath+"/"+os.path.basename(aFile)[:-4]+"RepeatAnnoRes.bed"
bAnnoResultFileName=outPath+"/"+os.path.basename(bFile)[:-4]+"RepeatAnnoRes.bed"
cAnnoResultFileName=outPath+"/"+os.path.basename(cFile)[:-4]+"RepeatAnnoRes.bed"
dAnnoResultFileName=outPath+"/"+os.path.basename(dFile)[:-4]+"RepeatAnnoRes.bed"
eAnnoResultFileName=outPath+"/"+os.path.basename(eFile)[:-4]+"RepeatAnnoRes.bed"
fAnnoResultFileName=outPath+"/"+os.path.basename(fFile)[:-4]+"RepeatAnnoRes.bed"

os.system('bedtools intersect -a %s -b %s -wao > %s' %(aFileBed,repAnnoFile,aAnnoResultFileName))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(bFileBed,repAnnoFile,bAnnoResultFileName))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(cFileBed,repAnnoFile,cAnnoResultFileName))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(dFileBed,repAnnoFile,dAnnoResultFileName))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(eFileBed,repAnnoFile,eAnnoResultFileName))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(fFileBed,repAnnoFile,fAnnoResultFileName))

# Annotation file result statistics
def repAnnoRes(annoResultFileName,tempFileName):
	annoResultFile=open(annoResultFileName,"r")
	tempFile=open(tempFileName,"w")
	for line in annoResultFile:
		if line.split()[13]=="0":
			tempFile.write(line.strip()+"\t"+"0"+"\t"+"0"+"\n")
		else:
			tempFile.write(line.strip()+"\t"+str(int(line.split()[13])/150)+"\t"+str(int(line.split()[13])/int(line.split()[12]))+"\n")

tempFileName=outPath+"/"+"temp.bed"
repAnnoRes(aAnnoResultFileName,tempFileName)
os.system('rm %s' %aAnnoResultFileName)
os.system('mv %s %s' %(tempFileName,aAnnoResultFileName))

repAnnoRes(bAnnoResultFileName,tempFileName)
os.system('rm %s' %bAnnoResultFileName)
os.system('mv %s %s' %(tempFileName,bAnnoResultFileName))

repAnnoRes(cAnnoResultFileName,tempFileName)
os.system('rm %s' %cAnnoResultFileName)
os.system('mv %s %s' %(tempFileName,cAnnoResultFileName))

repAnnoRes(dAnnoResultFileName,tempFileName)
os.system('rm %s' %dAnnoResultFileName)
os.system('mv %s %s' %(tempFileName,dAnnoResultFileName))

repAnnoRes(eAnnoResultFileName,tempFileName)
os.system('rm %s' %eAnnoResultFileName)
os.system('mv %s %s' %(tempFileName,eAnnoResultFileName))

repAnnoRes(fAnnoResultFileName,tempFileName)
os.system('rm %s' %fAnnoResultFileName)
os.system('mv %s %s' %(tempFileName,fAnnoResultFileName))

def repAnnoResSta(annoResultFileName,outputFileName,lineNum):
    annoResultFile=open(annoResultFileName,"r")
    outputFile=open(outputFileName,"w")
    resDict = {}
    for line in annoResultFile:
        line=line.split()
        if line[5] == '-1':
            if "nonRepeat" in resDict:
                resDict["nonRepeat"] += 1
            else:
                resDict["nonRepeat"] = 1
        else:
            repClass = line[9]
            if repClass in resDict:
                resDict[repClass] += 1
                if line[14] == "1":
                    if repClass+"IncludeRead" in resDict:
                        resDict[repClass+"IncludeRead"] += 1
                    else:
                        resDict[repClass+"IncludeRead"] = 1
                elif line[15] == "1":
                    if repClass+"IncludeByRead" in resDict:
                        resDict[repClass+"IncludeByRead"] += 1
                    else:
                        resDict[repClass+"IncludeByRead"] = 1
                else:
                    if repClass+"overlapRead" in resDict:
                        resDict[repClass+"overlapRead"] += 1
                    else:
                        resDict[repClass+"overlapRead"] = 1
            else:
                resDict[repClass] = 1
                if line[14] == "1":
                    if repClass+"IncludeRead" in resDict:
                        resDict[repClass+"IncludeRead"] += 1
                    else:
                        resDict[repClass+"IncludeRead"] = 1
                elif line[15] == "1":
                    if repClass+"IncludeByRead" in resDict:
                        resDict[repClass+"IncludeByRead"] += 1
                    else:
                        resDict[repClass+"IncludeByRead"] = 1
                else:
                    if repClass+"overlapRead" in resDict:
                        resDict[repClass+"overlapRead"] += 1
                    else:
                        resDict[repClass+"overlapRead"] = 1
    repeatNum = lineNum - resDict["nonRepeat"]
    outputFile.write("repeat"+"\t"+str(repeatNum )+"\n")
    for key in resDict.keys():
        outputFile.write(key+"\t"+str(resDict[key])+"\n")
		
aAnnoResultStaReportFileName=outPath+"/"+os.path.basename(aFile)[:-4]+"RepeatAnnoResStaReport.txt"
bAnnoResultStaReportFileName=outPath+"/"+os.path.basename(bFile)[:-4]+"RepeatAnnoResStaReport.txt"
cAnnoResultStaReportFileName=outPath+"/"+os.path.basename(cFile)[:-4]+"RepeatAnnoResStaReport.txt"
dAnnoResultStaReportFileName=outPath+"/"+os.path.basename(dFile)[:-4]+"RepeatAnnoResStaReport.txt"
eAnnoResultStaReportFileName=outPath+"/"+os.path.basename(eFile)[:-4]+"RepeatAnnoResStaReport.txt"
fAnnoResultStaReportFileName=outPath+"/"+os.path.basename(fFile)[:-4]+"RepeatAnnoResStaReport.txt"
	
repAnnoResSta(aAnnoResultFileName,aAnnoResultStaReportFileName,aFileBedLineNum)
repAnnoResSta(bAnnoResultFileName,bAnnoResultStaReportFileName,bFileBedLineNum)
repAnnoResSta(cAnnoResultFileName,cAnnoResultStaReportFileName,cFileBedLineNum)
repAnnoResSta(dAnnoResultFileName,dAnnoResultStaReportFileName,dFileBedLineNum)
repAnnoResSta(eAnnoResultFileName,eAnnoResultStaReportFileName,eFileBedLineNum)
repAnnoResSta(fAnnoResultFileName,fAnnoResultStaReportFileName,fFileBedLineNum)

### CGI annotation analysis
aAnnoResultOfCGI=outPath+"/"+os.path.basename(aFile)[:-4]+"CGIAnnoRes.bed"
bAnnoResultOfCGI=outPath+"/"+os.path.basename(bFile)[:-4]+"CGIAnnoRes.bed"
cAnnoResultOfCGI=outPath+"/"+os.path.basename(cFile)[:-4]+"CGIAnnoRes.bed"
dAnnoResultOfCGI=outPath+"/"+os.path.basename(dFile)[:-4]+"CGIAnnoRes.bed"
eAnnoResultOfCGI=outPath+"/"+os.path.basename(eFile)[:-4]+"CGIAnnoRes.bed"
fAnnoResultOfCGI=outPath+"/"+os.path.basename(fFile)[:-4]+"CGIAnnoRes.bed"

os.system('bedtools intersect -a %s -b %s -wao > %s' %(aFileBed,cgiAnnoFile,aAnnoResultOfCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(bFileBed,cgiAnnoFile,bAnnoResultOfCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(cFileBed,cgiAnnoFile,cAnnoResultOfCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(dFileBed,cgiAnnoFile,dAnnoResultOfCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(eFileBed,cgiAnnoFile,eAnnoResultOfCGI))
os.system('bedtools intersect -a %s -b %s -wao > %s' %(fFileBed,cgiAnnoFile,fAnnoResultOfCGI))


def cgiAnnoResSta(annoResultFileName,outputFileName,fileLineNum):
    annoResultFile=open(annoResultFileName,"r")
    outputFile=open(outputFileName,"w")
    nonCgiNum = 0
    cgiIncludeNum = 0
    cgiOverlapNum = 0
    for line in annoResultFile:
        line=line.split()
        if line[5] == '-1':
            nonCgiNum+=1
        else:
            if line[7]==150:
                cgiIncludeNum += 1
            else:
                cgiOverlapNum += 1
    cgiNum = fileLineNum - nonCgiNum
    outputFile.write("nonCgiNum"+"\t"+str(nonCgiNum)+"\n"+
                     "cgiNum"+"\t"+str(cgiNum)+"\n"+
								"cgiIncludeNum"+"\t"+str(cgiIncludeNum)+"\n"+
								"cgiOverlapNum"+"\t"+str(cgiOverlapNum)+"\n")

aCgiAnnoResultStaReportFileName=outPath+"/"+os.path.basename(aFile)[:-4]+"CgiAnnoResStaReport.txt"
bCgiAnnoResultStaReportFileName=outPath+"/"+os.path.basename(bFile)[:-4]+"CgiAnnoResStaReport.txt"
cCgiAnnoResultStaReportFileName=outPath+"/"+os.path.basename(cFile)[:-4]+"CgiAnnoResStaReport.txt"
dCgiAnnoResultStaReportFileName=outPath+"/"+os.path.basename(dFile)[:-4]+"CgiAnnoResStaReport.txt"
eCgiAnnoResultStaReportFileName=outPath+"/"+os.path.basename(eFile)[:-4]+"CgiAnnoResStaReport.txt"
fCgiAnnoResultStaReportFileName=outPath+"/"+os.path.basename(fFile)[:-4]+"CgiAnnoResStaReport.txt"


cgiAnnoResSta(aAnnoResultOfCGI,aCgiAnnoResultStaReportFileName,aFileBedLineNum)
cgiAnnoResSta(bAnnoResultOfCGI,bCgiAnnoResultStaReportFileName,bFileBedLineNum)
cgiAnnoResSta(cAnnoResultOfCGI,cCgiAnnoResultStaReportFileName,cFileBedLineNum)
cgiAnnoResSta(dAnnoResultOfCGI,dCgiAnnoResultStaReportFileName,dFileBedLineNum)
cgiAnnoResSta(eAnnoResultOfCGI,eCgiAnnoResultStaReportFileName,eFileBedLineNum)
cgiAnnoResSta(fAnnoResultOfCGI,fCgiAnnoResultStaReportFileName,fFileBedLineNum)


