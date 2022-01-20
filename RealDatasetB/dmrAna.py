

import argparse
import os
parser = argparse.ArgumentParser(description='The analysis of concordant and discordant DMR')
parser.add_argument('-a', '--afile', help='The file of DMR obtained by DSS,such as bismarbwt2_DMR.txt', required=True)
parser.add_argument('-b', '--bfile', help='The file of DMR obtained by DSS,such as bsmap_DMR.txt', required=True)
parser.add_argument('-c', '--cfile', help='The file of DMR obtained by DSS,such as bwameth_DMR.txt', required=True)
parser.add_argument('-d', '--dfile', help='The file of DMR obtained by DSS,such as walt_DMR.txt', required=True)
parser.add_argument('-r', '--repanno', help='Input annotation file of repeat sequence,such as ${species}rmskchr.bed', required=True)
parser.add_argument('-cg', '--cgianno', help='Input annotation file of CGI related area,such as CGIandNonCGI.bed', required=True)
parser.add_argument('-o', '--output', help='Path and the prefix of the file name,such as path+species', required=True)
args = parser.parse_args()

aFile = args.afile
bFile = args.bfile
cFile = args.cfile
dFile = args.dfile
repAnnoFile = args.repanno
cgiAnnoFile = args.cgianno
outPut =args.output

### 
aFileDmrBedSort=os.path.splitext(aFile)[0]+"_sort.bed"
bFileDmrBedSort=os.path.splitext(bFile)[0]+"_sort.bed"
cFileDmrBedSort=os.path.splitext(cFile)[0]+"_sort.bed"
dFileDmrBedSort=os.path.splitext(dFile)[0]+"_sort.bed"

def txtToBedOfDmr(inputFileName,outputFileName):
    inputFile=open(inputFileName,"r")
    dmrBedFileName=os.path.splitext(inputFileName)[0]+".bed"
    dmrBedFile=open(dmrBedFileName,"w")
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
            dmrBedFile.write(a+"\n")
    dmrBedFile.close()
    os.system("sort -t $'\t' -k 1,1 -k 2,2n -k 3,3n %s > %s" %(dmrBedFileName,outputFileName))
    os.system("rm -r %s" % dmrBedFileName)
		
txtToBedOfDmr(aFile,aFileDmrBedSort)
txtToBedOfDmr(bFile,bFileDmrBedSort)
txtToBedOfDmr(cFile,cFileDmrBedSort)
txtToBedOfDmr(dFile,dFileDmrBedSort)


### 一致性和非一致性DMR的划分
aMapperName=os.path.basename(aFile).split('_')[0]
bMapperName=os.path.basename(bFile).split('_')[0]
cMapperName=os.path.basename(cFile).split('_')[0]
dMapperName=os.path.basename(dFile).split('_')[0]

conAndNonconsisFile=outPut+"_conAndNonconsis.txt"

os.system("bedtools multiinter -header -i %s %s %s %s -names %s %s %s %s > %s" %(aFileDmrBedSort,bFileDmrBedSort,cFileDmrBedSort,dFileDmrBedSort,aMapperName,bMapperName,cMapperName,dMapperName,conAndNonconsisFile))

vennResFileName=outPut+"_vennRes.txt"
concordantFileName=outPut+"_concordant.bed"
discordantFileName=outPut+"_discordant.bed"

def resSta(inputFileName,outputFileName,concordantFileNamqe,discordantFileName):
    inputFile=open(inputFileName,"r")
    outputFile=open(outputFileName,"w")
    concordant=open(concordantFileName,"w")
    discordant=open(discordantFileName,"w")
    
    resDict={}
    for line in inputFile:
        if "start" in line:
            continue
        else:
            line=line.split()
            length=int(line[2])-int(line[1])
            if line[4] in resDict:
                resDict[line[4]] += length
            else:
                resDict[line[4]] = length
            
            if line[3]=="4":   # 一致性DMR区域
                concordant.write("\t".join(line[0:3])+"\n")
            else:              # 非一致性DMR区域
                discordant.write("\t".join(line[0:3])+"\n")       
    
    def oneMap(resDict,resClass):
        a=0
        for key in resDict.keys():
            if resClass in key:
                a += resDict[key]
        return(a)
    a = oneMap(resDict,'bismarkbwt2')
    b = oneMap(resDict,'bsmap')
    c = oneMap(resDict,'bwameth')
    d = oneMap(resDict,'walt')
    
    def twoMap(resDict,resClass1,resClass2):
        a=0
        for key in resDict.keys():
            if resClass1 in key and resClass2 in key:
                a += resDict[key]
        return(a)
    
    ab=twoMap(resDict,'bismarkbwt2','bsmap')
    ac=twoMap(resDict,'bismarkbwt2','bwameth')
    ad=twoMap(resDict,'bismarkbwt2','walt')
    bc=twoMap(resDict,'bsmap','bwameth')
    bd=twoMap(resDict,'bsmap','walt')
    cd=twoMap(resDict,'bwameth','walt')
    
    def threeMap(resDict,resClass1,resClass2,resClass3):
        a=0
        for key in resDict.keys():
            if resClass1 in key and resClass2 in key and resClass3 in key:
                a += resDict[key]
        return(a)
    abc=threeMap(resDict,'bismarkbwt2','bsmap','bwameth')
    abd=threeMap(resDict,'bismarkbwt2','bsmap','walt')
    acd=threeMap(resDict,'bismarkbwt2','bwameth','walt')
    bcd=threeMap(resDict,'bsmap','bwameth','walt')
    
    def fourMap(resDict,resClass1,resClass2,resClass3,resClass4):
        a=0
        for key in resDict.keys():
            if resClass1 in key and resClass2 in key and resClass3 in key and resClass4 in key:
                a += resDict[key]
        return(a)
    abcd=fourMap(resDict,'bismarkbwt2','bsmap','bwameth','walt')
    
    def allLenSta(resDict):
        a=0
        for key in resDict.keys():
            a += resDict[key]
        return(a)
    allLen=allLenSta(resDict)
    
    outputFile.write("aMapperName"+"\t"+"Bismark-bwt2"+"\n"+
                     "bMapperName"+"\t"+"BSMAP"+"\n"+
                     "cMapperName"+"\t"+"Bwa-meth"+"\n"+
                     "dMapperName"+"\t"+"Walt"+"\n"+
                     "a"+"\t"+str(a)+"\n"+
                     "b"+"\t"+str(b)+"\n"+
                     "c"+"\t"+str(c)+"\n"+
                     "d"+"\t"+str(d)+"\n"+
                     "ab"+"\t"+str(ab)+"\n"+
                     "ac"+"\t"+str(ac)+"\n"+
                     "ad"+"\t"+str(ad)+"\n"+
                     "bc"+"\t"+str(bc)+"\n"+
                     "bd"+"\t"+str(bd)+"\n"+
                     "cd"+"\t"+str(cd)+"\n"+
                     "abc"+"\t"+str(abc)+"\n"+
                     "abd"+"\t"+str(abd)+"\n"+
                     "acd"+"\t"+str(acd)+"\n"+
                     "bcd"+"\t"+str(bcd)+"\n"+
                     "abcd"+"\t"+str(abcd)+"\n"+
                     "allLen"+"\t"+str(allLen)+"\n")
resSta(conAndNonconsisFile,vennResFileName,concordantFileName,discordantFileName)                     

discordantFileName1=outPut+"_discordant1.bed"
os.system("bedtools merge -i %s > %s" %(discordantFileName,discordantFileName1))
os.system("mv %s %s" %(discordantFileName1,discordantFileName))

concorAndRepResFile=outPut+"_concordant_Rep_res.txt"
discorAndRepResFile=outPut+"_discordant_Rep_res.txt"
concorAndCgiResFile=outPut+"_concordant_Cgi_res.txt"
discorAndCgiResFile=outPut+"_discordant_Cgi_res.txt"

speciesName=outPut.split("/")[-1]

os.system("bedtools multiinter -header -i %s %s -names %s.concordant repeat > %s" %(concordantFileName,repAnnoFile,speciesName,concorAndRepResFile))
os.system("bedtools multiinter -header -i %s %s -names %s.discordant repeat > %s" %(discordantFileName,repAnnoFile,speciesName,discorAndRepResFile))
os.system("bedtools multiinter -header -i %s %s -names %s.concordant CGI > %s" %(concordantFileName,cgiAnnoFile,speciesName,concorAndCgiResFile))
os.system("bedtools multiinter -header -i %s %s -names %s.discordant CGI > %s" %(discordantFileName,cgiAnnoFile,speciesName,discorAndCgiResFile))


def singleFile(inputFileName,outputFile):
    inputFile=open(inputFileName,"r")

    resDict={}
    for line in inputFile:
        line=line.split()
        if "start" in line:
            sampleAndType=line[5]
            region=line[6]
            sampleRegion=",".join((sampleAndType,region))
            resDict[sampleAndType] = 0
            resDict[sampleRegion] = 0
            continue
        else:
            length=int(line[2])-int(line[1])
            if line[4] == sampleRegion:
                resDict[sampleRegion] += length
            elif line[4] == sampleAndType:
                resDict[sampleAndType] += length
    sample=sampleAndType.split('.')[0]
    Type=sampleAndType.split('.')[1]
    nonRegion="non_"+region
    outputFile.write(sample+"\t"+region+"\t"+Type+"\t"+str(resDict[sampleRegion])+"\n")
    outputFile.write(sample+"\t"+nonRegion+"\t"+Type+"\t"+str(resDict[sampleAndType])+"\n")

def resSta(inputFileName1,inputFileName2,inputFileName3,inputFileName4,outputFileName):
    outputFile=open(outputFileName,"w")
    outputFile.write("species"+"\t"+"region"+"\t"+"type"+"\t"+"length"+"\n")
    singleFile(inputFileName1,outputFile)
    singleFile(inputFileName2,outputFile)
    singleFile(inputFileName3,outputFile)
    singleFile(inputFileName4,outputFile)

repAndCgiReport=outPut+"_rep_cgi_report.txt" 
resSta(concorAndRepResFile,discorAndRepResFile,concorAndCgiResFile,discorAndCgiResFile,repAndCgiReport)

os.system("rm %s" % concordantFileName)
os.system("rm %s" % discordantFileName)
os.system("rm %s" % concorAndRepResFile)
os.system("rm %s" % discorAndRepResFile)
os.system("rm %s" % concorAndCgiResFile)
os.system("rm %s" % discorAndCgiResFile)
