import argparse
from subprocess import Popen,PIPE
parser = argparse.ArgumentParser(description='Regional division of CpG island: CGI, CGI Shores, CGI shelves')
parser.add_argument('-i', '--input', help='CpG island annotation file,end by .txt', required=True)
parser.add_argument('-o', '--output', help='Annotation file of CpG island related area(CGI CGI Shores CGI shelves),end by bed', required=True)
parser.add_argument('-s', '--chromSize', help='The size of each chromosome,end by .txt', required=True)
args = parser.parse_args()

inFile = open(args.input,'r')
outFile = open(args.output,'w')
chromSize = args.chromSize

##Dictionary of chrom Size
def creatChromSizeDict(chromSize):
  chromSizeDict={}
  with open(chromSize,"r") as f:
    for line in f:
      chromSizeDict[line.split()[0]] = int(line.split()[1])
  return(chromSizeDict)
chromSizeDict=creatChromSizeDict(chromSize)

##Calculate the maximum number of rows in a file
def calMaxNumRow(argsInput):
  com="wc -l %s" % args.input
  p=Popen(com,stdout=PIPE,shell=True,universal_newlines=True)
  maxLineNum=int(p.stdout.read().split()[0])
  return(maxLineNum)
maxLineNum=calMaxNumRow(args.input)

##Division of the last region of chromosome
def lastArea(chromSizeDict,lineChrom,lineChromEnd,outFile):
  chromEnd=chromSizeDict[lineChrom]
  interOldEnd=chromEnd-lineChromEnd
  if interOldEnd != 0:
    if interOldEnd > 2000:
      if interOldEnd >4000:
        cpgShoresStart=lineChromEnd
        cpgShoresEnd=lineChromEnd+2000
        cpgShelvesStart=lineChromEnd+2000
        cpgShelvesEnd=lineChromEnd+4000
        othersStart=lineChromEnd+4000
        othersEnd=chromEnd
        outFile.write(\
          lineChrom+"\t"+str(cpgShoresStart)+"\t"+str(cpgShoresEnd)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd-cpgShoresStart)+"\n"+\
          lineChrom+"\t"+str(cpgShelvesStart)+"\t"+str(cpgShelvesEnd)+"\t"+"cpgShelves"+"\t"+str(cpgShelvesEnd-cpgShelvesStart)+"\n"+\
          lineChrom+"\t"+str(othersStart)+"\t"+str(othersEnd)+"\t"+"others"+"\t"+str(othersEnd-othersStart)+"\n")
      else:
        cpgShoresStart=lineChromEnd
        cpgShoresEnd=lineChromEnd+2000
        cpgShelvesStart=lineChromEnd+2000
        cpgShelvesEnd=chromEnd
        outFile.write(\
          lineChrom+"\t"+str(cpgShoresStart)+"\t"+str(cpgShoresEnd)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd-cpgShoresStart)+"\n"+\
          lineChrom+"\t"+str(cpgShelvesStart)+"\t"+str(cpgShelvesEnd)+"\t"+"cpgShelves"+"\t"+str(cpgShelvesEnd-cpgShelvesStart)+"\n")
    else:
      cpgShoresStart=lineChromEnd
      cpgShoresEnd=chromEnd
      outFile.write(lineChrom+"\t"+str(cpgShoresStart)+"\t"+str(cpgShoresEnd)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd-cpgShoresStart)+"\n")

##Division of chromosome intermediate region
def interArea(newlineChromStart,newlineChromEnd,oldlineChromEnd,outFile,lineChrom):
  interNewOldLen=newlineChromStart-oldlineChromEnd
  if interNewOldLen != 0:
    if interNewOldLen > 4000:
      if interNewOldLen > 8000:
        cpgShoresStart1=oldlineChromEnd
        cpgShoresEnd1=oldlineChromEnd+2000
          
        cpgShelvesStart1=oldlineChromEnd+2000
        cpgShelvesEnd1=oldlineChromEnd+4000
        
        othersStart=oldlineChromEnd+4000
        othersEnd=newlineChromStart-4000
        
        cpgShelvesStart2=newlineChromStart-4000
        cpgShelvesEnd2=newlineChromStart-2000
      
        cpgShoresStart2=newlineChromStart-2000
        cpgShoresEnd2=newlineChromStart
        outFile.write(\
          lineChrom+"\t"+str(cpgShoresStart1)+"\t"+str(cpgShoresEnd1)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd1-cpgShoresStart1)+"\n"+\
          lineChrom+"\t"+str(cpgShelvesStart1)+"\t"+str(cpgShelvesEnd1)+"\t"+"cpgShelves"+"\t"+str(cpgShelvesEnd1-cpgShelvesStart1)+"\n"+\
          lineChrom+"\t"+str(othersStart)+"\t"+str(othersEnd)+"\t"+"others"+"\t"+str(othersEnd-othersStart)+"\n"+\
          lineChrom+"\t"+str(cpgShelvesStart2)+"\t"+str(cpgShelvesEnd2)+"\t"+"cpgShelves"+"\t"+str(cpgShelvesEnd2-cpgShelvesStart2)+"\n"+\
          lineChrom+"\t"+str(cpgShoresStart2)+"\t"+str(cpgShoresEnd2)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd2-cpgShoresStart2)+"\n"+\
          lineChrom+"\t"+str(newlineChromStart)+"\t"+str(newlineChromEnd)+"\t"+"cpgIsland"+"\t"+str(newlineChromEnd-newlineChromStart)+"\n")
      else:
        cpgShoresStart1=oldlineChromEnd
        cpgShoresEnd1=oldlineChromEnd+2000
      
        cpgShelvesStart1=oldlineChromEnd+2000
        cpgShelvesEnd1=newlineChromStart-2000
      
        cpgShoresStart2=newlineChromStart-2000
        cpgShoresEnd2=newlineChromStart
        outFile.write(\
          lineChrom+"\t"+str(cpgShoresStart1)+"\t"+str(cpgShoresEnd1)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd1-cpgShoresStart1)+"\n"+\
          lineChrom+"\t"+str(cpgShelvesStart1)+"\t"+str(cpgShelvesEnd1)+"\t"+"cpgShelves"+"\t"+str(cpgShelvesEnd1-cpgShelvesStart1)+"\n"+\
          lineChrom+"\t"+str(cpgShoresStart2)+"\t"+str(cpgShoresEnd2)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd2-cpgShoresStart2)+"\n"+\
          lineChrom+"\t"+str(newlineChromStart)+"\t"+str(newlineChromEnd)+"\t"+"cpgIsland"+"\t"+str(newlineChromEnd-newlineChromStart)+"\n")
    else:
      cpgShoresStart=oldlineChromEnd
      cpgShoresEnd=newlineChromStart
      outFile.write(\
        lineChrom+"\t"+str(cpgShoresStart)+"\t"+str(cpgShoresEnd)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd-cpgShoresStart)+"\n"+\
        lineChrom+"\t"+str(newlineChromStart)+"\t"+str(newlineChromEnd)+"\t"+"cpgIsland"+"\t"+str(newlineChromEnd-newlineChromStart)+"\n")
  else:
    outFile.write(lineChrom+"\t"+str(newlineChromStart)+"\t"+str(newlineChromEnd)+"\t"+"cpgIsland"+"\t"+str(newlineChromEnd-newlineChromStart)+"\n")

##Division of chromosome starting region
def startArea(newlineChromStart,newlineChromEnd,outFile,lineChrom):
  if newlineChromStart != 0:
    if newlineChromStart > 2000:
      if newlineChromStart > 4000:
        cpgShoresStart=newlineChromStart-2000
        cpgShoresEnd=newlineChromStart
        cpgShelvesStart=newlineChromStart-4000
        cpgShelvesEnd=newlineChromStart-2000
        othersStart=0
        othersEnd=newlineChromStart-4000
        outFile.write(\
          lineChrom+"\t"+str(othersStart)+"\t"+str(othersEnd)+"\t"+"others"+"\t"+str(othersEnd-othersStart)+"\n"+\
          lineChrom+"\t"+str(cpgShelvesStart)+"\t"+str(cpgShelvesEnd)+"\t"+"cpgShelves"+"\t"+str(cpgShelvesEnd-cpgShelvesStart)+"\n"+\
          lineChrom+"\t"+str(cpgShoresStart)+"\t"+str(cpgShoresEnd)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd-cpgShoresStart)+"\n"+\
          lineChrom+"\t"+str(newlineChromStart)+"\t"+str(newlineChromEnd)+"\t"+"cpgIsland"+"\t"+str(newlineChromEnd-newlineChromStart)+"\n")
      else:
        cpgShoresStart=newlineChromStart-2000
        cpgShoresEnd=newlineChromStart
        cpgShelvesStart=0
        cpgShelvesEnd=newlineChromStart-2000
        outFile.write(\
          lineChrom+"\t"+str(cpgShelvesStart)+"\t"+str(cpgShelvesEnd)+"\t"+"cpgShelves"+"\t"+str(cpgShelvesEnd-cpgShelvesStart)+"\n"+\
          lineChrom+"\t"+str(cpgShoresStart)+"\t"+str(cpgShoresEnd)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd-cpgShoresStart)+"\n"+\
          lineChrom+"\t"+str(newlineChromStart)+"\t"+str(newlineChromEnd)+"\t"+"cpgIsland"+"\t"+str(newlineChromEnd-newlineChromStart)+"\n")
    else:
      cpgShoresStart=0
      cpgShoresEnd=newlineChromStart
      outFile.write(\
          lineChrom+"\t"+str(cpgShoresStart)+"\t"+str(cpgShoresEnd)+"\t"+"cpgShores"+"\t"+str(cpgShoresEnd-cpgShoresStart)+"\n"+\
          lineChrom+"\t"+str(newlineChromStart)+"\t"+str(newlineChromEnd)+"\t"+"cpgIsland"+"\t"+str(newlineChromEnd-newlineChromStart)+"\n")
  else:
    outFile.write(lineChrom+"\t"+str(newlineChromStart)+"\t"+str(newlineChromEnd)+"\t"+"cpgIsland"+"\t"+str(newlineChromEnd-newlineChromStart)+"\n")
##Chromosomes without CGI
def noCGIChrom(annoFileChrom,chromSizeDict,outFile):
	for chrom in chromSizeDict.keys():
		if chrom in annoFileChrom:
			continue
		else:
			othersStart=0
			othersEnd=chromSizeDict[chrom]
			outFile.write(chrom+"\t"+str(othersStart)+"\t"+str(othersEnd)+"\t"+"others"+"\t"+str(othersEnd-othersStart)+"\n")
 
##Regional division of CpG island
def regDivision(inFile,outFile,chromSizeDict,maxLineNum):
  chrom=""
  newline=""
  oldline=""
  lineNum=0
  annoFileChrom=set()
  for line in inFile:
    lineNum += 1
    oldline=newline
    if oldline != "":
      oldlineChromEnd=int(oldline[3])
    newline=line.split()
    newlineChrom=newline[1]
    annoFileChrom.add(newlineChrom)
    if int(newline[2]) != 0:
        newlineChromStart=int(newline[2])-1
    else:
        newlineChromStart=0
    newlineChromEnd=int(newline[3])
    ##last chrom
    if lineNum==maxLineNum:
      if newlineChrom==chrom:
        interArea(newlineChromStart=newlineChromStart,newlineChromEnd=newlineChromEnd,oldlineChromEnd=oldlineChromEnd,outFile=outFile,lineChrom=newlineChrom)
        lastArea(chromSizeDict=chromSizeDict,lineChrom=newlineChrom,lineChromEnd=newlineChromEnd,outFile=outFile)
      else:
        oldlineChrom=oldline[1]
        lastArea(chromSizeDict=chromSizeDict,lineChrom=oldlineChrom,lineChromEnd=oldlineChromEnd,outFile=outFile)
        startArea(newlineChromStart=newlineChromStart,newlineChromEnd=newlineChromEnd,outFile=outFile,lineChrom=newlineChrom)
        lastArea(chromSizeDict=chromSizeDict,lineChrom=newlineChrom,lineChromEnd=newlineChromEnd,outFile=outFile)
      noCGIChrom(annoFileChrom,chromSizeDict,outFile)	
    else:                                                   ##same chrom 
      if newlineChrom==chrom:
        interArea(newlineChromStart=newlineChromStart,newlineChromEnd=newlineChromEnd,oldlineChromEnd=oldlineChromEnd,outFile=outFile,lineChrom=newlineChrom)
      else:                                                ## different chrom
        chrom=newlineChrom
        if oldline=="":                                    ##first line and first chrom
          startArea(newlineChromStart=newlineChromStart,newlineChromEnd=newlineChromEnd,outFile=outFile,lineChrom=newlineChrom)
        else:                                              ##old chrom and new chrom
          oldlineChrom=oldline[1]
          lastArea(chromSizeDict=chromSizeDict,lineChrom=oldlineChrom,lineChromEnd=oldlineChromEnd,outFile=outFile)
          startArea(newlineChromStart=newlineChromStart,newlineChromEnd=newlineChromEnd,outFile=outFile,lineChrom=newlineChrom)

regDivision(inFile,outFile,chromSizeDict,maxLineNum)