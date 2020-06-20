#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import os
import sys
import linecache
import random
import numpy as np
import math

def sampCt(sampsCtF,minCtValue,maxCtValue):
	sampleCt29dic = {}
	sampsCtFls = open(sampsCtF).readlines()
	for Person_SIDFl in sampsCtFls:
		if "Patient" not in Person_SIDFl and Person_SIDFl != "\n":
			sample_SID = Person_SIDFl.split("\t")[0]
			Ct = Person_SIDFl.split("\t")[4].split("\n")[0]
			if Ct != 'NA' and float(Ct) >= float(minCtValue) and float(Ct) < float(maxCtValue) :
				sampleCt29dic[sample_SID] = Ct
	return sampleCt29dic

def outAllSampFStat(outAllSampF,sampleCt25dic,flag):	
	CtPosiSampDic = {}
	effDic = {}
	geneDic = {}
	outAllSampFls = open(outAllSampF).readlines()
	for outAllSampFl in outAllSampFls:
		if "sample" not in outAllSampFl and outAllSampFl != "\n":
			linekey = outAllSampFl.split("\n")[0].split("\t")
			sample = linekey[0].split("/")[-1]
			if sample in sampleCt25dic.keys():
				posi = linekey[1]
				Freq = linekey[17]
				Annotation = linekey[4]
				Gene_Name = linekey[5]

				if int(posi) < 29782 and int(posi) > 86 and int(posi) != 12436:
					if Freq != "NO" and Freq != "NA" and float(Freq) >= 5 and float(Freq) < 95:						
						if posi not in CtPosiSampDic:				
							CtPosiSampDic[posi] = [sample]
						else :					
							CtPosiSampDic[posi].append(sample)

				if sample + "_" +  posi not in effDic:
					effDic[sample + "_" + posi] = Annotation
				if posi  not in geneDic:
					geneDic[posi] = Annotation

	Ct25singlePosL = []
	shareNumDic = {}
	CtsinglePosD = CtPosiSampDic

	print "jj"
	print CtsinglePosD
	for Ct25Posi in CtPosiSampDic.keys():
		shareNumDic[Ct25Posi] = len(CtPosiSampDic[Ct25Posi])
		if flag == "single":			
			if len(CtPosiSampDic[Ct25Posi]) == 1:
				Ct25singlePosL.append(int(Ct25Posi))
			else:
				CtsinglePosD.pop(Ct25Posi)
		if flag == "shared":			
			if len(CtPosiSampDic[Ct25Posi]) > 1:
				Ct25singlePosL.append(int(Ct25Posi))
			else:
				CtsinglePosD.pop(Ct25Posi)


	Ct25singlePosL.sort()
	return Ct25singlePosL,CtsinglePosD,shareNumDic,effDic,geneDic,




def iSNVlieDic(iSNV):
	Head = linecache.getline(iSNV,3)
	HeadLst = Head.split("\t")
	lieDic = {}
	Count = 0
	for Headsamp in HeadLst:
		lieDic[Headsamp] = Count
		Count += 1
	return lieDic

def tongjisnpeff(sampsLst,iSNV,lieDic,effDic,shareNumDic,lowFreq,upFreq,singlePos):
	sysDic = {}
	missDic = {}
	nonsenseDic = {}
	for x in range(10,100,5):
		sysDic[x] = 0
		missDic[x] = 0

		
##统计
	iSNVls = open(iSNV).readlines()
	NumCount = 0
	for iSNVl in iSNVls:
		if "#" not in iSNVl  and iSNVl != "\n":
			lineLst = iSNVl.split("\t")
			posi = lineLst[0]
			if int(posi) in  singlePos: #				
				for sample in sampsLst:
					if sample in lieDic:						
						lieIndex = lieDic[sample]
						Freq = lineLst[lieIndex]					
						if Freq != "NO" and Freq != "NA"  and float(Freq) > float(lowFreq)  and float(Freq) <= float(upFreq) :
							if sample + "_" + posi in effDic:
								eff = effDic[sample + "_" + posi]
								NumCount += 1
	
	
							for x in range(10,100,5):		
								if float(Freq)*100  >= x-5 and float(Freq)*100 < x and eff == "synonymous_variant":
									sysDic[x] += 1						
								if float(Freq)*100 >= x-5  and float(Freq)*100 < x and (eff == "missense_variant" or eff == "stop_gained"):
									missDic[x] += 1


	synSum = 0
	print sysDic
	for KEY in sysDic.keys():
		synSum += sysDic[KEY]
		print "syn" + "\t" +str(KEY)  + "\t" + str(sysDic[KEY])

	missSum = 0
	for KEY in missDic.keys():
		missSum += missDic[KEY]
		print "miss" + "\t" +str(KEY)  + "\t" + str(missDic[KEY] )


##比例
	synlie = {}
	for KEY in sysDic.keys():
		synlie[KEY] = round(float(sysDic[KEY])/synSum,2)

	misslie = {}
	for KEY in missDic.keys():
		misslie[KEY] = round(float(missDic[KEY])/missSum,2)



	return [synlie,misslie]
def outHead():	
	headLst = []
	headLst.append("bin(<=)")
	headLst.append("xzhou")
	headLst.append("eff")
	headLst.append("propotion")
	head = "\t".join(headLst)
	return head




def outTuLine(synlie,misslie):
	OUTLST = []
	for x in range(10,100,5):		
		#syn
		outLst = []
		outLst.append(str(x)) 
		outLst.append("syn")
		outLst.append(str(synlie[x]))
		outline = "\t".join(outLst)
		OUTLST.append(outline)

		#miss
		outLst = []
		outLst.append(str(x)) 
		outLst.append("miss")
		outLst.append(str(misslie[x]))
		outline = "\t".join(outLst)
		OUTLST.append(outline)

	out = "\n".join(OUTLST)
	return out






