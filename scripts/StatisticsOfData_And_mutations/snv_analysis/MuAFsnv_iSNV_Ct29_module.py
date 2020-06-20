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
	for Ct25Posi in CtPosiSampDic.keys():
		shareNumDic[Ct25Posi] = len(CtPosiSampDic[Ct25Posi])
		if flag == "share-":		
			if len(CtPosiSampDic[Ct25Posi]) > 1:
				Ct25singlePosL.append(int(Ct25Posi))	
			else:
				CtsinglePosD.pop(Ct25Posi)

		if flag == "single-":		
			if len(CtPosiSampDic[Ct25Posi]) == 1:
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


def iSNV_shareSnv(iSNV,lieDic,sampsLst,lowFreq,upFreq):
	sampSNPdic = {}
	shareNumDic = {}
	iSNVls = open(iSNV).readlines()
	for iSNVl in iSNVls:
		if "#" not in iSNVl  and iSNVl != "\n":
			lineLst = iSNVl.split("\t")
			posi = lineLst[0]			
			if int(posi) > 86 and int(posi) < 29782 and int(posi) != 12436:
				count = 0
				sampSNPdic[posi] = []
				for sample in sampsLst:
					lieIndex = lieDic[sample]
					Freq = lineLst[lieIndex]					
					if Freq != "NO" and Freq != "NA"  and float(Freq) >= lowFreq  and float(Freq) < upFreq:
						count += 1
						sampSNPdic[posi].append(float(Freq))
				shareNumDic[posi] = count

	outlist = []
	for position in sampSNPdic.keys():	
		count = len(sampSNPdic[position])
		flagShare = ''			
		if count == 1:
			flagShare = "single"
		if count > 1 and count <= 10:
			flagShare = "share"			
		if flagShare != '' :		
			for Freq in sampSNPdic[position]:
				OUTLst = []
			
				#OUTLst.append(flag)
				OUTLst.append(position)
				OUTLst.append(str(count))
				OUTLst.append(str(Freq))	
				OUTLst.append(flagShare)			
				OUT = "\t".join(OUTLst)
				#if flagShare == "single":
				
				outlist.append(OUT)
	
	out = "\n".join(outlist)
	return out



def boxplot_outHead():	
	headLst = []
	headLst.append("snpeffFname")
	headLst.append("pos")
	headLst.append("eff")
	headLst.append("Freq")
	headLst.append("flag")
	boxplot_head = "\t".join(headLst)
	return boxplot_head

def boxplot_syn_miss(sampsLst,lieDic,iSNV,effDic,shareNumDic,lowFreq,upFreq,singlePos,flag):
	FreqLst = []
	iSNVls = open(iSNV).readlines()
	for iSNVl in iSNVls:
		if "#" not in iSNVl  and iSNVl != "\n":
			lineLst = iSNVl.split("\t")
			posi = lineLst[0]
			if int(posi) in  singlePos: #				
				for sample in sampsLst:
					if sample in  lieDic:
						lieIndex = lieDic[sample]
						Freq = lineLst[lieIndex]						
						if Freq != "NO" and Freq != "NA"  and float(Freq) > float(lowFreq)  and float(Freq) <= float(upFreq) :
							if sample + "_" + posi in effDic:
								eff = effDic[sample + "_" + posi]	
							if eff in ["synonymous_variant","missense_variant","stop_gained"]:
								if eff == "synonymous_variant":
									ann = flag + "Syn"
									clor = "Syn"
								if eff == "missense_variant"  or eff == "stop_gained":
									ann = flag + "Nonsyn"
									clor = "Nonsyn"

								hanglst = []								
								hanglst.append(sample)
								hanglst.append(posi)
								hanglst.append(ann)
								hanglst.append(Freq)
								hanglst.append(clor)
								hang = "\t".join(hanglst)
								FreqLst.append(hang)
	boxout = "\n".join(FreqLst)
	return boxout




