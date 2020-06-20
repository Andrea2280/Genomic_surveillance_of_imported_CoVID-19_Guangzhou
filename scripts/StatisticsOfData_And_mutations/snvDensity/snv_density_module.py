#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import os,sys,linecache
import numpy as np

def iSNVFileColsDic(iSNV):
	Head = linecache.getline(iSNV,3)
	HeadLst = Head.split("\t")
	ColsDic = {}
	Count = 0
	for Headsamp in HeadLst:
		ColsDic[Headsamp] = Count
		Count += 1
	return ColsDic


def readsampsLst(iSNV):
	samps = linecache.getline(iSNV,3)
	sampsLst = samps.split("\n")[0].split("\t")[2:-1]
	return sampsLst



def readCt(Person_SIDF):
	sampleCtdic = {}
	Person_SIDFls = open(Person_SIDF).readlines()
	for Person_SIDFl in Person_SIDFls:
		if "Patient" not in Person_SIDFl and Person_SIDFl != "\n":
			sample_SID = Person_SIDFl.split("\t")[0]
			Ct = Person_SIDFl.split("\t")[4].split("\n")[0]
			sampleCtdic[sample_SID] = Ct		
	return sampleCtdic

def iSNV_density(iSNV,ColsDic,sampsLst,sampleCtdic,minFreq,maxFreq,minCtValue,maxCtValue):
	sampiSNVdic = {}
	for sample in sampsLst:
		sampiSNVdic[sample] = 0

	iSNVls = open(iSNV).readlines()
	for iSNVl in iSNVls:
		if "#" not in iSNVl  and iSNVl != "\n":
			lineLst = iSNVl.split("\t")
			for sample in sampsLst:
				lieIndex = ColsDic[sample]
				sampFreq = lineLst[lieIndex]
				if sampFreq != "NO" and sampFreq != "NA" and float(sampFreq) >= float(minFreq) and  float(sampFreq) < float(maxFreq) :
					sampiSNVdic[sample] += 1 

	outlist = []
	sampCount = 0
	for sample in sampsLst:			
		Ct = sampleCtdic[sample]
		if Ct != "NA"  and  float(Ct)  >= float(minCtValue)  and  float(Ct)  < float(maxCtValue)   :	
			sampCount += 1	
			iSNVCount = sampiSNVdic[sample]
			iSNVdensity = round(float(iSNVCount)/29,3)
			OUTLst = []
			OUTLst.append(sample)
			OUTLst.append(Ct)
			OUTLst.append(str(iSNVCount))
			OUTLst.append(str(iSNVdensity))
			OUT = "\t".join(OUTLst)
			outlist.append(OUT)
	out = "\n".join(outlist)
	print  out
	print 
	print "sample Num(Ct:" + str(minCtValue)  + "-" + str(maxCtValue) + ") :  "  +  str(sampCount)
	return out




def outHead():	
	headLst = []
	headLst.append("samp")
	headLst.append("Ct")
	headLst.append("iSNVCount")
	headLst.append("iSNVdensity")
	head = "\t".join(headLst)
	print head
	return head




	


