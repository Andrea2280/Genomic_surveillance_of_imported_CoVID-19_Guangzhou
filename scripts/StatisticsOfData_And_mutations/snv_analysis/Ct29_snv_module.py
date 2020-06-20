#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys,os
import numpy as np


def sample_ID(iSNVTable):
	iSNVTablels = open(iSNVTable,'r').readlines()
	for iSNVTablel in iSNVTablels:
		if "pos	" in iSNVTablel:
			sampIDs = iSNVTablel.split("\n")[0].split("\t")[2:-1]			
	return sampIDs 


def charRank_geneRegion(gene,posi):
	refGeneDic = {"orf1a":"266..13468","orf1b":"13468..21555","S":"21563..25384"\
	,"ORF3a":"25393..26220","E":"26245..26472","M":"26523..27191",\
	"ORF6":"27202..27387","ORF7a":"27394..27759","ORF8":"27894..28259",\
	"N":"28274..29533","ORF10":"29558..29674",}
	for refGene in refGeneDic:
		geneIDStrt = int(refGeneDic[refGene].split("..")[0])
		geneIDEnd = int(refGeneDic[refGene].split("..")[1])
		if int(posi)  >= geneIDStrt  and int(posi)  <= geneIDEnd:
			charRank = (int(posi)-geneIDStrt )%3 + 1
			break
	return charRank

def RegionLen():
	refGeneDic = {"orf1a":"266..13468","orf1b":"13468..21555","S":"21563..25384"\
	,"ORF3a":"25393..26220","E":"26245..26472","M":"26523..27191",\
	"ORF6":"27202..27387","ORF7a":"27394..27759","ORF8":"27894..28259",\
	"N":"28274..29533","ORF10":"29558..29674",}
	geneLenD = {}
	codingLen = 0
	for refGene in refGeneDic:
		geneIDStrt = int(refGeneDic[refGene].split("..")[0])
		geneIDEnd = int(refGeneDic[refGene].split("..")[1])
		genelen = geneIDEnd - geneIDStrt + 1
		geneLenD[refGene] = genelen
		codingLen += genelen

	NoncodingLen = 29903-codingLen
	geneLenD['Coding'] = codingLen
	geneLenD['Non-Coding'] = NoncodingLen
	return geneLenD



def sampCt(sampsCtF):
	sampleCt29dic = {}
	sampsCtFls = open(sampsCtF).readlines()
	for Person_SIDFl in sampsCtFls:
		if "Patient" not in Person_SIDFl and Person_SIDFl != "\n":
			sample_SID = Person_SIDFl.split("\t")[0]
			Ct = Person_SIDFl.split("\t")[4].split("\n")[0]
			if Ct != 'NA' and float(Ct) < 29 :
				sampleCt29dic[sample_SID] = Ct
	return sampleCt29dic

	


def outAllSampFStat(outAllSampF,sampleCt29dic,single_shareFlag):	
	print sampleCt29dic
	Ct29PosiSampDic = {}
	outAllSampFls = open(outAllSampF).readlines()
	for outAllSampFl in outAllSampFls:
		if "sample" not in outAllSampFl and outAllSampFl != "\n":
			linekey = outAllSampFl.split("\n")[0].split("\t")
			sample = linekey[0].split("/")[-1]
			
			if sample in sampleCt29dic.keys():
				posi = linekey[1]
				Freq = linekey[17]
				if int(posi) < 29782 and int(posi) > 86 and int(posi) != 12436:
					if Freq != "NO" and Freq != "NA" and float(Freq) < 95:						
						if posi not in Ct29PosiSampDic:				
							Ct29PosiSampDic[posi] = [sample]
						else :					
							Ct29PosiSampDic[posi].append(sample)

	Ct29singlePosL = []
	Ct29singlePosD = Ct29PosiSampDic
	for Ct29Posi in Ct29PosiSampDic.keys():
		if single_shareFlag == "single":			
			if len(Ct29PosiSampDic[Ct29Posi]) == 1:
				Ct29singlePosL.append(int(Ct29Posi))
			else:
				Ct29singlePosD.pop(Ct29Posi)
		if single_shareFlag == "share":
			if len(Ct29PosiSampDic[Ct29Posi]) > 1:
				Ct29singlePosL.append(int(Ct29Posi))
			else:
				Ct29singlePosD.pop(Ct29Posi)

	Ct29singlePosL.sort()
	return Ct29singlePosL,Ct29singlePosD



def Htu(Htu_snvDic,Annotation,HGVS_c):
	if Annotation in ['missense_variant','stop_gained']:
		AnnFlag = "Nonsyn"
	elif Annotation in ['synonymous_variant']:
		AnnFlag = "Syn"
	else:
		AnnFlag = "Non-Coding"

	if HGVS_c not in Htu_snvDic:
		Htu_snvDic[HGVS_c] = {}
		if AnnFlag not in Htu_snvDic[HGVS_c]:
			Htu_snvDic[HGVS_c][AnnFlag] = 1
	else:
		if AnnFlag not in Htu_snvDic[HGVS_c]:
			Htu_snvDic[HGVS_c][AnnFlag] = 1
		else:
			Htu_snvDic[HGVS_c][AnnFlag] += 1
	return Htu_snvDic						


		
def AllSampFStat(outAllSampF,Ct29singlePosL,Ct29singlePosD):
	Htu_snvDic = {}
	outAllSampFls = open(outAllSampF).readlines()
	for outAllSampFl in outAllSampFls:
		if "sample" not in outAllSampFl and outAllSampFl != "\n":			
			linekey = outAllSampFl.split("\n")[0].split("\t")
			sample = linekey[0].split("/")[-1]
			posi = linekey[1]
			if int(posi) in Ct29singlePosL and sample in Ct29singlePosD[posi] :	##single 			
				ref = linekey[2]
				alle = linekey[3]
				if ref == 'T':
					ref = "U"
				if alle == 'T':
					alle = "U"

				Annotation = linekey[4]
				Gene_Name = linekey[5]
				Feature_Type = linekey[6]
				Transcript_BioType = linekey[7]
				rank = linekey[8]
				rankall = linekey[9]
				HGVS_c_1 = linekey[10]
				HGVS_c = ref + ">" + alle
				HGVS_p = linekey[11]
				cDNA_pos = linekey[12]
				cDNA_posall = linekey[13]
				AA_pos = linekey[14]
				AA_posall = linekey[15]
				charRank = linekey[16]
				Freq = linekey[17]				
				if Gene_Name == "orf1ab":
					if int(posi) >= 266 and int(posi) <= 13468:
						Gene_Name = "orf1a"
					if int(posi) >= 13468 and int(posi) <= 21555:
						Gene_Name = "orf1b"
				if float(Freq) < 95:   
					Htu_snvDic = Htu(Htu_snvDic,Annotation,HGVS_c)

	return Htu_snvDic


def BaseChgnOut(Htu_Dic,geneLenD,H_outF):  ##isnv/kb
	H_outFO = open(H_outF,'a')
	H_outFO.write("baseChange" +  "\t" + "Nonsyn_Num" +  "\t" + "Syn_Num" +  "\t" + "Non_Coding_Num"  + "\n")
	H_outFO.close()
	print Htu_Dic

	outLst = []
	for Htu_key in Htu_Dic.keys():
		lineLst = []
		lineLst.append(Htu_key)		
		for Charposi in ['Nonsyn','Syn',"Non-Coding"]:
			if Charposi in  Htu_Dic[Htu_key]:
				CharposiCout = Htu_Dic[Htu_key][Charposi]
			else:
				CharposiCout = 0
			lineLst.append(str(CharposiCout))				
			#lineLst.append(str(round(float(CharposiCout)/geneLenD[Charposi],3)))
		outline = "\t".join(lineLst)
		outLst.append(outline)
	out = "\n".join(outLst)
	print out
	H_outFO = open(H_outF,'a')
	H_outFO.write(out+ "\n")
	H_outFO.close()



def RplotsingleData(outAllSampF,Ct29singlePosL,Ct29singlePosD,RplotsingleDataF,Flag):	
	outAllSampFls = open(outAllSampF).readlines()
	outsinglesnpeffFO = open(RplotsingleDataF,'a')
	outsinglesnpeffFO.write("posi" + "\t" + "sample" + "\t" + "Flag"  + "\t" + "Freq"  +"\t" + "BaseChange" + "\n")
	outsinglesnpeffFO.close()
	for outAllSampFl in outAllSampFls:
		if "sample" not in outAllSampFl and outAllSampFl != "\n":
			linekey = outAllSampFl.split("\n")[0].split("\t")
			sample = linekey[0].split("/")[-1]
			posi = linekey[1]
			Freq = linekey[17]
			ref = linekey[2]
			alle = linekey[3]
			if ref == 'T':
				ref = "U"
			if alle == 'T':
				alle = "U"
			cChange=  ref + ">" + alle
			print cChange
			if int(posi) in Ct29singlePosL and  sample in Ct29singlePosD[posi]:								
				outsinglesnpeffFO = open(RplotsingleDataF,'a')
				outsinglesnpeffFO.write(posi + "\t" + sample + "\t" + Flag  + "\t" + Freq  +"\t" + cChange + "\n")
				outsinglesnpeffFO.close()
				









