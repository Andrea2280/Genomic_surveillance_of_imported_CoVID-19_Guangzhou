#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import sys,os,linecache
from optparse import OptionParser
import numpy as np

import MuAFSpectrum_Ct29_module  as MuAFSpectrumMOD
#1



# The script is used to view the MuAF Spectrum of singleton and shared iSNV sites from the iSNVtable.
################################################################################################################################
#
# For help as a standalone program type: python MuAFSpectrum_Ct29.py   -h
#
# Examples:
#	   python MuAFSpectrum_Ct29.py  -i  ../../../data/iSNV_with_SNP.all.txt  -t  ../../../data/sampCtValuefile.txt  -e  ../../../data/allsampMutatAnno.txt  -o ../../../data  
#		
# From "iSNV_with_SNP.txt"	
################################################################################################################################



def main(): 
	for flag in ['single','shared']:

		sampleCtdic = MuAFSpectrumMOD.sampCt(sampsCtF,minCtValue,maxCtValue)
		singleXins = MuAFSpectrumMOD.outAllSampFStat(allSampsAnno,sampleCtdic,flag)
		sampsLst = sampleCtdic.keys()
		PosiPosi = singleXins[0]
		PosiDict = singleXins[1]
		shareNumDic = singleXins[2]
		effDic = singleXins[3]
		lieDic = MuAFSpectrumMOD.iSNVlieDic(iSNV)

		tongji = MuAFSpectrumMOD.tongjisnpeff(sampsLst,iSNV,lieDic,effDic,shareNumDic,lowFreq,upFreq,PosiPosi)
		synlie = tongji[0]
		misslie = tongji[1]
		print synlie	
		Head = MuAFSpectrumMOD.outHead()
		out = MuAFSpectrumMOD.outTuLine(synlie,misslie)
	
		lineplot_outF = lineplot_outP + "Fig4de-MuAFSpectrum_Ct29--" + flag  + ".txt"
		lineplot_outFexit = os.path.exists(lineplot_outF)
		if lineplot_outFexit == True:
			os.remove(lineplot_outF)
		lineplot_outFO = open(lineplot_outF,'a')
		lineplot_outFO.write(Head + "\n")
		lineplot_outFO.write(out + "\n")
		lineplot_outFO.close()



if __name__ == '__main__':
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputfile",
	                  dest = "inputfile",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to input iSNV Table [required]")
	parser.add_option("-t","--sampCtfile",
	                  dest = "sampCtfile",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to Ct value Table [required]")
	parser.add_option("-e","--effsnp",
	                  dest = "effsnp",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to All Samp snpeff stat file [required]")

	parser.add_option("-o","--lineplot_outP",
	                  dest = "lineplot_outP",
	                  default = "",
	                  metavar = "path",
	                  help = "path of output files  [required]")

	parser.add_option("-m","--minMuAF ",
	                  dest = "minMuAF",
	                  default = "0.05",
	                  metavar = "",
	                  help = "Minimum MuAF for the iSNV cite to be counted.[optional]")

	parser.add_option("-M","--maxMuAF",
	                  dest = "maxMuAF",
	                  default = "0.95",
	                  metavar = "",
	                  help = "Maximum MuAF for the iSNV cite to be counted.[optional]")

	parser.add_option("-c","--minCtValue",
	                  dest = "minCtValue",
	                  default = "10",
	                  metavar = "",
	                  help = "Minimum Ct Value  for the samples to be counted,Ct value >= minCtValue.[optional]")

	parser.add_option("-C","--maxCtValue",
	                  dest = "maxCtValue",
	                  default = "29",
	                  metavar = "",
	                  help = "Maximum Ct Value  for the samples to be counted,Ct value < maxCtValue.[optional]")



	(options,args) = parser.parse_args()

	iSNV           = os.path.abspath(options.inputfile)
	sampsCtF       = os.path.abspath(options.sampCtfile)	
	allSampsAnno   = os.path.abspath(options.effsnp)
	lineplot_outP  = options.lineplot_outP
	lowFreq        = options.minMuAF
	upFreq         = options.maxMuAF
	minCtValue     = options.minCtValue
	maxCtValue     = options.maxCtValue

	
	
 


	main()