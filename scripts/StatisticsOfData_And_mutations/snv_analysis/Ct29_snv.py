#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys,os
from optparse import OptionParser
import Ct29_snv_module  as geneM
#
# This script is used to count the base substitution types for singleton and shared iSNV sites and to generate files showing their frequency distribution 
################################################################################################################################
#
# For help as a standalone program type: python Ct29_snv.py  -h
#
# Examples:
#	   1、python Ct29_snv.py   -i  ../../../data/iSNV_with_SNP.all.txt  -t  ../../../data/sampCtValuefile.txt  -e  ../../../data/allsampMutatAnno.txt  -g share -o ../../../result/
#      2、python Ct29_snv.py   -i  ../../../data/iSNV_with_SNP.all.txt  -t  ../../../data/sampCtValuefile.txt  -e  ../../../data/allsampMutatAnno.txt  -g single -o ../../../result/
#		
#	
################################################################################################################################

def main():
	SAMPIDs = geneM.sample_ID(iSNVTable)
	sampleCtdic = geneM.sampCt(sampsCtF)
	print "The Number of samples that met the requirements :"
	print len(sampleCtdic)


	singleXins = geneM.outAllSampFStat(outAllSampF,sampleCtdic,single_shareFlag)
	singlePosL = singleXins[0]
	singlePosD = singleXins[1]
	geneLenD = geneM.RegionLen()
	Statouts = geneM.AllSampFStat(outAllSampF,singlePosL,singlePosD)
	Htu_snvDic = Statouts
	
	H_snv_outF = outputpath + "/Fig4fg-" + single_shareFlag + "--snv-effStat.txt" 
	print H_snv_outF
	Fexit = os.path.exists(H_snv_outF)
	if Fexit == True:
		os.remove(H_snv_outF)
	geneM.BaseChgnOut(Htu_snvDic,geneLenD,H_snv_outF)


	RplotsingleDataF = outputpath  + "/Fig4hi-" + single_shareFlag + "-RplotFreqboxDataF.txt" 
	Fexit = os.path.exists(RplotsingleDataF)
	print  RplotsingleDataF
	if Fexit == True:
		os.remove(RplotsingleDataF)
	
	geneM.RplotsingleData(outAllSampF,singlePosL,singlePosD,RplotsingleDataF,single_shareFlag)



if __name__ == "__main__":
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
	                  help = "Path to Ct value table [required]")
	parser.add_option("-e","--effsnp",
	                  dest = "effsnp",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to All Samp snpeff stat file [required]")

	parser.add_option("-g","--single_shareFlag",
	                  dest = "single_shareFlag",
	                  default = "",
	                  action = "store",
	                  type = "string",
	                  help = "Flag of snv-sample type :[single | share]  [required]")

	parser.add_option("-o","--outputpath",
	                  dest = "outputpath",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to save output files in [required]")

	(options,args) = parser.parse_args()

	iSNVTable       = os.path.abspath(options.inputfile)
	sampsCtF       = os.path.abspath(options.sampCtfile)	
	outAllSampF      = os.path.abspath(options.effsnp)
	single_shareFlag = options.single_shareFlag
	outputpath       = os.path.abspath(options.outputpath)
	if (not os.path.exists(outputpath)):
		os.mkdir(outputpath)
	main()












