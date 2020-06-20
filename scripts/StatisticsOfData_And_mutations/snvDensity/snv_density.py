#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import os,sys,time
import numpy as np
import snv_density_module
from optparse import OptionParser

# The role of this script is to count snv density from the iSNVFile stat file "samples.statistics.txt"
################################################################################################################################
#
# For help as a standalone program type: python  snv_density.py   -h
#
# Examples:
#	1、python snv_density.py -i  ../../../data/iSNV_with_SNP.all.txt   -t  ../../../data/sampCtValuefile.txt  -o  ../../../data/SFig6a-snvDensity.Stat.txt
#   2、python snv_density.py  -i ../../../data/iSNV_with_SNP.all.txt   -t  ../../../data/sampCtValuefile.txt  -o  ../../../data/SFig6b-Ct29_snvDensity.Stat.txt  -C 29
#
#
#parameters:
#		-i  Path to input statistic file 'iSNVFile_with_SNP.all.txt'  generated from 'iSNVFile_calling.sh' .
#		-c  Input file of samples' Ct values
#		-o  Output file 
#		-m Minimum frequency  of the iSNV to be counted,Freq >= minFreq.
#		-M Maximum frequency  of the iSNV to be counted, Freq <maxFreq.
#		-c Minimum Ct Value  for the samples to be counted,Ct value >= minCtValue.
#		-C Maximum Ct Value  for the samples to be counted,Ct value < maxCtValue.
#	
################################################################################################################################


def main():
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputfile",
	                  dest = "inputfile",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to input statistic file 'iSNVFile_with_SNP.all.txt'  generated from 'iSNVFile_calling.sh'   [required]")
	parser.add_option("-t","--sampCtfile",
	                  dest = "sampCtfile",
	                  default = "",
	                  metavar = "file",
	                  help = "Input file of samples' Ct values  [required]")
	parser.add_option("-o","--outputfile",
	                  dest = "outputfile",
	                  default = "",
	                  metavar = "file",
	                  help = "File to save output statistic of SNP density [required]")
	parser.add_option("-m","--minFreq",
	                  dest = "minFreq",
	                  default = "0.05",
	                  metavar = "float",
	                  help = "Minimum frequency  of the iSNV to be counted,Freq >= minFreq.[optional]")
	parser.add_option("-M","--maxFreq",
	                  dest = "maxFreq",
	                  default = "0.95",
	                  metavar = "float",
	                  help = "Maximum frequency  of the iSNV to be counted, Freq <maxFreq.[optional]")
	parser.add_option("-c","--minCtValue",
	                  dest = "minCtValue",
	                  default = "15",
	                  metavar = "",
	                  help = "Minimum Ct Value  for the samples to be counted,Ct value >= minCtValue.[optional]")

	parser.add_option("-C","--maxCtValue",
	                  dest = "maxCtValue",
	                  default = "40",
	                  metavar = "",
	                  help = "Maximum Ct Value  for the samples to be counted,Ct value < maxCtValue.[optional]")
	
	
	(options,args) = parser.parse_args()

	iSNVFile     = os.path.abspath(options.inputfile)
	sampsCtF   = os.path.abspath(options.sampCtfile)	
	outF           = os.path.abspath(options.outputfile)
	minFreq     = options.minFreq
	maxFreq     = options.maxFreq
	minCtValue     = options.minCtValue
	maxCtValue     = options.maxCtValue


	startTime = time.time()
	ColsDic =  snv_density_module.iSNVFileColsDic(iSNVFile)
	sampsLst =  snv_density_module.readsampsLst(iSNVFile)
	sampleCtdic =  snv_density_module.readCt(sampsCtF)
	if (os.path.exists(outF)):
		os.remove(outF)
	print 
	Head =  snv_density_module.outHead()
	out = snv_density_module.iSNV_density(iSNVFile,ColsDic,sampsLst,sampleCtdic,minFreq,maxFreq,minCtValue,maxCtValue)

	defFO = open(outF,'a')
	defFO.write(Head + "\n")
	defFO.write(out + "\n")
	defFO.close()

	endTime = time.time()
	sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")


if __name__ == '__main__':
	main()

