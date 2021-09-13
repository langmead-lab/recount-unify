import re, sys
import glob, os
import random
import numpy as np

# Script is currently hard coded for three TSCV lists and all user parameters

# Discussed Pipeline:
# - Use suffix not prefix (both implemented)
# - Three Part Encoding
# DATA SOURCE | HASH | COUNT
# DATA SOURCE - Stores the repo location or collection location. Only relevant for distributed collections (NOT IMPLEMENTED)
# HASH - Prefix or Suffix of two alpha-numerics or three numerics (6 bits per alpha numeric, 4 bits per numeric)
# COUNT - Unique number, starting with 0. As new number added to PREFIX, this count increases. (NOT IMPLEMENTED)
# OUTPUT: A BINARY NUMBER ENCODED: "<PREFIX><COUNT>"
# PYTHON SHOULD BE ABLE TO CONVERT THAT TO LITERAL NUMBER BY PREPENDING "0b<PREFIX><COUNT>"

# Requires: module load python/2.7-anaconda 

TCGA="./tcga_samples.tsv"
GTEX="./gtex_samples.tsv"
SRA="./sra_samples.tsv"
sampleSize=5

#hashSize = Number of bits for prefix
hashSize=12 

if (hashSize % 6 != 0) or (hashSize % 4 != 0):
	print "ERROR! Hash size cannot encode both numeric and alphanumeric inputs!"

tcga_list=[]
gtex_list=[]
sra_list=[]


def strip_suffix(s, suf):
	if s.endswith(suf):
		return s[:len(s)-len(suf)]
	return s

def strip_prefix(s, pre):
	if s.startswith(pre):
		return s[len(pre):]
	return s

def strip_ends(s, b):
	s = strip_prefix(s,b)
	s = strip_suffix(s,b)
	return s

# NAME INDEX USED:
#1 GTEX RUN
#25 TCGA UUID
def readNames(fName, name_index):
	count = 0
	name_list=[]
	with open(fName) as myFile:
		for line in myFile:
			if line:
				spline = line.split("\t")
				name = spline[name_index]
				# FIRST LINE OF TSV FILE IS HEADER SO WE SKIP
				if count == 0:
					count+=1
				else:
					name_list.append(name)
	return name_list

def saveNames(outList, outName):
	with open(outName, 'w') as myFile:
		for line in outList:
			myFile.write(line+"\n")

#Converts numerics into binary of fixed length 4
#Converts alphanumeric to ascii representation (default ORD for Python 2.X), fixed length 6
def convertStringToBin(string):
	outNum=""
	#Numeric
	if len(string)==hashSize/4:
		for l in string:
			bi=bin(int(l))[2:].zfill(4)
			outNum+=str(bi)
	#AlphaNumeric
	else:
		for l in string:
			bi=format(ord(l), 'b')[1:].zfill(6)
			outNum+=str(bi)
	return outNum
	
def alphaParsePrefix(string):
	#Remove suffixes (second removal comes from the TCGA files that dont have standard UUIDS)
	noEnds = strip_ends(string,".tar.gz")
	noEnds = strip_ends(noEnds,"_rnaseq_fastq.tar")	
	noEnds = strip_ends(noEnds,"TCGA")
	noEnds = strip_ends(noEnds,"UNCID")
	
	#Split string on [^a-zA-Z0-9] ('-', '_', and '.' are all necessary) 
	spline = re.split('\W+|_', noEnds)
	
	#if SRRID then spline will have only one group and be numeric
	if len(spline)==1 and (spline[0].startswith('SRR') or spline[0].startswith('ERR') or spline[0].startswith('DRR')):
		parse=spline[0]
		parse=parse[3:]

		#numeric requires only 4 bits to encode.
		numChars = hashSize/4
		if len(parse)>=numChars:
			parse=parse[:numChars]
		else: 
			print "ERROR - hashSize ({}) [{}] larger than string ({}) for numerics.".format(hashSize,numChars, len(parse))
	else:
		#alphanumeric requires 6 bits to encode
		numChars = hashSize/6
		reString=""
		for subs in spline:
			reString+=subs

		if len(reString)>=numChars:
			parse=reString[:numChars]
		else:
			print "ERROR - hashSize ({}) [{}] larger than string ({}) for alphanumerics".format(hashSize,numChars, len(reString))
	
	#print string
	#print parse
	return parse
	
def alphaParseSuffix(string):
	#Remove suffixes (second removal comes from the TCGA files that dont have standard UUIDS)
	noEnds = strip_ends(string,".tar.gz")
	noEnds = strip_ends(noEnds,"_rnaseq_fastq.tar")	
	noEnds = strip_ends(noEnds,"TCGA")
	noEnds = strip_ends(noEnds,"UNCID")

	#Split string on '-' and '_'
	spline = re.split('\W+|_', noEnds)
	
	#if SRRID then spline will have only one group and be numeric
	if len(spline)==1 and (spline[0].startswith('SRR') or spline[0].startswith('ERR') or spline[0].startswith('DRR')):
		parse=spline[-1]
		parse=parse[3:]
		
		#numeric requires only 4 bits to encode.
		numChars = hashSize/4
		if len(parse)>=numChars:
			parse=parse[-numChars:]
		else: 
			print "ERROR - hashSize ({}) [{}] larger than string ({}) for numerics.".format(hashSize,numChars, len(parse))
	else:
		#alphanumeric requires 6 bits to encode
		numChars = hashSize/6
		reString=""
		for subs in spline:
			reString+=subs

		if len(reString)>=numChars:
			parse=reString[-numChars:]	
		else:
			print "ERROR - hashSize ({}) [{}] larger than string ({}) for alphanumerics".format(hashSize,numChars, len(reString))
	
	return parse	

#Format Specific alphaParse (Test method of splitting TCGA from SRA labeling) [Data source]
def addDataSource(string, function):
	spline = re.split('\W+|_', string)
	if len(spline)==1:
		parse = "0"+function(string)
	else:
		parse = "1"+function(string)
	return parse

def binResults(function):
	outBin={}
	for name in tcga_rand:
		parse=function(name)
		if parse in outBin:
			outBin[parse]+=1
		else:
			outBin[parse]=1
	for name in gtex_rand:
		parse=function(name)
		if parse in outBin:
			outBin[parse]+=1
		else:
			outBin[parse]=1
	for name in sra_rand:
		parse=function(name)
		if parse in outBin:
			outBin[parse]+=1
		else:
			outBin[parse]=1
	return outBin

def encodeResults(function):
	outBin={}
	for name in tcga_rand:
		parse=function(name)
		print "{}: {:12s}".format(name,convertStringToBin(parse))
	for name in gtex_rand:
		parse=function(name)
		print "{}: {:12s}".format(name,convertStringToBin(parse))
	for name in sra_rand:
		parse=function(name)
		print "{}: {:12s}".format(name,convertStringToBin(parse))

def analyzeResults(outBin):
	dist=[]
	for k in outBin.keys():
		dist.append(outBin[k])

	print len(outBin.keys())
	print np.mean(dist)
	print np.std(dist)

#s is user string
#function is either alphaParsePrefix or alphaParseSuffix
def encodeString(s, function):
	return convertStringToBin(function(s))

if __name__ == '__main__':
    # MAIN
    #Parse names from later tsv files
    tcga_list=readNames(TCGA,25)
    gtex_list=readNames(GTEX,1)
    sra_list=readNames(SRA,1)

    #Select random subset (sampleSize)
    tcga_rand=random.sample(tcga_list,sampleSize)
    gtex_rand=random.sample(gtex_list,sampleSize)
    sra_rand=random.sample(sra_list,sampleSize)

    # Print num in subset versus total number in set
    print "{}, {}".format(len(tcga_list),len(tcga_rand))
    print "{}, {}".format(len(gtex_list),len(gtex_rand))
    print "{}, {}".format(len(sra_list),len(sra_rand))

    #save samples (each run randomly reselects so this gives a 'backup' of the most recent run and overwrites each time)
    saveNames(tcga_rand, "tcga_sample.txt")
    saveNames(gtex_rand, "gtex_sample.txt")
    saveNames(sra_rand, "sra_sample.txt")

    #encoding strategy 1: Unique Header Parse (2 alphanumeric)
    print "PREFIX"
    prefix = encodeResults(alphaParsePrefix)
    print "SUFFIX"
    suffix = encodeResults(alphaParseSuffix)

    print int(encodeString(tcga_list[0], alphaParseSuffix), 2)
    print int(encodeString(gtex_list[0], alphaParseSuffix), 2)
    print int(encodeString(sra_list[0], alphaParseSuffix), 2)

    #analyzeResults(outBin)
    #analyzeResults(fsap_outBin)
