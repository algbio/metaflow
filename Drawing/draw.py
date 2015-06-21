import argparse
import subprocess
import os

# set input parameters
parser = argparse.ArgumentParser(description="This script takes as input MetaFlow's output csv file and runs graphlan to produce a tree-of-life represenation of the species. \n Ver. 1.0")

group = parser.add_argument_group("Required arguments")
group.add_argument("-i", dest="inputFileName", help="the abundance.CSV output of MetaFlow", required=True)
group.add_argument("-o", dest="outputFileName", help="output fileName", required=True)

args = parser.parse_args()

inputFilename = args.inputFileName;
outputFileName = args.outputFileName;

with open('FlatTaxonomy.txt') as f1:
    speciesFull = f1.read().splitlines()

with open(inputFilename) as f2:
    results = f2.read().splitlines()

outputFile = open(outputFileName,'w')


for species in results[1:]:
	speciesName = species.split('\t')[1]
	speciesAbundance = species.split('\t')[6]
	speciesPath = ""
	for sF in speciesFull:
		tempName = speciesName.replace('_',' ')
		tempSF = sF.replace('_',' ')
		if tempSF.startswith(tempName):
			speciesPathList = sF.split(';')
			# removing trailing ;
			speciesPathList = speciesPathList[:-1] 
			for clade in reversed(speciesPathList):
				speciesPath += clade + '.'
			# removing trailing .
			speciesPath = speciesPath[:-1]
			break;
	if speciesPath != "":
		outputFile.write(speciesPath + '\n')
	else:
		outputFile.write('NOT FOUND\n')

outputFile.close()
