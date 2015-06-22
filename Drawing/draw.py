import argparse
import subprocess
import os



def writeAbundance(tree, file):
	if not tree:
		return
	newTree=dict()
	sepIdx=-1	
	for clade in tree:		
		abundance=tree[clade]
		if "." in clade:
			sepIdx= clade.rfind(".")
			newClade=clade[0:sepIdx]
			if newClade in newTree:
				newTree[newClade]=newTree[newClade]+abundance
			else:
				newTree[newClade]=abundance
		else:
			sepIdx= -1					
		file.write(clade[(sepIdx+1):len(clade)]+"\t"+str(abundance)+"\n")		
	writeAbundance(newTree,file)

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


fullTree=dict()
for species in results[1:]:
	speciesName = species.split('\t')[1]
	speciesAbundance = species.split('\t')[6].strip()
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
		if speciesPath.count(".")==5:	##Include only species with 7 levels. More than 7 levels Should be fixed by hand.
			fullTree[speciesPath]=float(speciesAbundance)
	else:
		outputFile.write('NOT FOUND\n')
outputFile.close()

test=open("Test_Example","w")


fixedFullTree=dict()
for clade in fullTree:			## Put the genera in the right location in the path
	splitClade=clade.split("\t")
	path=splitClade[0]
	genera=splitClade[1]
	sepIdx= clade.rfind(".")
	fullpath=path[0:sepIdx]+"."+genera+"."+path[sepIdx+1:len(path)]
	fixedFullTree[fullpath]=fullTree[clade]
	print (fullpath)
writeAbundance(fixedFullTree,test)	## Recursive call for writing the sum on every clade
test.close()
