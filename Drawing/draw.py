import argparse
import subprocess
import os
import collections
import random
import operator
import copy



## Sort the tree by abundance on class level to match the colors
def sortByAbundance(tree, level, currentLevel):
	if currentLevel<level:
		return tree;
	treeDict=dict()
	for clade in tree:
		abundance=tree[clade]
		sepIdx= clade.rfind(".")
		newClade=clade[0:sepIdx]
		#print newClade+ "\t"+str(abundance)
		if newClade in treeDict:
			treeDict[newClade]=treeDict[newClade]+abundance
		else:
			treeDict[newClade]=abundance
	currentLevel=currentLevel-1
	return sortByAbundance (treeDict, level, currentLevel)

### Generates annot1 file
def writeAnnot1(tree, file, level):
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
		file.write(clade[(sepIdx+1):len(clade)]+"\tclade_marker_size\t"+str((abundance*100)+5)+"\n")	
	level=level-1	
	writeAnnot1(newTree,file, level)


## Generates annot2 file
def writeAnnot2(tree, file):
	colors=["#FF0A00", "#FF5F00","#FF8700","#FF9B00","#FFB400","#FFC800","#FFFF00"] ## Colors from red to yellow. Red is the most abundant and yellow is the least abundant. 7 colors because we have only 7 classes in our sample.
	levels={"Species":4, "Genera":3, "Family":2, "Order":1, "Class":0}
	colorIdx=0
	cladeColorDic=dict()
	writtenClades=set()	
	speciesToWrite={'Bacteroidia'}
	cladesToCut={'Faecalibacterium':' G1:Faecalibacterium', 'Parabacteroides':' G2:Parabacteroides','unclassified Bifidobacteriales (miscellaneous)':' G3:un.Bifidobacteriales', 'unclassified Bifidobacteriales':' F1:unclassified Bifidobacteriales'}
	
	start=2		# From class level to species level
	end=7
	letters=map(chr, range(65, 100)) ## Letters list to replace the species names
	letterIdx=0	## Idx over the letters list
	for branch in tree:
		color=None
		includeSpecies=0
		clades=branch[0].split(".")
		abundance=branch[1]
		font=8
		## Check if to write the species in Bacteroidia. If not, it won't be writtin
		for clade in speciesToWrite:
			if clade in branch[0]:
				includeSpecies=1
		
		## Setting the color
		if clades[start] in cladeColorDic:
			color=cladeColorDic[clades[start]]
		else:
			color=colors[colorIdx]
			colorIdx=colorIdx+1
			cladeColorDic[clades[start]]=color

		## looping over the clades in each branch	
		for idx, clade in enumerate(clades[start:end]):
			if font>4:
				font=font-1
			if clade not in writtenClades:
				cladeAnnot="_"
				if idx!=levels["Species"]: ## Level higher than species
					if clade in cladesToCut:
						cladeAnnot=cladesToCut[clade]
					else:
						cladeAnnot=clade						
				else:		## Species Level
					if includeSpecies==1:
						cladeAnnot=letters[letterIdx]+":"+clade[0:1]+"."+clade[clade.find('_')+1:len(clade)]
						letterIdx=letterIdx+1					
				# print str(idx)+"\t"+cladeAnnot
				file.write(clade+"\tannotation\t"+cladeAnnot+"\n")
				file.write(clade+"\tannotation_font_size\t"+str(font)+"\n")
				file.write(clade+"\tannotation_background_color\t"+color+"\n")		
				writtenClades.add(clade)


# set input parameters
parser = argparse.ArgumentParser(description="This script takes as input MetaFlow's output csv file and runs graphlan to produce a tree-of-life represenation of the species. \n Ver. 1.0")

group = parser.add_argument_group("Required arguments")
group.add_argument("-i", dest="inputFileName", help="the abundance.CSV output of MetaFlow", required=True)
group.add_argument("-g", dest="pathToGraPhlAn", help="the path to where GraPhlAn is located ", required=True)

args = parser.parse_args()

inputFilename = args.inputFileName;
graphlanPath = args.pathToGraPhlAn;
outputFileName = os.path.splitext(inputFilename)[0] + ".txt";

dir_to_stript = os.path.dirname(os.path.abspath(__file__))
with open(dir_to_stript + '/FlatTaxonomy.txt') as f1:
    speciesFull = f1.read().splitlines()

with open(inputFilename) as f2:
    results = f2.read().splitlines()

outputFile = open(outputFileName,'w')


treeDict=dict()
for species in results[1:]:
	speciesName = species.split('\t')[1].replace("."," ")
	speciesAbundance = float(species.split('\t')[6].strip())
	speciesPath = ""
	for sF in speciesFull:
		tempName = speciesName.split('_')[0]
		tempSF = sF.replace('_',' ')
		if tempSF.startswith(tempName):
			speciesPathList = sF.split('\t')[1].split(';')
			# removing trailing ;
			speciesPathList = speciesPathList[:-1] 
			for clade in reversed(speciesPathList):
				if "group" in clade:
					continue
				speciesPath += clade + '.'
			# removing trailing .
			speciesPath = speciesPath[:-1]+"."+speciesName
				
			break;
	if speciesPath != "":
		treeDict[speciesPath]=speciesAbundance		
	else:
		print 'NOT FOUND: '+speciesName+'\n'
		outputFile.write('NOT FOUND:'+speciesName+'\n')




sortedByKeyTree=collections.OrderedDict(sorted(treeDict.items()))
#sortedTree=sorted(treeDict.keys())
for branch in sortedByKeyTree.keys():
	outputFile.write(branch + '\n')
outputFile.close()

stage1File=open("annot_1.txt", "w")
writeAnnot1(sortedByKeyTree, stage1File, 7)
stage1File.close()	

		


sortedByAbund=sortByAbundance(treeDict, 4, 7)
sortedByValTuple= sorted(treeDict.items(), key=operator.itemgetter(1), reverse=True)
sortedTreeByAbund=dict()
for atuple in sortedByValTuple:
	branch=atuple[0]
	for clade in sortedByAbund:
		if clade in branch:
			sortedTreeByAbund[branch]=sortedByAbund[clade]
			break
sortedByValTuple= sorted(sortedTreeByAbund.items(), key=operator.itemgetter(1), reverse=True)
stage2File=open("annot_2.txt", "w")
writeAnnot2(sortedByValTuple, stage2File)
stage2File.close()	


######################### LAUNCHING GRAPHLAN ######################
subprocess.call("python " + graphlanPath + "/graphlan_annotate.py --annot annot_1.txt " + outputFileName + " " + os.path.splitext(inputFilename)[0] + "_1.xml", shell = True)
subprocess.call("python " + graphlanPath + "/graphlan_annotate.py --annot annot_2.txt " + os.path.splitext(inputFilename)[0] + "_1.xml " + os.path.splitext(inputFilename)[0] + "_2.xml", shell = True)

subprocess.call("python " + graphlanPath + "/graphlan.py " + os.path.splitext(inputFilename)[0] + "_2.xml " + os.path.splitext(inputFilename)[0] + ".svg --size 3.5", shell = True)

os.remove("annot_1.txt");
os.remove("annot_2.txt");
os.remove(outputFileName);
os.remove(os.path.splitext(inputFilename)[0] + "_1.xml");
os.remove(os.path.splitext(inputFilename)[0] + "_2.xml");


