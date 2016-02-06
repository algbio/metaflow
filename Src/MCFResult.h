/*
 * MCFResult.h
 *
 *  Created on: 16.1.2014
 *      Author: sobih
 */
#ifndef MCFRESULT_H_
#define MCFRESULT_H_


#include "MCFHeaders.h"
#include "MCFGenomeInfo.h"
#include "MCFLogger.h"

extern MCFLogger mcfLogger;
/*
 * A struct contains the information related to the reads
 */
struct ReadInfo{
	int numOfMappedReads; 			// Total number of reads that have Blast hit in the graph. The value is updated when reads are removed or assigned to unknown.
	int numOfSingleHitReads;		// Number of reads with single hit in the sample.
	int totalNumOfRemovedReads;		// Total number
	int totalNumOfReads;			// Total numbe of the sample reads (This includes reads mapped and not mapped)
	int numOfHitsForSharedReads;	// Total number of hits for all shared reads. This value is for heuristic trimming of the graph if the number of shared hits larger than MAX_NUMBER_OF_ARCS.
	double 	averageReadLength;		// The average length of the reads
	int chunkLength;
	ReadInfo(){
		numOfMappedReads=0;
		numOfSingleHitReads=0;
		totalNumOfRemovedReads=0;
		totalNumOfReads=0;
		numOfHitsForSharedReads=0;
		averageReadLength=0;
		chunkLength=0;
	}
};

struct MCFResult{
	vector<GenomeInfo> 	genomeInfoVector;
	ReadInfo			readInfoSt;
	string 				lgfFile;
	string 				ncbiRefGenomesFile;
	int numOfGenomes;
	int chunkSize;
	int maxScore;
	int minScore;
	MCFResult(string lgfFile, string ncbiRefGenomesFile){
		this->lgfFile=lgfFile;
		this->ncbiRefGenomesFile=ncbiRefGenomesFile;
		numOfGenomes=0;
		chunkSize=-2;
		maxScore=minScore=0;
	}
	void copyResult(MCFResult& mcfResultSt){
		this->genomeInfoVector.insert(genomeInfoVector.end(),mcfResultSt.genomeInfoVector.begin(),mcfResultSt.genomeInfoVector.end());
		this->readInfoSt.averageReadLength=mcfResultSt.readInfoSt.averageReadLength;
	}
	void initGenomes(){
		for(vector<GenomeInfo>::iterator itr=genomeInfoVector.begin();itr!=genomeInfoVector.end();itr++)
		{
			GenomeInfo& gInfo=*itr;
			gInfo.initGenome();
		}
	}

	/*
	 * The function returns a map of the NCBI genomes.
	 * @return: a map(key, value) contains the NCBI genomes. key: the genome name. value: a pair(genomeId, genomeLength). genomeId is a sequence (0,1,...)
	 */
	map<string, pair<int,int> > getNCBIRefMap()
	{
		map<string, pair<int,int> > ncbiMap;
		ifstream ncbiRefFile((char*)ncbiRefGenomesFile.c_str());
		string line;
		vector<string> splittedLine;
		int genomeId=0;
		while (getline(ncbiRefFile, line)){
			split(splittedLine, line, "\t");
			pair<int,int> genomeIdLengthPair=make_pair(genomeId, atoi(splittedLine[1].c_str()));
			if (ncbiMap.find(splittedLine[0]) != ncbiMap.end())
			{
				cout << "duplicated " << genomeId << ":" << splittedLine[0] << " " << splittedLine[1] << endl;
			}
			ncbiMap[splittedLine[0]]=genomeIdLengthPair;
			genomeId++;
		}
		ncbiRefFile.close();
		return ncbiMap;
	}

/*
 * This functions is called only once at the beginning to initialise the structs with the genomes exists in the NCBI genomes dictionary.
 */
void initGenomes(
			ListDigraph& 					g,
			ListDigraph::NodeMap<int>& 		inWhichGenome,
			ListDigraph::NodeMap<string>& 	nodeLabel)
	{
		mcfLogger.log("******************************initGenomes Start********************");
		map<string, pair<int,int> > ncbiRefMap= getNCBIRefMap(); // returns a map of ncbi reference genomes
		// Creates the structs for saving the genomes information (including genome length, chunks, abundance...etc)
		int mapSize=ncbiRefMap.size();
		for(int i=0;i<mapSize;i++){
			GenomeInfo gInfo;
			genomeInfoVector.push_back(gInfo);
		}
		mcfLogger.log("NCBI Map Size:"+getStringValue(mapSize));
		// setting the genome info
		for(map<string, pair<int,int> >::iterator itr=ncbiRefMap.begin();itr!=ncbiRefMap.end();itr++)
		{
			GenomeInfo& gInfo=genomeInfoVector[itr->second.first];
			gInfo.genomeName=itr->first; 
			gInfo.genomeId=itr->second.first;
			gInfo.genomeLength=itr->second.second;
		}
		//int chunkNum;
		// If there is a hit for the genome in the lgf file, create the chunks information for the genome. Otherwise assign -1 to the genome Id. (-1)
		//means the genomes doesn't exist in the mapping file.
		// It adds a chunk to gInfo for each chunk exists in the lgf file.
		for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
		{
			if (inWhichGenome[node] > -1)
			{
				GenomeInfo& gInfo=genomeInfoVector[inWhichGenome[node]];
				if(gInfo.numOfChunks>0)
					continue;
				gInfo.numOfChunks=(gInfo.genomeLength/chunkSize)+1; // Adding 1 chunk, for rounding error
				for(int i=0;i<gInfo.numOfChunks;i++){
					//chunkNum=GenomeInfo::getChunkNum(nodeLabel[node]); // GenomeId_ChunkNum
					gInfo.optimumAbundanceSt.chunksAbundanceVector.push_back(0);
					gInfo.optimumAbundanceSt.chunksMinAbundanceVector.push_back(0);
				}
			}
		}
		for(vector<GenomeInfo >::iterator genomesIterator=genomeInfoVector.begin();genomesIterator!=genomeInfoVector.end();genomesIterator++)
		{
			GenomeInfo& gInfo=*genomesIterator;
			if(gInfo.numOfChunks==0){ // If the number of chunks is zero, means the genomes doesn't exist in the lfg file and is excluded from abundance calculation
				gInfo.genomeId=-1;
			}else{
				gInfo.hasHits=true;
			}
		}
		mcfLogger.log("******************************initGenomes End**********************");
		mcfLogger.log("*******************************************************************");
	}

	/*
	 * A function for calculating the expected chunk abundance for the genomes.
	 */
	void calculateExpectedAbundanceInChunk(){
		for(vector<GenomeInfo >::iterator genomesIterator=genomeInfoVector.begin();genomesIterator!=genomeInfoVector.end();genomesIterator++){
			GenomeInfo& gInfo=*genomesIterator;
			if(gInfo.genomeId!=-1 && gInfo.numOfChunks!=0){
				gInfo.calcExpectedAbundanceInChunk();
	//			cout<<"Genome: "<<gInfo.genomeName<<" Optimum Abundance: "<<gInfo.optimumAbundance<<" ExpectedChunkAbundance: "<<gInfo.expectedChunkAbudance<<endl;
			}
		}
	}

	void eraseGenome(set<int>& genomesToErase){
		for(set<int>::iterator itr=genomesToErase.begin();itr!=genomesToErase.end();itr++){
			genomeInfoVector[*itr].genomeId=-1;
		}
	}

	void updateAbundance(
			vector<ListDigraph::Arc> 				arcsToErase,
			ListDigraph& 							g,
			ListDigraph::NodeMap<int>& 				inWhichGenome,
			ListDigraph::NodeMap<string>& 			nodeLabel
		)
	{
		initAbundance(g,inWhichGenome,nodeLabel);
		for (vector<ListDigraph::Arc>::iterator arcItr = arcsToErase.begin(); arcItr != arcsToErase.end(); ++arcItr)
		{
	//		cout<<"Arc= "<<ListDigraph::id(*arcItr)<<endl;
			ListDigraph::Node targetNode=g.target(*arcItr);
			GenomeInfo& gInfo=genomeInfoVector[inWhichGenome[targetNode]];
			gInfo.optimumAbundanceSt.genomeAbundance--;
			gInfo.optimumAbundanceSt.chunksAbundanceVector[gInfo.getChunkNum(nodeLabel[targetNode])]--;
		}
		calculateExpectedAbundanceInChunk();
	}
	/*
	 * This function used for initialising the abundance of genomes. The following are the steps:
	 * 1) Initialise the abundance for all genomes to zero. Same for the fields: numOfMappedReads, numOfHitsForSharedReads. Check why the field numOfHitsForChunksSharedReads is not initialised to zero (Is it a bug?)
	 * 2) Set the genomes abundance and chunk abundance.
	 * 3) Calculate the expected chunk abundance for each genome.
	 *
	 */
	void initAbundance(
			ListDigraph& 							g, // will be updated
			ListDigraph::NodeMap<int>& 				inWhichGenome,
			ListDigraph::NodeMap<string>& 			nodeLabel
		)
	{
		mcfLogger.log("***********************initAbundance start**************************");
		mcfLogger.log("********************************************************************");
//		set<int> readGenomesSet;
//		set<int> genomesInternalSharedReadsSet;
		int chunkNum=0;
		readInfoSt.numOfMappedReads=0;			// Setting mapped reads info to zero
		readInfoSt.numOfHitsForSharedReads=0;
		initGenomes();// reinitialise the abundance of the genomes (Setting all abundances to zero)
		for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
		{
			if (inWhichGenome[node] == -1)
			{
				readInfoSt.numOfMappedReads++;
				if(countOutArcs(g,node) == 1)
				{
					ListDigraph::OutArcIt arc(g,node);
					GenomeInfo& gInfo=genomeInfoVector[inWhichGenome[g.target(arc)]];
					gInfo.optimumAbundanceSt.genomeAbundance++;
					gInfo.optimumAbundanceSt.genomeMinAbundance++; // If there is a single hit only, increase the genome minimum abundance.
					chunkNum=GenomeInfo::getChunkNum(nodeLabel[g.target(arc)]);
					gInfo.optimumAbundanceSt.chunksAbundanceVector[chunkNum]++;
					gInfo.optimumAbundanceSt.chunksMinAbundanceVector[chunkNum]++;
				}
				// If the read is mapped to more than one genome, increase the abundance of all genomes the read mapped to.
				// If the read is mapped to more than one chunk in the same genome, increase the abundance of the genome only by one. And increase the abundance of the first chunk found ignoring the other chunks.
				else
				{
					map<int,vector<int> > genomesSelectedChunks;
					for (ListDigraph::OutArcIt arc(g,node); arc != INVALID; ++arc)
					{
						GenomeInfo& gInfo=genomeInfoVector[inWhichGenome[g.target(arc)]];
						if(genomesSelectedChunks.find(gInfo.genomeId)==genomesSelectedChunks.end()){
							vector<int> chunksVector;
							genomesSelectedChunks.insert(pair<int, vector<int> >(gInfo.genomeId,chunksVector));
							gInfo.optimumAbundanceSt.genomeAbundance++;
						}
						chunkNum=GenomeInfo::getChunkNum(nodeLabel[g.target(arc)]);
						genomesSelectedChunks.find(gInfo.genomeId)->second.push_back(chunkNum);
//						if(!readGenomesSet.count(gInfo.genomeId)){
//							gInfo.optimumAbundanceSt.genomeAbundance++;
//							chunkNum=GenomeInfo::getChunkNum(nodeLabel[g.target(arc)]);
//							gInfo.optimumAbundanceSt.chunksAbundanceVector[chunkNum]++;
//							readGenomesSet.insert(gInfo.genomeId);// Reads will be duplicated in different places (No Problem, will be filtered out)
//						}else{
//							gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads++; // Increases the number of genome shared arcs if the read is mapped to more than one chunk in the same genome. Check the description of the field.
//							genomesInternalSharedReadsSet.insert(gInfo.genomeId);
//						}
					}
					for(map<int,vector<int> >::iterator genomesItr=genomesSelectedChunks.begin();genomesItr!=genomesSelectedChunks.end();genomesItr++){
						vector<int>& chunksVector=genomesItr->second;
						GenomeInfo& gInfo=genomeInfoVector[genomesItr->first];
						if(chunksVector.size()>1){
							gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads+=chunksVector.size(); // Increases the number of genome shared arcs if the read is mapped to more than one chunk in the same genome. Check the description of the field.
						}
						chunkNum=chunksVector[rand()%chunksVector.size()];
						gInfo.optimumAbundanceSt.chunksAbundanceVector[chunkNum]++;
					}
					if(genomesSelectedChunks.size()>1){
						readInfoSt.numOfHitsForSharedReads+=genomesSelectedChunks.size();// Increases the number of total arcs between shared reads. Check the description of the field.
					}
					genomesSelectedChunks.clear();
					// Add one to numOfHitsForChunksSharedReads for all genomes have a read shared by more than one chunk in the same genome.
//					for(set<int>::iterator itr=genomesInternalSharedReadsSet.begin();itr!=genomesInternalSharedReadsSet.end();itr++){
//						genomeInfoVector[*itr].optimumAbundanceSt.numOfHitsForChunksSharedReads++;
//					}
//					genomesInternalSharedReadsSet.clear();
//					readGenomesSet.clear();
//					readInfoSt.numOfHitsForSharedReads+=countOutArcs(g,node); // Increases the number of total arcs between shared reads. Check the description of the field.
				}
			}
		}
		mcfLogger.log("Total Number of all shared reads= "+getStringValue(readInfoSt.numOfHitsForSharedReads));
		calculateExpectedAbundanceInChunk(); // Calculating the expected chunk abundance of the genoms.
		mcfLogger.log("**************************initAbundance End*************************");
		mcfLogger.log("********************************************************************");
	}
	void copyLocalAbundance()
	{
		for(vector<GenomeInfo >::iterator genomesIterator=genomeInfoVector.begin();genomesIterator!=genomeInfoVector.end();genomesIterator++)
		{
			GenomeInfo& gInfo=*genomesIterator;
			gInfo.copyLocalAbundance();
		}
	}
	int getNumberOfGenomes(){
		int numOfGenomes=0;
		for (vector<GenomeInfo>::iterator itr= genomeInfoVector.begin(); itr != genomeInfoVector.end(); ++itr)
		{
			GenomeInfo& gInfo=*itr;
			if(gInfo.genomeId!=-1)
				numOfGenomes++;
		}
		return numOfGenomes;
	}

};



#endif /* MCFRESULT_H_ */
