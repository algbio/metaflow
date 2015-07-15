
#ifndef __MCFmapper__MCFutils__
#define __MCFmapper__MCFutils__

#include "MCFHeaders.h"
#include "MCFResult.h"
#include "MCFUtils_Temp.h"

void initDigraph(
	MCFResult& 						mcfResultSt,
	ListDigraph& 					g,
	ListDigraph::NodeMap<int>& 		inWhichGenome,
	ListDigraph::NodeMap<string>& 	nodeLabel,
	ListDigraph::ArcMap<int>& 		arcWeight
);
void readDigraph(
	MCFResult& 						mcfResultSt,
	ListDigraph& 					g,
	ListDigraph::NodeMap<int>& 		inWhichGenome,
	ListDigraph::NodeMap<string>& 	nodeLabel,
	ListDigraph::ArcMap<int>& 		arcWeight
);


void trimDigraph(
	MCFResult& 								mcfResultSt,// will be updated
	ListDigraph& 							g, // will be updated
	ListDigraph::NodeMap<int>& 				inWhichGenome,
	ListDigraph::NodeMap<string>& 			nodeLabel
);
void eraseArcs(
	vector<ListDigraph::Arc>& arcsToErase,
	ListDigraph& g
);
void eraseArcs(
	ListDigraph& 							g,
	ListDigraph::NodeMap<int>& 				inWhichGenome,
	ListDigraph::NodeMap<string>& 			nodeLabel,
	vector<ListDigraph::Arc>& 				arcsToErase
);
void eraseNodes(
	vector<ListDigraph::Node> nodesToErase,
	ListDigraph& g
);
void eraseNodes(
	set<int>& 					genomesToErase,
	ListDigraph& 				g,
	ListDigraph::NodeMap<int>& 	inWhichGenome
);
void printRunningTime(clock_t& startTime, string msg);
void printNodes(
	ListDigraph& 					g,
	ListDigraph::NodeMap<int>& 		inWhichGenome,
	vector<string>& 				genomeNames,
	ListDigraph::NodeMap<string>& 	nodeLabel
);
void saveChunksResults(MCFResult& mcfResultSt);
void saveChunksResults(MCFResult& mcfResultSt, string fileName);
void saveAbundanceResults(MCFResult& mcfResultSt, string fileName);
void saveAbundanceResults(MCFResult& mcfResultSt);
void saveNaiveBlastAbundanceResults(
	string lgfFile,
	string ncbiRefGenomesFile
);
void printGenomeInfo(vector<GenomeInfo>& genomeInfoVector);
void printGraphInfo(ListDigraph& g);
void printRelativeAbundance(
	MCFResult& 	mcfresult
);

bool isExceedingMaxNumOfArcs(int numOfSharedArcs);
void eraseArcs(
	MCFResult& 						mcfResultSt,
	ListDigraph& 					g,
	ListDigraph::NodeMap<int>& 		inWhichGenome,
	ListDigraph::NodeMap<string>& 	nodeLabel,
	ListDigraph::ArcMap<int>& 		arcWeight,
	int 							genomeId
);
int calcNumOfReadsWithLowPenalty(GenomeInfo& gInfo, string chunkName, bool removeSingleHitReads);

void removeOutOfBoundaryHits(
	MCFResult& mcfResultSt,
	ListDigraph& 					g,
	ListDigraph::NodeMap<int>& 		inWhichGenome,
	ListDigraph::NodeMap<string>& 	nodeLabel
);

#endif /* defined(__MCFmapper__MCFutils__) */
