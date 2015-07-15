
#ifndef __MCFmapper__MCFgenetic__
#define __MCFmapper__MCFgenetic__

#include "MCFFlowSolver.h"
#include "MCFGenomeInfo.h"


vector<int> GA_findAbundances(
	int 								numberOfGenomes, 
	ListDigraph& 						g, 
	ListDigraph::NodeMap<string>& 		nodeLabel, 
	ListDigraph::ArcMap<double>& 			arcWeight,
	ListDigraph::NodeMap<int>& 			inWhichGenome,
	vector<int64_t> 					genomeLengths,		
	vector<int>& 						chunksPerGenome,
	vector<int>&						minimumAbundancesPerGenome,
	ListDigraph::NodeMap<int>&			minimumAbundancesPerChunk	
);
void findMaxAbundance(
	MCFResult& 							mcfResultSt,
	ListDigraph& 						g,
	ListDigraph::NodeMap<int>& 			inWhichGenome,
	ListDigraph::NodeMap<string>& 		nodeLabel,
	ListDigraph::ArcMap<int>&			arcWeight,
	FlowNetworkConfig& 					flowNetworkConfigSt
);
void greedy_findAbundances(
	MCFResult& 							mcfResultSt,
	ListDigraph& 						g,
	ListDigraph::NodeMap<int>& 			inWhichGenome,
	ListDigraph::NodeMap<string>& 		nodeLabel, 
	ListDigraph::ArcMap<int>&			arcWeight
);

bool stopMutation(vector<int> previousAbundance, vector<GenomeInfo>& genomeInfoVector, int64_t previousCost, int64_t currentCost, int nOfMutations, int maxNumOfMutations);
vector<int> copyPreviousAbundance(vector<GenomeInfo>& genomeInfoVector);
#endif /* defined(__MCFmapper__MCFgenetic__) */
