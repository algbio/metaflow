
#ifndef __MCFmapper__MCFflowsolver__
#define __MCFmapper__MCFflowsolver__

#include "MCFHeaders.h"
#include "MCFResult.h"
#include "MCFUtils.h"

struct GraphErasedData{
	vector<ListDigraph::Node> nodesToErase;
	vector<ListDigraph::Arc> arcsToErase;
	void clear(){
		nodesToErase.clear();
		arcsToErase.clear();
	}
};

struct FlowNetworkConfig{
	int genomeId;
	bool useUnknownNode;
	bool removeSingleHitReads;
	FlowNetworkConfig(int genomeId, bool useUnkownNode, int removeSingleHitReads){
		this->genomeId=genomeId;
		this->useUnknownNode=useUnkownNode;
		this->removeSingleHitReads=removeSingleHitReads;
	}
};

int64_t solveWithGivenAbundances(
	MCFResult&	 					mcfResultSt,
	ListDigraph& 					g,
	ListDigraph::NodeMap<int>& 		inWhichGenome,
	ListDigraph::NodeMap<string>& 	nodeLabel,
	ListDigraph::ArcMap<int>& 		arcWeight,
	GraphErasedData& 				graphErasedDataSt,
	FlowNetworkConfig&				flowNetworkConfigSt
);

#endif /* defined(__MCFmapper__MCFflowsolver__) */
