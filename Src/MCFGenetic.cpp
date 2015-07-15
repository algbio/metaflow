#include "MCFGenetic.h"
#include "MCFLogger.h"
#include "MCFUtils_Temp.h"
//////////////////////////////////////////////////////////////////////////
extern MCFConfig configSt;
extern MCFLogger mcfLogger;
void findMaxAbundance(
		MCFResult& 							mcfResultSt,
		ListDigraph& 						g,
		ListDigraph::NodeMap<int>& 			inWhichGenome,
		ListDigraph::NodeMap<string>& 		nodeLabel,
		ListDigraph::ArcMap<int>&			arcWeight,
		FlowNetworkConfig&					flowNetworkConfigSt
	){
	mcfLogger.log("***********************findMaxAbundance Start********************");
	GraphErasedData graphErasedDataSt;

	int chunkNum;
	for(vector<GenomeInfo >::iterator genomesIterator=mcfResultSt.genomeInfoVector.begin();genomesIterator!=mcfResultSt.genomeInfoVector.end();genomesIterator++)
	{
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId!=-1){
			mcfLogger.log("Finding Max Abundance For: "+gInfo.genomeName);
			mcfLogger.log("Before trimming: Number of shared arcs: "+getStringValue(gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads)+"  Max Number Of Arcs= "+configSt.getStrProperty(MAX_NUMBER_OF_ARCS));
			if(isExceedingMaxNumOfArcs(gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads)){
				eraseArcs(mcfResultSt,g,inWhichGenome,nodeLabel,arcWeight,gInfo.genomeId);
			}
			mcfLogger.log("After trimming: Number of shared arcs: "+getStringValue(gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads));

			while(true){
				flowNetworkConfigSt.genomeId=gInfo.genomeId;
				solveWithGivenAbundances(mcfResultSt, g, inWhichGenome, nodeLabel, arcWeight, graphErasedDataSt, flowNetworkConfigSt);
				if(graphErasedDataSt.arcsToErase.size()==0)
					break;

				// Only in the final step, we assign the reads to unkown.
				if(flowNetworkConfigSt.useUnknownNode) {
					for(vector<ListDigraph::Arc>::iterator arcsItr=graphErasedDataSt.arcsToErase.begin(); arcsItr!=graphErasedDataSt.arcsToErase.end();arcsItr++){
						chunkNum=GenomeInfo::getChunkNum(nodeLabel[g.target(*arcsItr)]);
						gInfo.optimumAbundanceSt.genomeAbundance--;
						gInfo.optimumAbundanceSt.chunksAbundanceVector[chunkNum]--;
					}
					gInfo.calcExpectedAbundanceInChunk();
				}
				eraseArcs(graphErasedDataSt.arcsToErase, g);
//				if(flowNetworkConfigSt.useUnknownNode) // Performance improvement: Use the removed arcs instead of looping over the whole flow network.
//					mcfResultSt.initAbundance(g,inWhichGenome, nodeLabel);
			}
		}
	}
	mcfResultSt.initAbundance(g,inWhichGenome, nodeLabel);
	mcfLogger.log("***********************findMaxAbundance End**********************");
	mcfLogger.log("*****************************************************************");
}

void greedy_findAbundances(
	MCFResult& 				mcfResultSt,
	ListDigraph& 						g, 
	ListDigraph::NodeMap<int>& 			inWhichGenome,
	ListDigraph::NodeMap<string>& 		nodeLabel, 
	ListDigraph::ArcMap<int>&			arcWeight
)
{
	int64_t currentCost;

	mcfLogger.log("******************greedy_findAbundances Start********************");
	int64_t previousCost=10000000000000;
	vector<int> previousAbundance=copyPreviousAbundance(mcfResultSt.genomeInfoVector);
	int nOfMutations=0;
	GraphErasedData graphErasedDataSt;
	FlowNetworkConfig flowNetworkConfigSt(-1, false, true); // Don't use Unknown node, Remove single hit reads
	mcfLogger.log("Before trimming: Number of shared arcs: "+getStringValue(mcfResultSt.readInfoSt.numOfHitsForSharedReads)+"  Max Number Of Arcs= "+configSt.getStrProperty(MAX_NUMBER_OF_ARCS));
	if(isExceedingMaxNumOfArcs(mcfResultSt.readInfoSt.numOfHitsForSharedReads)){
		eraseArcs(mcfResultSt,g,inWhichGenome,nodeLabel,arcWeight,-1);
	}
	mcfLogger.log("After trimming: Number of shared arcs: "+getStringValue(mcfResultSt.readInfoSt.numOfHitsForSharedReads));
	clock_t startTime=clock();
	double maxRunningTime=configSt.getDoubleProperty(MAX_RUNNING_TIME);
	while(true){
			mcfLogger.log("Total number of reads: "+getStringValue(mcfResultSt.readInfoSt.numOfMappedReads));
			currentCost = solveWithGivenAbundances(mcfResultSt, g, inWhichGenome, nodeLabel, arcWeight, graphErasedDataSt, flowNetworkConfigSt);
			printGenomeInfo(mcfResultSt.genomeInfoVector);
			if(stopMutation(previousAbundance,mcfResultSt.genomeInfoVector,previousCost,currentCost, nOfMutations, configSt.getIntProperty(MAX_NUMBER_OF_MUTATION_LOOPS))){
				mcfLogger.log("Breaking out of the loop");
				break;
			}
			if(((double)(clock()-startTime)/CLOCKS_PER_SEC/60 )> maxRunningTime){
				mcfLogger.log("Exceeding max running time. Breaking out of the loop");
				break;
			}
			mcfResultSt.updateAbundance(graphErasedDataSt.arcsToErase,g,inWhichGenome, nodeLabel);
//			eraseNodes(graphErasedDataSt.nodesToErase,g);
			trimDigraph(mcfResultSt, g, inWhichGenome, nodeLabel);
//			mcfResultSt.initAbundance(g,inWhichGenome, nodeLabel);
			previousCost=currentCost;
			previousAbundance=copyPreviousAbundance(mcfResultSt.genomeInfoVector);
			nOfMutations++;
	}
	eraseArcs(graphErasedDataSt.arcsToErase,g);
	mcfResultSt.initAbundance(g,inWhichGenome, nodeLabel);
	mcfLogger.log("*********************greedy_findAbundances End********************");
	mcfLogger.log("********************************************************************");
}

bool stopMutation(vector<int> previousAbundance, vector<GenomeInfo>& genomeInfoVector, int64_t previousCost, int64_t currentCost, int nOfMutations, int maxNumOfMutations){
	if(nOfMutations>maxNumOfMutations)
		return true;
	int64_t costDifference=absolute(previousCost-currentCost);
	if(costDifference<configSt.getIntProperty(MAX_COST_DIFFERENCE))
		return true;
	int64_t readsDifference=0;
	GenomeInfo gInfo;
	int genomeId=0;
	for(vector<GenomeInfo >::iterator genomesIterator=genomeInfoVector.begin();genomesIterator!=genomeInfoVector.end();genomesIterator++){
		gInfo=*genomesIterator;
		readsDifference+=abs(gInfo.optimumAbundanceSt.genomeAbundance-previousAbundance[genomeId]);
		genomeId++;
	}
	if(readsDifference<configSt.getIntProperty(MAX_READS_DIFFERENCE))
		return true;
	return false;
}

vector<int> copyPreviousAbundance(vector<GenomeInfo>& genomeInfoVector){
	vector<int> previousAbundance;
	GenomeInfo gInfo;
	for(vector<GenomeInfo >::iterator genomesIterator=genomeInfoVector.begin();genomesIterator!=genomeInfoVector.end();genomesIterator++){
		previousAbundance.push_back(gInfo.optimumAbundanceSt.genomeAbundance);
	}
	return previousAbundance;
}
