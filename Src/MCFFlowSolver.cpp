#include "MCFGenetic.h"
//#define NON_UNIFORM_PENALITY 100
#include "MCFLogger.h"
extern MCFConfig configSt;
extern MCFLogger mcfLogger;
int64_t updateGenomesOptimumAbundance(
	ListDigraph& 								g,
	ListDigraph& 								flowNetwork,
	ListDigraph::NodeMap<int>& 					inWhichGenome,
	ListDigraph::NodeMap<string>& 				nodeLabel,
	ListDigraph::NodeMap<ListDigraph::Node>& 	graphNodetoFlowNode,
	ListDigraph::NodeMap<ListDigraph::Node>& 	flowNodetoGraphNode,
	ListDigraph::ArcMap<int64_t>& 				flowMap,
	ListDigraph::ArcMap<int64_t>& 				costMap,
	ListDigraph::Node							unkownNode,
	ListDigraph::Node							expectedAbundanceNode,
	MCFResult& 									mcfResultSt,
	GraphErasedData& 							graphErasedDataSt,
	FlowNetworkConfig&							flowNetworkConfigSt
)
{
	int64_t cost=0;
//	if(flowNetworkConfigSt.genomeId==-1)
//		mcfResultSt.initGenomes();
//	else{
//		mcfResultSt.genomeInfoVector[flowNetworkConfigSt.genomeId].initGenome();
//	}
	graphErasedDataSt.clear();
	int whichGenome;
	for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
	{
		if(countInArcs(g,node)==0 && countOutArcs(g,node)==0)
			continue;
		whichGenome=inWhichGenome[node];
		if(flowNetworkConfigSt.genomeId!=-1){
			if (whichGenome > -1 && flowNetworkConfigSt.genomeId!=whichGenome)
				continue;
			if(whichGenome== -1){
				bool exist=false;
				for (ListDigraph::OutArcIt arc(g,node); arc != INVALID; ++arc){
					if(inWhichGenome[g.target(arc)] ==flowNetworkConfigSt.genomeId){
						exist=true;
						break;
					}
				}
				if(!exist)
					continue;
			}
		}

		ListDigraph::Node flNode = graphNodetoFlowNode[node];
		// If it's a read node
		if (whichGenome == -1)
		{
			if(flowNetworkConfigSt.removeSingleHitReads && countOutArcs(g,node) == 1){
				continue;
			}
			for (ListDigraph::OutArcIt arc(flowNetwork,flNode); arc != INVALID; ++arc)
			{ // go through all in-arcs of flNode
				ListDigraph::Node flowTargetNode = flowNetwork.target(arc) ;
				ListDigraph::Node graphChunkNode= flowNodetoGraphNode[flowTargetNode];
				if(!flowNetworkConfigSt.useUnknownNode){
					if ( (flowMap[arc] == 0) && (flowTargetNode != INVALID) ){
						graphErasedDataSt.arcsToErase.push_back(findArc(g,node, graphChunkNode));
					}
				}else{
//					if ( (flowMap[arc] > 0) && (flowTargetNode != INVALID) && (flowTargetNode==unkownNode)){
//						for (ListDigraph::OutArcIt gArc(g,node); gArc != INVALID; ++gArc){
//							graphErasedDataSt.arcsToErase.push_back(gArc);
//						}
//					}
					if ( (flowMap[arc] ==0) && (flowTargetNode != INVALID) && (flowTargetNode!=unkownNode)){
						graphErasedDataSt.arcsToErase.push_back(findArc(g,node, graphChunkNode));
					}
				}
			}
		}
	}
    return cost;
}

int64_t solveWithGivenAbundances(
	MCFResult& 								mcfResultSt,
	ListDigraph& 							g, 
	ListDigraph::NodeMap<int>& 				inWhichGenome,
	ListDigraph::NodeMap<string>& 			nodeLabel, 
	ListDigraph::ArcMap<int>& 				arcWeight,
	GraphErasedData& 						graphErasedDataSt,
	FlowNetworkConfig&						flowNetworkConfigSt
) 
{
	clock_t tStart=clock();
	/****************************************/
	// 		Creating the flow network
	/****************************************/
	mcfLogger.log("***********************solveWithGivenAbundances Start*****************");
	string species="ALL";
	if(flowNetworkConfigSt.genomeId!=-1)
		species=mcfResultSt.genomeInfoVector[flowNetworkConfigSt.genomeId].genomeName;
	mcfLogger.log("Constructing the flow network for "+species+" ...");
	// Alpha takes values between [0,1]
	double alpha=configSt.getDoubleProperty("ALPHA");
	double beta=1-alpha;
//	int lowReadsLenientCost=beta*configSt.getIntProperty(LOW_READS_LENIENT_COST);
//	int lowReadsRegCost=beta*configSt.getIntProperty(LOW_READS_REGULAR_COST);
//	int extraReadsLenientCost=beta*configSt.getIntProperty(EXTRA_READS_LENIENT_COST);
//	int extraReadsRegCost=beta*configSt.getIntProperty(EXTRA_READS_REGULAR_COST);
//	int unknownCost=configSt.getIntProperty(UNKWON_COST);
	int lowReadsLenientCost=alpha*mcfResultSt.minScore;
	int lowReadsRegCost=alpha*mcfResultSt.maxScore;
	int extraReadsLenientCost=alpha*mcfResultSt.minScore;
	int extraReadsRegCost=alpha*mcfResultSt.maxScore;
	int unknownCost=0.75*mcfResultSt.maxScore+0.25*mcfResultSt.minScore;
	string useUnknownNode="True";
	if(!flowNetworkConfigSt.useUnknownNode)
		useUnknownNode="False";
	mcfLogger.log("UseUnkownNode= "+useUnknownNode+"\t lowReadsLenientCost="+getStringValue(lowReadsLenientCost)+
			"\tlowReadsRegCost="+getStringValue(lowReadsRegCost)+"\tunknownCost="+getStringValue(unknownCost));
	ListDigraph flowNetwork;
    ListDigraph::ArcMap<int64_t> lowerMap(flowNetwork), upperMap(flowNetwork);
    ListDigraph::ArcMap<int64_t> costMap(flowNetwork);
    ListDigraph::NodeMap<int64_t> supplyMap(flowNetwork);
    
    //correspondence between nodes in flowNetwork and nodes in g
    ListDigraph::NodeMap<ListDigraph::Node> flowNodetoGraphNode(flowNetwork);
    ListDigraph::NodeMap<ListDigraph::Node> graphNodetoFlowNode(g);

	// adding star source and sink
	ListDigraph::Node s_star = flowNetwork.addNode(); 
	ListDigraph::Node t_star = flowNetwork.addNode(); 
	ListDigraph::Arc a = flowNetwork.addArc(s_star,t_star); /**/ lowerMap.set(a,0); upperMap.set(a,int64_t_MAX); costMap.set(a,0);
	

	//adding expected total abundance node
	int64_t totalExpectedAbundance=0;
	ListDigraph::Node expectedAbundanceNode = flowNetwork.addNode();


	//adding 'unknown' chunk
	ListDigraph::Node unKnownNode;
	if(flowNetworkConfigSt.useUnknownNode){
		unKnownNode= flowNetwork.addNode();
		a = flowNetwork.addArc(unKnownNode,t_star); /**/ lowerMap.set(a,0); upperMap.set(a,mcfResultSt.readInfoSt.numOfMappedReads); costMap.set(a,0);
	}
	GenomeInfo gInfo;
	int whichGenome;

	for (ListDigraph::NodeIt node(g); node != INVALID; ++node) 
	{
		if(countInArcs(g,node)==0 && countOutArcs(g,node)==0)
			continue;
		whichGenome=inWhichGenome[node];
		if (whichGenome < 0) // a read node
		{
			if(flowNetworkConfigSt.genomeId!=-1){ // Construct the network for only one genome
				bool exist=false;
				for (ListDigraph::OutArcIt arc(g,node); arc != INVALID; ++arc){
					if(inWhichGenome[g.target(arc)] == mcfResultSt.genomeInfoVector[flowNetworkConfigSt.genomeId].genomeId){
						exist=true;
						break;
					}
				}
				if(!exist) // There is no arcs from this read to the genome we calculate the abundance for.
					continue;
			}
			if(flowNetworkConfigSt.removeSingleHitReads && countOutArcs(g,node) == 1){ // Don't include reads with single hits in the graph (to minimise the graph size)
				continue;
			}
			// add readNode to flowNetwork
			ListDigraph::Node readNode = flowNetwork.addNode(); /**/ supplyMap.set(readNode,0);
			flowNodetoGraphNode.set(readNode,node);
			graphNodetoFlowNode.set(node,readNode);

			// Add an arc from s_star to the read node and push 1
			a = flowNetwork.addArc(s_star,readNode); /**/ lowerMap.set(a,1); upperMap.set(a,1); costMap.set(a,0);

			if(flowNetworkConfigSt.useUnknownNode){
				// add an arc from the readNode to the unknownNode
				a = flowNetwork.addArc(readNode, unKnownNode); /**/ lowerMap.set(a,0); upperMap.set(a,1);costMap.set(a,unknownCost);
			}
		} else // a genome chunk node
		{
			if(flowNetworkConfigSt.genomeId!=-1 && whichGenome!=flowNetworkConfigSt.genomeId){
				continue; // construct the network for only one genome
			}

			// add chunkNode to flowNetwork
			ListDigraph::Node chunkNode = flowNetwork.addNode(); /**/ supplyMap.set(chunkNode,0);
			flowNodetoGraphNode.set(chunkNode,node);
			graphNodetoFlowNode.set(node,chunkNode);
			gInfo=mcfResultSt.genomeInfoVector[whichGenome];
			string chunkName=nodeLabel[node];
			int64_t expectedAbundanceInChunk = gInfo.getChunkExpecetedAbundance(chunkName, flowNetworkConfigSt.removeSingleHitReads);
			totalExpectedAbundance+=expectedAbundanceInChunk;
//			cout<<"Genome=" << genomeAbundances[inWhichGenome[node]]<<" Chunks="<<chunksPerGenome[inWhichGenome[node]]<<" expected="<<expectedAbundanceInChunk<<endl;


			// pushing the expectedAbundance to the sink with zero cost
			a = flowNetwork.addArc(chunkNode,t_star); /**/ lowerMap.set(a,expectedAbundanceInChunk); upperMap.set(a,expectedAbundanceInChunk); costMap.set(a,0);

			if(configSt.getboolProperty(SPLIT_ARCS)){
				int numOfLowPenaltyReads=calcNumOfReadsWithLowPenalty(gInfo,chunkName, flowNetworkConfigSt.removeSingleHitReads);
				//adding arcs to penalize for coverage lower than required "expected" abundance
				a = flowNetwork.addArc(expectedAbundanceNode,chunkNode); /**/ lowerMap.set(a,0); upperMap.set(a,numOfLowPenaltyReads); costMap.set(a,lowReadsLenientCost);
				a = flowNetwork.addArc(expectedAbundanceNode,chunkNode); /**/ lowerMap.set(a,0); upperMap.set(a,expectedAbundanceInChunk-numOfLowPenaltyReads); costMap.set(a,lowReadsRegCost);
				//adding arcs to penalize for mappings exceeding the expected abundance
				a = flowNetwork.addArc(chunkNode,t_star); /**/ lowerMap.set(a,0); upperMap.set(a,numOfLowPenaltyReads); costMap.set(a,extraReadsLenientCost);
				a = flowNetwork.addArc(chunkNode,t_star); /**/ lowerMap.set(a,0); upperMap.set(a,expectedAbundanceInChunk+countInArcs(g,node)-numOfLowPenaltyReads); costMap.set(a,extraReadsRegCost);
			}
			else{
				a = flowNetwork.addArc(expectedAbundanceNode,chunkNode); /**/ lowerMap.set(a,0); upperMap.set(a,expectedAbundanceInChunk); costMap.set(a,lowReadsRegCost);
				a = flowNetwork.addArc(chunkNode,t_star); /**/ lowerMap.set(a,0); upperMap.set(a,expectedAbundanceInChunk); costMap.set(a,extraReadsRegCost);
			}
		}
	}

	// adding the arcs from reads to chunks
	for (ListDigraph::ArcIt arc(g); arc != INVALID; ++arc)
	{
		if(flowNetworkConfigSt.genomeId!=-1 && inWhichGenome[g.target(arc)]!=mcfResultSt.genomeInfoVector[flowNetworkConfigSt.genomeId].genomeId){
			continue; // include only the arcs for a specific genome if genomeId!=-1
		}
		if(flowNetworkConfigSt.removeSingleHitReads && countOutArcs(g,g.source(arc)) == 1){
			continue; // exclude single hits reads
		}
		// It's minimisation solver, and blast score is reversed (the lower, the better).
		a = flowNetwork.addArc(graphNodetoFlowNode[g.source(arc)],graphNodetoFlowNode[g.target(arc)]); /**/ lowerMap.set(a,0); upperMap.set(a,1); costMap.set(a,beta*arcWeight[arc]);
	}

	int supply_s_star=totalExpectedAbundance+mcfResultSt.readInfoSt.numOfMappedReads;
	mcfLogger.log( "supply_s_star = "+getStringValue(supply_s_star));

	a = flowNetwork.addArc(s_star,expectedAbundanceNode); /**/ lowerMap.set(a,totalExpectedAbundance); upperMap.set(a,totalExpectedAbundance); costMap.set(a,0);
	a = flowNetwork.addArc(expectedAbundanceNode,t_star); /**/ lowerMap.set(a,0); upperMap.set(a,totalExpectedAbundance); costMap.set(a,0);

	supplyMap.set(s_star,supply_s_star);
	supplyMap.set(t_star,-supply_s_star);

	/****************************************/
	// 		Running the min-cost flow engine
	/****************************************/

	
	
	// Network Simplex Algorithm
    NetworkSimplex<ListDigraph, int64_t> ns(flowNetwork);
    ns.lowerMap(lowerMap).upperMap(upperMap).costMap(costMap).supplyMap(supplyMap);
    printRunningTime(tStart, "Constructing the flow network took: ");
    tStart=clock();
    printGraphInfo(flowNetwork);
    mcfLogger.log("Running the flow engine...");
    NetworkSimplex<ListDigraph, int64_t>::ProblemType res = ns.run(NetworkSimplex<ListDigraph, int64_t>:: CANDIDATE_LIST); //ALTERING_LIST  FIRST_ELIGIBLE BEST_ELIGIBLE BLOCK_SEARCH CANDIDATE_LIST
    if (res != NetworkSimplex<ListDigraph, int64_t>::OPTIMAL)
    {
    	mcfLogger.log("ERROR: flow not found");
    }
    ListDigraph::ArcMap<int64_t> flowMap(flowNetwork); ns.flowMap(flowMap);    
    printRunningTime(tStart, "Running the flow engine & finding the minimum cost flow took: ");

    /****************************************/
	// 		Extracting the original cost from the flow
	/****************************************/
    tStart=clock();
    int64_t cost = ns.totalCost();
    mcfLogger.log("Objective function: " +getStringValue((int)cost));
    updateGenomesOptimumAbundance(
    			g,
    			flowNetwork,
    			inWhichGenome,
    			nodeLabel,
				graphNodetoFlowNode,
				flowNodetoGraphNode,
				flowMap,
				costMap,
				unKnownNode,
				expectedAbundanceNode,
				mcfResultSt,
				graphErasedDataSt,
				flowNetworkConfigSt
				);
    mcfLogger.log("**************************solveWithGivenAbundances End*****************");
    mcfLogger.log("***********************************************************************");
    return cost;
}
