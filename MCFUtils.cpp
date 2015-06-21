#include "MCFUtils.h"
#include "MCFLogger.h"

extern MCFLogger mcfLogger;
extern MCFConfig configSt;

/*
 * The function does the following steps:
 * 1) Fills the lemon objects with the mapping details in the lgf file.
 * 2) Initialises the mcfResultSt genomes that have Blast hits in the lemon graph.
 * 3) Initialises the mcfResultSt abundance information for these genomes.
 * @param mcfResultSt: Struct holding the genomes and their abundance information.
 * @param g,inWhichGenome, nodeLabel, arcWeight: Lemon library objects.
 *
 */
void initDigraph(
		MCFResult& 						mcfResultSt,
		ListDigraph& 					g,
		ListDigraph::NodeMap<int>& 		inWhichGenome,
		ListDigraph::NodeMap<string>& 	nodeLabel,
		ListDigraph::ArcMap<int>& 		arcWeight
	){
	readDigraph(mcfResultSt,g,inWhichGenome,nodeLabel, arcWeight);
	mcfResultSt.initGenomes(g, inWhichGenome, nodeLabel);
	removeOutOfBoundaryHits(mcfResultSt,g,inWhichGenome,nodeLabel);// Usless function ?
	mcfResultSt.initAbundance(g,inWhichGenome, nodeLabel);
}
void saveNaiveBlastAbundanceResults(
		string lgfFile,
		string ncbiRefGenomesFile
	){
	mcfLogger.log("****************************saveNaiveBlastAbundanceResults Start*************************");
	ListDigraph g;   //
	ListDigraph::ArcMap<int> arcWeight(g); // A map holds the acrs weights (Blast score)
	ListDigraph::NodeMap<string> nodeLabel(g); // A map holds the nodes label (represents the genomes chunks) in the form (GenomeId_ChunkNumber)
	ListDigraph::NodeMap<int> inWhichGenome(g); // a map holds the genomes ids in the form (0,1,2,...)arcMap(arcWeight, arcWeightCopy).run();
	MCFResult naiveBlastResultSt(lgfFile,ncbiRefGenomesFile);
	readDigraph(naiveBlastResultSt,g,inWhichGenome,nodeLabel, arcWeight);
//	naiveBlastResultSt.chunkSize
	vector<ListDigraph::Arc> arcsToErase;
	int arcToKeep=0;
	for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
	{
		if (inWhichGenome[node] == -1)
		{

			for (ListDigraph::OutArcIt arc(g,node); arc != INVALID; ++arc){
				arcsToErase.push_back(arc);
			}
			arcToKeep=rand()%arcsToErase.size();
			arcsToErase.erase(arcsToErase.begin()+arcToKeep);
			eraseArcs(arcsToErase, g);
			arcsToErase.clear();
		}
	}
	naiveBlastResultSt.initGenomes(g, inWhichGenome, nodeLabel);
	naiveBlastResultSt.initAbundance(g,inWhichGenome, nodeLabel);

	double totalAbundance=0;
	for(vector<GenomeInfo >::iterator genomesIterator=naiveBlastResultSt.genomeInfoVector.begin();genomesIterator!=naiveBlastResultSt.genomeInfoVector.end();genomesIterator++){
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId!=-1){
			gInfo.optimumAbundanceSt.absoluteAbundance=gInfo.optimumAbundanceSt.genomeAbundance*naiveBlastResultSt.readInfoSt.averageReadLength/(double)gInfo.genomeLength;
			totalAbundance+=gInfo.optimumAbundanceSt.absoluteAbundance;
		}
	}
	for(vector<GenomeInfo >::iterator genomesIterator=naiveBlastResultSt.genomeInfoVector.begin();genomesIterator!=naiveBlastResultSt.genomeInfoVector.end();genomesIterator++){
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId!=-1){
			if((gInfo.optimumAbundanceSt.absoluteAbundance/totalAbundance)<configSt.getDoubleProperty(NAIVE_BLAST_MIN_ABUNDANCE))
				gInfo.genomeId=-1;
		}
	}
	saveAbundanceResults(naiveBlastResultSt,naiveBlastResultSt.lgfFile+".NaiveBlast");
	mcfLogger.log("****************************saveNaiveBlastAbundanceResults End*************************");
}
/*
 * This function removes the arcs that has low BLAST score (Hits that their score is less than the best hit score).
 * @param g,inWhichGenome, arcWeight : Lemon objects
 *
 */
void removeLowScoreArcs(
		ListDigraph& 					g,
		ListDigraph::NodeMap<int>& 		inWhichGenome,
		ListDigraph::ArcMap<int>& 		arcWeight
	){
	int max_score_diff=configSt.getIntProperty(MAX_SCORE_DIFFERENCE);
	if(max_score_diff > -1){ // If negative value do nothing
		vector<ListDigraph::Arc> arcsToErase;
		for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
		{
			double minWeight=10000000;
			if (inWhichGenome[node] == -1)
			{
				// Finding the max blast score.
				for (ListDigraph::OutArcIt arc(g,node); arc != INVALID; ++arc){
					if(minWeight<arcWeight[arc])
						minWeight=arcWeight[arc];
				}
				minWeight+=max_score_diff;
				// adding the low score arcs to arcsToErase vector.
				for (ListDigraph::OutArcIt arc(g,node); arc != INVALID; ++arc){
					if(minWeight<arcWeight[arc])
						arcsToErase.push_back(arc);
				}
			}
		}
		// remove low score arcs
		eraseArcs(arcsToErase,g);
	}
}
/*
 * Initilasing the lemon library objects from lgf file. This includes reading the lgf file and removing acrs with low blast score.
 * It also updates the mcfResultSt with the number of mapped read and average readLength.
 * @param mcfResultSt: holds the graph file name.
 * @param g,inWhichGenome,inWhichGenome,arcWeight: Lemon library objects
 */
void readDigraph(
		MCFResult& 						mcfResultSt,
		ListDigraph& 					g,
		ListDigraph::NodeMap<int>& 		inWhichGenome,
		ListDigraph::NodeMap<string>& 	nodeLabel,
		ListDigraph::ArcMap<int>& 		arcWeight
	){
	try {
			digraphReader(g, mcfResultSt.lgfFile). // read the graph into g
				arcMap("weight", arcWeight).			 // read the 'weight' arc map into edgeWeight
				nodeMap("label", nodeLabel).			 // read the 'label' node map into nodeLabel
				nodeMap("genome", inWhichGenome).			 // read the 'genome' node map into inWhichGenome
				attribute("number_of_genomes", mcfResultSt.numOfGenomes). // Not used
				attribute("number_of_mapping_reads", mcfResultSt.readInfoSt.numOfMappedReads).
				attribute("avg_read_length", mcfResultSt.readInfoSt.averageReadLength).
				attribute("max_score", mcfResultSt.maxScore).
				attribute("min_score", mcfResultSt.minScore).
				run();
				removeLowScoreArcs(g,inWhichGenome,arcWeight);
//				printRunningTime(tStart, "Finished reading data. Reading data took: ");
		} catch (Exception& error) { // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return ;
		}
}



set<int> getGenomesToErase(vector<GenomeInfo>& genomeInfoVector, int averageReadLength){
	set<int> genomesToErase;
	for(vector<GenomeInfo >::iterator genomesIterator=genomeInfoVector.begin();genomesIterator!=genomeInfoVector.end();genomesIterator++){
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId==-1)
			continue;
		if(gInfo.isLowCovered(averageReadLength, configSt))
			genomesToErase.insert(gInfo.genomeId);
	}
	return genomesToErase;
}
void eraseNodes(vector<ListDigraph::Node> nodesToErase, ListDigraph& g){
	for (vector<ListDigraph::Node>::iterator vit = nodesToErase.begin(); vit != nodesToErase.end(); ++vit)
	{
		g.erase(*vit);
	}
}
void eraseNodes(set<int>& genomesToErase, ListDigraph& g, ListDigraph::NodeMap<int>& 	inWhichGenome){
	int genomeId;
	vector<ListDigraph::Node> nodesToErase;
	for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
	{
		genomeId=inWhichGenome[node];
		if (genomeId > -1 && genomesToErase.count(genomeId))
		{
			nodesToErase.push_back(node);
		}
	}
	eraseNodes(nodesToErase,g);
	nodesToErase.clear();
	for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
	{
		if (inWhichGenome[node] == -1 && countOutArcs(g,node) == 0)
		{
			nodesToErase.push_back(node);
		}
	}
	eraseNodes(nodesToErase,g);
}
void trimDigraph(
	MCFResult& 								mcfResultSt,// will be updated
	ListDigraph& 							g, // will be updated
	ListDigraph::NodeMap<int>& 				inWhichGenome,
	ListDigraph::NodeMap<string>& 			nodeLabel
)
{
	mcfLogger.log("****************************trimDigraph Start*************************");
	mcfLogger.log("We have "+ getStringValue(mcfResultSt.getNumberOfGenomes())+ " Genomes.\t"+getStringValue(mcfResultSt.readInfoSt.numOfMappedReads)+" Reads.\t "+getStringValue(countNodes(g))+" Nodes.");
	set<int> genomesToErase=getGenomesToErase(mcfResultSt.genomeInfoVector, mcfResultSt.readInfoSt.averageReadLength);
	if(genomesToErase.size()>0){
		eraseNodes(genomesToErase, g, inWhichGenome);
		mcfResultSt.eraseGenome(genomesToErase);
		mcfResultSt.initAbundance( g,inWhichGenome,nodeLabel);
	}
	mcfLogger.log("After trimming, we are left with " + getStringValue(mcfResultSt.getNumberOfGenomes())+ " Genomes.\t"+getStringValue(mcfResultSt.readInfoSt.numOfMappedReads)+" Reads.\t "+ getStringValue(countNodes(g))+" Nodes.");
	mcfLogger.log("*******************************trimDigraph End*************************");
	mcfLogger.log("***********************************************************************");
}
void printRunningTime(clock_t& startTime, string msg){
	mcfLogger.log(msg +getStringValue((double)(clock()-startTime)/CLOCKS_PER_SEC/60 )+" min");
	startTime=clock();
}

void printNodes(ListDigraph& g, ListDigraph::NodeMap<int>& inWhichGenome, vector<string>& genomeNames, ListDigraph::NodeMap<string>& 	nodeLabel){
	map<string, set<int> > readGenomeMap;
	string genomeName, read;
	int genomeNum;
	for(ListDigraph::ArcIt arc(g);arc!=INVALID;++arc){
		read=nodeLabel[g.source(arc)];
		genomeNum=inWhichGenome[g.target(arc)];
		if(!readGenomeMap.count(read)){
			readGenomeMap.insert(pair<string,set<int> >(read,set<int>()));
		}
		readGenomeMap.find(read)->second.insert(genomeNum);
	}
	map<string, set<int> >::iterator mapItr=readGenomeMap.begin();
	set<int> genomeNumSet;
	set<int>::iterator firstItr;
	set<int>::iterator secondItr;
	map<string, int>finalMap;
	string key;
	int matrixSize=genomeNames.size();
	int gephiMatrix[matrixSize][matrixSize];
	for(int i=0;i <matrixSize; i++){
		for(int j=0;j <matrixSize; j++){
			gephiMatrix[i][j]=0;
		}
	}
	int row=0,column=0;
	while(mapItr!=readGenomeMap.end()){
		genomeNumSet=mapItr->second;
		firstItr=genomeNumSet.begin();
		while(firstItr!=genomeNumSet.end()){
			secondItr=genomeNumSet.begin();
			while(secondItr!=genomeNumSet.end()){
				if(*firstItr == *secondItr){
					++secondItr;
					continue;
				}
				row=*firstItr;
				column=*secondItr;
//				if(*firstItr > *secondItr){
//					row=*secondItr;
//					column=*firstItr;
//				}
				gephiMatrix[row][column]=gephiMatrix[row][column]+1;
//				gephiMatrix[column][row]=gephiMatrix[column][row]+1;
				++secondItr;
			}
			++firstItr;
		}
		++mapItr;
	}
	ofstream file("../../metagenomics/blast/101_WS_Log.Clusters");
	for(int i=0;i <matrixSize; i++){
		file << ";" <<genomeNames[i];
	}
	for(int i=0;i <matrixSize; i++){
		file <<	genomeNames[i];
		for(int j=0;j <genomeNames.size(); j++){
			file << ";" << gephiMatrix[i][j];
		}
	}
	file.close();
	mcfLogger.log("Finished Counting! ");
//	while(mapItr!=readGenomeMap.end()){
//		genomeNumSet=mapItr->second;
//		firstItr=genomeNumSet.begin();
//		while(true){
//			if(firstItr==genomeNumSet.end()){
//				break;
//			}
//			secondItr=firstItr;
//			++secondItr;
//			while(secondItr!=genomeNumSet.end()){
//				if(finalMap.count(genomeNames[*firstItr]+";"+genomeNames[*secondItr])){
//					key=genomeNames[*firstItr]+";"+genomeNames[*secondItr];
//				}else if(finalMap.count(genomeNames[*secondItr]+";"+genomeNames[*firstItr])){
//					key=genomeNames[*secondItr]+";"+genomeNames[*firstItr];
//				}else{
//					if(genomeNames[*firstItr].compare(genomeNames[*secondItr])<0)
//						key=genomeNames[*firstItr]+";"+genomeNames[*secondItr];
//					else
//						key=genomeNames[*secondItr]+";"+genomeNames[*firstItr];
//					finalMap.insert(pair<string, int>(key,0));
//				}
//				finalMap.find(key)->second=finalMap.find(key)->second+1;
//				++secondItr;
//			}
//			++firstItr;
//		}
//		++mapItr;
//	}
//	map<string, int>::iterator finalMapItr=finalMap.begin();
//
//	while(finalMapItr!=finalMap.end()){
//		file << finalMapItr->first<<endl;
////		file << finalMapItr->first<< ";"<<finalMapItr->second <<endl;
//		++finalMapItr;
//	}


}
void saveChunksResults(MCFResult& mcfResultSt, string fileName){
	mcfLogger.log("***********************saveGenomeAbundance Start*************************");
	ofstream abundanceWriter((char*)(fileName+".dist.csv").c_str());
		mcfLogger.log("Saving abundance results to file: "+fileName);
		vector<int> ::iterator chunkAbundanceItr;
		int maxNumOfChunks=0;
		for(vector<GenomeInfo >::iterator genomesIterator=mcfResultSt.genomeInfoVector.begin();genomesIterator!=mcfResultSt.genomeInfoVector.end();genomesIterator++){
			GenomeInfo& gInfo=*genomesIterator;
			if(gInfo.genomeId!=-1){
				maxNumOfChunks=max(maxNumOfChunks, gInfo.numOfChunks);
			}
		}
		vector<string> abundance;
		for(int i=0;i<maxNumOfChunks;i++){
			abundance.push_back(""+getStringValue(i));
		}
		string header="";
		for(vector<GenomeInfo >::iterator genomesIterator=mcfResultSt.genomeInfoVector.begin();genomesIterator!=mcfResultSt.genomeInfoVector.end();genomesIterator++){
				GenomeInfo& gInfo=*genomesIterator;
				int chunkNum=0;
				if(gInfo.genomeId!=-1){
					header=header+"\t"+gInfo.genomeName;
					for(vector<int>::iterator chunkAbundanceItr=gInfo.optimumAbundanceSt.chunksAbundanceVector.begin();
							chunkAbundanceItr!=gInfo.optimumAbundanceSt.chunksAbundanceVector.end();
							chunkAbundanceItr++){
						abundance[chunkNum]= abundance[chunkNum]+"\t"+getStringValue(*chunkAbundanceItr);
						chunkNum++;
					}
					for(int i=gInfo.numOfChunks;i<maxNumOfChunks;i++){
						abundance[i]=abundance[i]+"\tx";
					}
				}
		}
		abundanceWriter<<header<<endl;
		for(int i=0;i<maxNumOfChunks;i++){
			abundanceWriter<<abundance[i]<<endl;
		}
		abundanceWriter.close();
		mcfLogger.log("**************************saveGenomeAbundance End*********************");
		mcfLogger.log("**********************************************************************");
}
void saveChunksResults(MCFResult& mcfResultSt){
	string fileName=mcfResultSt.lgfFile;
	saveChunksResults(mcfResultSt,fileName);

}
void saveAbundanceResults(MCFResult& mcfResultSt, string fileName){
	mcfLogger.log("*****************************************************************************************");
	mcfLogger.log("*****************************************************************************************");
	mcfLogger.log("**********************************saveAbundanceResults Start*****************************");
	ofstream abundanceFile((char*)(fileName+".abundance.csv").c_str());
	double totalAbundance=0;
	for(vector<GenomeInfo >::iterator genomesIterator=mcfResultSt.genomeInfoVector.begin();genomesIterator!=mcfResultSt.genomeInfoVector.end();genomesIterator++){
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId!=-1){
			gInfo.optimumAbundanceSt.absoluteAbundance=gInfo.optimumAbundanceSt.genomeAbundance*mcfResultSt.readInfoSt.averageReadLength/(double)gInfo.genomeLength;
			totalAbundance+=gInfo.optimumAbundanceSt.absoluteAbundance;
		}
	}
	abundanceFile<<"Genome_Id\tGenome_Name\tGenome_Length\tNum_Of_Chunks\tNum_Of_Mapped_Reads\tAbsolute_Abundance\tRelative_Abundance"<<endl;
	for(vector<GenomeInfo >::iterator genomesIterator=mcfResultSt.genomeInfoVector.begin();genomesIterator!=mcfResultSt.genomeInfoVector.end();genomesIterator++){
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId!=-1){
			gInfo.optimumAbundanceSt.relativeAbundance=gInfo.optimumAbundanceSt.absoluteAbundance/totalAbundance;
			abundanceFile<<gInfo.genomeId<<"\t"<<gInfo.genomeName<<"\t"<<gInfo.genomeLength<<"\t"<<gInfo.numOfChunks<<"\t"<<gInfo.optimumAbundanceSt.genomeAbundance<<"\t"<<gInfo.optimumAbundanceSt.absoluteAbundance<<"\t"<<gInfo.optimumAbundanceSt.relativeAbundance<<endl;
		}
	}
	abundanceFile.close();
	mcfLogger.log("Abundance file: "+fileName+".abundance.csv is saved.");
	mcfLogger.log("******************************saveAbundanceResults End **********************************");
}
void saveAbundanceResults(MCFResult& mcfResultSt){
	string fileName= mcfResultSt.lgfFile;
	saveAbundanceResults(mcfResultSt, fileName);
}
void printRelativeAbundance(MCFResult& mcfresult){
	double totalAbundance=0;
	mcfLogger.log("*******************************printRelativeAbundance Start********************************");
	mcfLogger.log("********************************Final Abundance Results***********************************");

	for(vector<GenomeInfo >::iterator genomesIterator=mcfresult.genomeInfoVector.begin();genomesIterator!=mcfresult.genomeInfoVector.end();genomesIterator++){
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId!=-1){
			gInfo.optimumAbundanceSt.relativeAbundance=gInfo.optimumAbundanceSt.genomeAbundance*mcfresult.readInfoSt.averageReadLength/(double)gInfo.genomeLength;
			totalAbundance+=gInfo.optimumAbundanceSt.relativeAbundance;
		}
	}
	mcfLogger.log("Genome_Id\tGenome_Name\tOptimum_Abundance\tRelative_Abundance");
	for(vector<GenomeInfo >::iterator genomesIterator=mcfresult.genomeInfoVector.begin();genomesIterator!=mcfresult.genomeInfoVector.end();genomesIterator++){
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId!=-1){
			gInfo.optimumAbundanceSt.relativeAbundance/=totalAbundance;
			mcfLogger.log(getStringValue(gInfo.genomeId)+"\t"+gInfo.genomeName+"\t"+getStringValue(gInfo.optimumAbundanceSt.genomeAbundance)+"\t"+getStringValue(gInfo.optimumAbundanceSt.relativeAbundance));
		}
	}
	mcfLogger.log("********************************************************************************************");
}



void printGenomeInfo(vector<GenomeInfo>& genomeInfoVector){
	vector<GenomeInfo >::iterator genomesIterator=genomeInfoVector.begin();
	mcfLogger.log("****************************printGenomeInfo Start******************");
	mcfLogger.log("*******************************Genomes Info************************");
	mcfLogger.log("Genome_Id\tGenome_Name\tGenome_Length\tNum_Of_Chunks\tOptimum_Abundance\tExpectedAbundanceInChunk\tMin_Abundance\tNum_Of_Shared_Reads");
	for(vector<GenomeInfo >::iterator genomesIterator=genomeInfoVector.begin();genomesIterator!=genomeInfoVector.end();genomesIterator++){
		GenomeInfo& gInfo=*genomesIterator;
		if(gInfo.genomeId!=-1){
			mcfLogger.log(getStringValue(gInfo.genomeId)+"\t"+gInfo.genomeName+"\t"+getStringValue(gInfo.genomeLength)+"\t"+
			getStringValue(gInfo.numOfChunks)+"\t"+getStringValue(gInfo.optimumAbundanceSt.genomeAbundance)+"\t"+
			getStringValue(gInfo.optimumAbundanceSt.expectedChunkAbudance)+"\t"+
			getStringValue(gInfo.optimumAbundanceSt.genomeMinAbundance)+"\t"+
			getStringValue(gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads));
//			for(int i=0;i<gInfo.chunksOptimumAbundanceVector.size();i++){
//				cout<< gInfo.chunksOptimumAbundanceVector[i]<<"  ";
//			}
		}
//			cout<<gInfo.genomeId<<"\t"<<gInfo.genomeName<<"\t"<<gInfo.genomeLength<<"\t"<<gInfo.numOfChunks<<"\t"<<gInfo.optimumAbundance<<"\t"<<gInfo.minAbundance<<"\t"<<gInfo.maxAbundance<<endl;
	}
	mcfLogger.log("****************************printGenomeInfo End******************");
	mcfLogger.log("*****************************************************************");

}


/*
 * Removing acrs from ListDigraph.
 * @param arcsToErase: a vector of acrs to erase.
 * @param g: the graph object
 *
 */
void eraseArcs(vector<ListDigraph::Arc>& arcsToErase, ListDigraph& g)
{
	for (vector<ListDigraph::Arc>::iterator arcItr = arcsToErase.begin(); arcItr != arcsToErase.end(); ++arcItr)
	{
		g.erase(*arcItr);
	}
}
void eraseArcs(
		ListDigraph& 							g,
		ListDigraph::NodeMap<int>& 				inWhichGenome,
		ListDigraph::NodeMap<string>& 			nodeLabel,
		vector<ListDigraph::Arc>& 				arcsToErase)
{
	eraseArcs(arcsToErase,g);
	mcfLogger.log("Removing arcs finished. ");
}
void printGraphInfo(ListDigraph& g){
	mcfLogger.log("****************************printGraphInfo Start*******************");
	mcfLogger.log("***************Graph Information**************");
	mcfLogger.log("Total Number Of Node: "+getStringValue(countNodes(g)));
	mcfLogger.log("Total Number Of Arcs: "+getStringValue(countArcs(g)));
	mcfLogger.log("****************************printGraphInfo End*********************");
	mcfLogger.log("*******************************************************************");
}


int calcNumOfReadsWithLowPenalty(GenomeInfo& gInfo, string chunkName, bool removeSingleHitReads){
	int lowPenaltyReads=configSt.getIntProperty(NUM_OF_READS_WITH_LOWER_COST);
	if(removeSingleHitReads){
		int chunkNum=gInfo.getChunkNum (chunkName);
		if(gInfo.optimumAbundanceSt.chunksMinAbundanceVector[chunkNum] > gInfo.optimumAbundanceSt.expectedChunkAbudance){
			//lowPenaltyReads=lowPenaltyReads - (difference between the number of single hit reads in a chunk & the expected abundance in that chunk)
			// if the difference less than zero return zero
			lowPenaltyReads=lowPenaltyReads-(gInfo.optimumAbundanceSt.chunksMinAbundanceVector[chunkNum] - gInfo.optimumAbundanceSt.expectedChunkAbudance);
			return (lowPenaltyReads>0? lowPenaltyReads:0);
		}
		else
			return lowPenaltyReads;
	}else{
		return lowPenaltyReads;
	}
}
bool isExceedingMaxNumOfArcs(int numOfSharedArcs){
	if(numOfSharedArcs > configSt.getIntProperty(MAX_NUMBER_OF_ARCS))
		return true;
	return false;
}
void eraseArcs(
	MCFResult& 						mcfResultSt,
	ListDigraph& 					g,
	ListDigraph::NodeMap<int>& 		inWhichGenome,
	ListDigraph::NodeMap<string>& 	nodeLabel,
	ListDigraph::ArcMap<int>& 		arcWeight,
	int 							genomeId)
{
	int nodeGenomeId=0;
	int score=0;
	vector<ListDigraph::Arc> allArcsToErase;
	vector<ListDigraph::Arc> arcsToErase;
	for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
	{
		if (inWhichGenome[node] == -1)
		{
			if(countOutArcs(g,node) > 1)
			{
				arcsToErase.clear();
				score=10000000;
				for (ListDigraph::OutArcIt arc(g,node); arc != INVALID; ++arc)
				{
					nodeGenomeId=inWhichGenome[g.target(arc)];
					if(genomeId!=-1 && nodeGenomeId!=genomeId)
						continue;
					if(score==arcWeight[arc])
						continue;
					else if(score>arcWeight[arc]){
						score=arcWeight[arc];
						arcsToErase.insert(arcsToErase.begin(),arc);
					}
					else{
						arcsToErase.push_back(arc);
					}
				}
				if(arcsToErase.size()>0){
					arcsToErase.erase(arcsToErase.begin());
				}
				if(arcsToErase.size()>0){
					allArcsToErase.insert(allArcsToErase.end(),arcsToErase.begin(),arcsToErase.end());
					if(genomeId!=-1){
						GenomeInfo& gInfo=mcfResultSt.genomeInfoVector[genomeId];
						gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads -= arcsToErase.size()+1;
//						mcfLogger.log("Total Number of genome "+gInfo.genomeName+" shared reads= "+getStringValue(gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads));
						if(! isExceedingMaxNumOfArcs(gInfo.optimumAbundanceSt.numOfHitsForChunksSharedReads))
							break;
					}else{
						mcfResultSt.readInfoSt.numOfHitsForSharedReads -= arcsToErase.size()+1;
//						mcfLogger.log("Total Number of all shared reads= "+getStringValue(mcfResultSt.readInfoSt.numOfHitsForSharedReads));
						if(! isExceedingMaxNumOfArcs(mcfResultSt.readInfoSt.numOfHitsForSharedReads))
							break;
					}
				}
			}
		}
	}
	eraseArcs(allArcsToErase,g);
}
void removeOutOfBoundaryHits(
	MCFResult& mcfResultSt,
	ListDigraph& 					g,
	ListDigraph::NodeMap<int>& 		inWhichGenome,
	ListDigraph::NodeMap<string>& 	nodeLabel
){
	int chunkNum;
	vector<ListDigraph::Arc> outOfBoundaryArcs;
	map<string, int> erasedArcsMap;
	for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
	{
		if (inWhichGenome[node] == -1)
		{
			for (ListDigraph::OutArcIt arc(g,node); arc != INVALID; ++arc)
			{
				GenomeInfo& gInfo=mcfResultSt.genomeInfoVector[inWhichGenome[g.target(arc)]];
				chunkNum=GenomeInfo::getChunkNum(nodeLabel[g.target(arc)]);
				if(chunkNum > gInfo.optimumAbundanceSt.chunksAbundanceVector.size()-1){
					outOfBoundaryArcs.push_back(arc);
					if(erasedArcsMap.find(gInfo.genomeName)==erasedArcsMap.end()){
						erasedArcsMap[gInfo.genomeName]=0;
					}
					erasedArcsMap[gInfo.genomeName]=erasedArcsMap[gInfo.genomeName]+1;
				}
			}
		}
	}
	eraseArcs(outOfBoundaryArcs, g);
	int numberOfErasedNodes=0;
	vector<ListDigraph::Node> nodesToEraseVector;
	for (ListDigraph::NodeIt node(g); node != INVALID; ++node)
	{
		if (inWhichGenome[node] == -1)
		{
			if(countOutArcs(g,node)==0){
				nodesToEraseVector.push_back(node);
				numberOfErasedNodes++;
			}
		}
	}
	mcfLogger.log("*******************Out Of Boundary Hits********************");
	mcfLogger.log("Total Number of removed Reads= "+getStringValue(numberOfErasedNodes));
	eraseNodes(nodesToEraseVector, g);
	for(map<string, int>::iterator itr=erasedArcsMap.begin();itr!=erasedArcsMap.end();itr++){
		mcfLogger.log(itr->first+"\t"+getStringValue(itr->second));
	}
	mcfLogger.log("*******************Out Of Boundary Hits End********************");
}
