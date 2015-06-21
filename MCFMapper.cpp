
#include "MCFHeaders.h"
#include "MCFConfig.h"
#include "MCFFlowSolver.h"
#include "MCFGenetic.h"
#include "MCFLogger.h"

using namespace lemon;
//using namespace boost;
using namespace std;
//using namespace algorithm;

namespace lemon{
	Random rnd;
}
MCFLogger mcfLogger;
MCFConfig configSt;
clock_t tStart;
/*
 * The entry function for running the tool.
 * @param mcfFinalResultSt: the Struct aggregating all results.
 * @return return 0 in case of success. Should return -1 in case of any exceptions (not implemented)
 */
int runMCFMapper(MCFResult& mcfFinalResultSt)
{
		mcfLogger.log("*********************************runMCFMapper Start***********************************");

		if(configSt.getIntProperty(CALC_NAIVE_BLAST)==1)
			saveNaiveBlastAbundanceResults(mcfFinalResultSt.lgfFile, mcfFinalResultSt.ncbiRefGenomesFile);

		//******* Lemon library objects******//
		ListDigraph g;   //
		ListDigraph::ArcMap<int> arcWeight(g); // A map holds the acrs weights (Blast score)
		ListDigraph::NodeMap<string> nodeLabel(g); // A map holds the nodes label (represents the genomes chunks) in the form (GenomeId_ChunkNumber)
		ListDigraph::NodeMap<int> inWhichGenome(g); // a map holds the genomes ids in the form (0,1,2,...)


		// mcfResultSt: a Struct for handling the reads, and the species and their abundance.
		MCFResult mcfResultSt(mcfFinalResultSt.lgfFile, mcfFinalResultSt.ncbiRefGenomesFile);
		mcfResultSt.chunkSize=configSt.getIntProperty(CHUNK_SIZE);

		initDigraph(mcfResultSt,g,inWhichGenome,nodeLabel, arcWeight);
		printRunningTime(tStart, "Finished graph initialization. Graph initialization took: ");
		printGenomeInfo(mcfResultSt.genomeInfoVector);
		printGraphInfo(g);
		saveAbundanceResults(mcfResultSt,mcfResultSt.lgfFile+".Step0");
		saveChunksResults(mcfResultSt,mcfResultSt.lgfFile+".Step0");

		trimDigraph(mcfResultSt, g, inWhichGenome, nodeLabel);
		printRunningTime(tStart, "Finished trimming graph. Trimming graph took: ");
		printGenomeInfo(mcfResultSt.genomeInfoVector);
		printGraphInfo(g);
		saveAbundanceResults(mcfResultSt,mcfResultSt.lgfFile+".Step1");
		saveChunksResults(mcfResultSt,mcfResultSt.lgfFile+".Step1");

		FlowNetworkConfig flowNetworkConfigSt(-1,false, true); // Don't use Unknown node, Don't remove single hit
		findMaxAbundance(mcfResultSt,g,inWhichGenome,nodeLabel,arcWeight,flowNetworkConfigSt);
		printRunningTime(tStart, "Finished finding Max abundance. Finding Max Abundance took: ");
		printGenomeInfo(mcfResultSt.genomeInfoVector);
		printGraphInfo(g);
		saveAbundanceResults(mcfResultSt,mcfResultSt.lgfFile+".Step2");
		saveChunksResults(mcfResultSt,mcfResultSt.lgfFile+".Step2");
		mcfResultSt.copyLocalAbundance();

		greedy_findAbundances(mcfResultSt, g, inWhichGenome, nodeLabel, arcWeight);
		printRunningTime(tStart, "Finished finding optimum abundance. Finding Optimum Abundance took: ");
		printGenomeInfo(mcfResultSt.genomeInfoVector);
		printGraphInfo(g);
		saveAbundanceResults(mcfResultSt,mcfResultSt.lgfFile+".Step3");
		saveChunksResults(mcfResultSt,mcfResultSt.lgfFile+".Step3");

		flowNetworkConfigSt.useUnknownNode=true;
		flowNetworkConfigSt.removeSingleHitReads=false; // Use Unknown node, Don't remove single hit
		findMaxAbundance(mcfResultSt,g,inWhichGenome,nodeLabel,arcWeight,flowNetworkConfigSt);
		printRunningTime(tStart, "Finished finding Max abundance. Finding Max Abundance took: ");
		printGenomeInfo(mcfResultSt.genomeInfoVector);
		printGraphInfo(g);
		saveAbundanceResults(mcfResultSt,mcfResultSt.lgfFile+".Step4");
		saveChunksResults(mcfResultSt,mcfResultSt.lgfFile+".Step4");
		mcfResultSt.copyLocalAbundance();
		trimDigraph(mcfResultSt, g, inWhichGenome, nodeLabel);
		printRunningTime(tStart, "Finished trimming graph. Trimming graph took: ");
		printGenomeInfo(mcfResultSt.genomeInfoVector);
		printGraphInfo(g);

		mcfFinalResultSt.copyResult(mcfResultSt);

		mcfLogger.log("***********************************runMCFMapper End***********************************");
		mcfLogger.log("**************************************************************************************");
		return 0;
}
int runMCFMapper(string clustersFilesNames){
	double aveReadLength=0;
	map<string,string> lgfFilesMap;
	ifstream clustersFile(clustersFilesNames.c_str());
	if (!clustersFile.is_open()) { // check for successful opening
		mcfLogger.log("Error. Can't open clusters file! " );
		return -1;
	}

	string line;
	vector<string> splittedLine;
	string filePath="";
	split(splittedLine, clustersFilesNames, "/");
	for(int i=0;i<splittedLine.size()-1;i++){
		filePath+=splittedLine[i]+"/";
	}
	while(getline(clustersFile,line)){
		split(splittedLine, line, "\t");
		lgfFilesMap.insert(pair<string, string>(filePath+splittedLine[0],filePath+splittedLine[1]));
	}
	map<string,string>::iterator lgfFilesItr=lgfFilesMap.begin();
	MCFResult mcfFinalResult("","");
	while(lgfFilesItr!=lgfFilesMap.end()){
		string lgfFile=lgfFilesItr->first;
		mcfFinalResult.lgfFile=lgfFile;
		string genomeNamesFile=lgfFilesItr->second.c_str();
		mcfFinalResult.ncbiRefGenomesFile=genomeNamesFile;
		mcfLogger.log(lgfFile+"\t"+lgfFilesItr->second);
		if(runMCFMapper(mcfFinalResult)!=0){
			return -1;
		}
		lgfFilesItr++;
	}
	saveChunksResults(mcfFinalResult);
	saveAbundanceResults(mcfFinalResult);
	printRelativeAbundance(mcfFinalResult);
	return 0;
}
int runMCFMapper(string lgfFile, string ncbiRefFile){
	mcfLogger.log("*************************runMCFMapper Start*******************");
	MCFResult mcfFinalResult(lgfFile, ncbiRefFile);
	mcfLogger.log("LGF File: "+lgfFile);
	mcfLogger.log("NCBI File: "+ncbiRefFile);
	// mcfFinalResult is the object aggregating the results from the clustered samples. In the new implementation, there is only one sample.
	int result=runMCFMapper(mcfFinalResult);
	saveChunksResults(mcfFinalResult);
	saveAbundanceResults(mcfFinalResult);
	printRelativeAbundance(mcfFinalResult);
	mcfLogger.log("**********************runMCFMapper End************************");
	mcfLogger.log("**************************************************************");
	return result;
}

int main(int argc, char **argv) {
	mcfLogger.log("*************************main Start*******************");
	/*
	 * Always pass two parameters. One parameter is an old implementation, in which reads were clustered into different samples to speed up the tool.
	 * argc[1]: The lgf mapping file.
	 * argc[2]: NCBI reference genome file. Currently we use the file (NCBI_Ref_Genome.txt) inside the folder NCBI.
	 */
	if (argc<2) {
		mcfLogger.log("Parameters: [MappingFile] [ReferenceGenome]" );
		return 0;
	}
	clock_t startTime=clock();
	int result=0;
	if(argc==3){
		result=runMCFMapper(argv[1],argv[2]);
	}
	else{
		result= runMCFMapper(argv[1]);
	}
	mcfLogger.log("**********************main End************************");
	mcfLogger.log("*****************************************************");
	printRunningTime(startTime, "Finding Optimal Abundance Total Running Time: ");
	mcfLogger.log("Finished!!");
	mcfLogger.close();
	return result;
}



/*
g++ MCFmapper.cpp MCFutils.cpp MCFgenetic.cpp MCFflowsolver.cpp -I lemon/include -L lemon/lib -lemon -o MCFmapper -O3 -pedantic -Wall -Wno-long-long
./MCFmapper mapping.lgf genome-names.txt
g++ MCFmapper.cpp MCFutils.cpp MCFgenetic.cpp MCFflowsolver.cpp -I /cs/fs/home/tomescu/lemon/include -L /cs/fs/home/tomescu/lemon/lib -lemon -o MCFmapper -O3 -pedantic -Wall -Wno-long-long
*/
