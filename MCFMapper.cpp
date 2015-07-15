
#include "MCFHeaders.h"
#include "MCFConfig.h"
#include "MCFFlowSolver.h"
#include "MCFGenetic.h"
#include "MCFLogger.h"
#include "OptionParser.h"

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

	// command line argument parser
	string usage = "\n  %prog OPTIONS";
	const string version = "%prog 0.9.1\nCopyright (C) 2013-2015\n"
		"License GPLv2+: GNU GPL version 2 "
		"<http://gnu.org/licenses/gpl.html>.\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.";
	const string desc = "MetaFlow is a tool for community profiling of a metagenomic sample."
		" It reports the known species present in a metagenomics sample and their abundances "
		"(relative to the known reference genomes).";
	const string epilog = "";
	
	optparse::OptionParser parser = optparse::OptionParser()
    	.usage(usage)
    	.version(version)
    	.description(desc)
    	.epilog(epilog);

	parser.add_option("-m", "--mappings") .type("string") .dest("m") .set_default("") .help("The read mapping file, in LGF format. This can be produced from BLAST alignments (tabular format with format=6) with the enclosed Python script BLAST_TO_LGF.py. See the manual for more details.");
	parser.add_option("-g", "--genome") .type("string") .dest("g") .set_default("") .help("Genome file. You can use: 'NCBI/NCBI_Ref_Genome.txt'. See the manual for more details.");
	parser.add_option("-c", "--config") .type("string") .dest("c").set_default("") .help("The config file. You can use/modifiy: 'metaflow.config'. See the manual for more details.");

	optparse::Values& options = parser.parse_args(argc, argv);

	string mappings_file = (string) options.get("m");
	string genome_file = (string) options.get("g");
	string config_file = (string) options.get("c");
	
	if (mappings_file == "")
	{
		cerr << "Parameter -m|--mappings should not be empty. Run ./metaflow -h|--help for details." << endl;
		return EXIT_FAILURE;
	}
	if (genome_file == "")
	{
		cerr << "Parameter -g|--genome should not be empty. Run ./metaflow -h|--help for details." << endl;
		return EXIT_FAILURE;
	}
	if (config_file == "")
	{
		cerr << "Parameter -c|--config should not be empty. Run ./metaflow -h|--help for details." << endl;
		return EXIT_FAILURE;
	}

	configSt = MCFConfig(config_file);

	mcfLogger.log("*************************main Start*******************");
	clock_t startTime=clock();
	int result=0;
	result=runMCFMapper(mappings_file.c_str(),genome_file.c_str());

	mcfLogger.log("**********************main End************************");
	mcfLogger.log("*****************************************************");
	printRunningTime(startTime, "Finding Optimal Abundance Total Running Time: ");
	mcfLogger.log("Finished!!");
	mcfLogger.close();
	return result;

}
