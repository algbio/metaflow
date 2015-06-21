
#ifndef CONGIF_H_
#define CONGIF_H_

#include "MCFHeaders.h"
#include "MCFLogger.h"
#include "MCFUtils_Temp.h"


#define CHUNK_SIZE "CHUNK_SIZE"
#define ALPHA "ALPHA"
#define TRIMMING_PERCENTAGE "TRIMMING_PERCENTAGE"
#define SPLIT_ARCS	"SPLIT_ARCS"
#define NUM_OF_READS_WITH_LOWER_COST "NUM_OF_READS_WITH_LOWER_COST"
#define REQUIRED_MIN_ABUNDANCE "REQUIRED_MIN_ABUNDANCE"
#define REQUIRED_AVERAGE_CHUNKS_COVERAGE "REQUIRED_AVERAGE_CHUNKS_COVERAGE"
#define REQUIRED_MAX_PER_OF_EMPTY_CHUNKS "REQUIRED_MAX_PER_OF_EMPTY_CHUNKS"
#define MAX_NUMBER_OF_MUTATION_LOOPS "MAX_NUMBER_OF_MUTATION_LOOPS"
#define MAX_COST_DIFFERENCE "MAX_COST_DIFFERENCE"
#define MAX_READS_DIFFERENCE "MAX_READS_DIFFERENCE"
#define MAX_SCORE_DIFFERENCE "MAX_SCORE_DIFFERENCE"
#define MAX_RUNNING_TIME "MAX_RUNNING_TIME"
#define MAX_NUMBER_OF_ARCS "MAX_NUMBER_OF_ARCS"

#define CALC_NAIVE_BLAST "CALC_NAIVE_BLAST"
#define NAIVE_BLAST_MIN_ABUNDANCE "NAIVE_BLAST_MIN_ABUNDANCE"
#define int64_t_MAX std::numeric_limits<int64_t>::max()
//#define UNKWON_COST "UNKWON_COST"
//#define FITLER_READS "FITLER_READS" // For use unknown node
//#define EXTRA_READS_REGULAR_COST "EXTRA_READS_REGULAR_COST"
//#define EXTRA_READS_LENIENT_COST "EXTRA_READS_LENIENT_COST"
//#define LOW_READS_REGULAR_COST "LOW_READS_REGULAR_COST"
//#define LOW_READS_LENIENT_COST "LOW_READS_LENIENT_COST"
//#define GENERATE_CHUNKS_COVERAGE_FILE	"GENERATE_CHUNKS_COVERAGE_FILE"
//#define GENERATE_INITIAL_MAPPING_RESULTS_FILE "GENERATE_INITIAL_MAPPING_RESULTS_FILE"
//#define NCBI_REFERENCE_FILE_PATH "NCBI_REFERENCE_FILE_PATH"
//#define NCBI_TAXONOMY_FILE_PATH "NCBI_TAXONOMY_FILE_PATH"
//#define OUTPUT_PATH	"OUTPUT_PATH"

extern MCFLogger mcfLogger;
struct MCFConfig{
	map<string, string> configParams;
	MCFConfig(){
		ifstream confFile("MCF.config");
		string line;
		vector<string>splittedLine;
		mcfLogger.log("Configuration Parameters: ");
		while (getline(confFile, line)){
			if(line[0]!='#'){
				split(splittedLine, line, "\t");
				configParams[splittedLine[0]]=splittedLine[1];
				mcfLogger.log(splittedLine[0]+"\t"+splittedLine[1]);
			}
		}
	}
	string getStrProperty(string propName){
		string propertyValue=configParams[propName];
		//mcfLogger.log("Property ["+propName+"]="+propertyValue);
		return propertyValue;
	}
	long int getIntProperty(string propName){
		return atoi(getStrProperty(propName).c_str());
	}
	double getDoubleProperty(string propName){
		return atof(getStrProperty(propName).c_str());
	}
	bool getboolProperty(string propName){
		return (atoi(getStrProperty(propName).c_str())!= 0);
	}
};

#endif /* CONGIF_1_H_*/

