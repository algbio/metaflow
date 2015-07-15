#ifndef MCFGENOMEINFO_H_
#define MCFGENOMEINFO_H_

#include "MCFHeaders.h"
#include "MCFConfig.h"
#include "MCFLogger.h"
#include "MCFUtils_Temp.h"
extern MCFLogger mcfLogger;
extern MCFConfig configSt;
/*
 * A struct contains the abundance information for each genome.
 */
struct AbundanceInfo{
	int genomeMinAbundance;					// The number of single hit reads assigned to the genome.
	int genomeAbundance;					// The total number of reads assigned to a genome.
	int expectedChunkAbudance;				// The mean coverage expected for all chunks to achieve uniform distribution.
	double absoluteAbundance;				// (number of reads*read length)/genome length
	double relativeAbundance;				// abunance relative to the species in the sample
	int numOfHitsForChunksSharedReads;		// Total number of arcs for reads mapped to more than one chunk in the same genome. This is used for heuristic trimming of graph by removing arcs if the number of arcs is higher than MAX_NUMBER_OF_ARCS.
	vector<int> chunksAbundanceVector;		// A vector contains the actual abundance of each chunk in the genome.
	vector<int> chunksMinAbundanceVector;	// A vector contains the number of single hits for each chunk.
	AbundanceInfo(){
		genomeAbundance=0;
		genomeMinAbundance=0;
		expectedChunkAbudance=0;
		relativeAbundance=0;
		absoluteAbundance=0;
		numOfHitsForChunksSharedReads=0;
	}
};

/*
 * A struct contains the information for a single genome.
 */
struct GenomeInfo{
	int genomeId;			// a number represents the line number for the genome information in the NCBI reference genome file.
							// This number if used for mapping between mapping information in the lemon graph object and the genomeInfo struct.
	bool hasHits; 			// This parameter used for printing statistics about genomes that were originally in the sample and were excluded by the MCF.
							//If the genome has genomeId !=-1 and hasHits =false, it means it was in the sample but identified as an outlier by the our tool.
	string genomeName;		// The actual name of the species.
	int genomeLength;		// The length of the genome.
	int numOfChunks ;		// The total number of chunks in the genome.
	AbundanceInfo optimumAbundanceSt; // Contains the abundance information for the genome. This information is updated during running the MCF till it reaches the final optimal abundance of the species. See the abundance information Struct.
	AbundanceInfo localOptAbundanceSt; //
	GenomeInfo(){
		genomeId=-1;
		genomeLength=0;
		numOfChunks=0;
		hasHits=false;
	}

	/*
	 * Copying the abundance information to localOptAbundanceSt
	 */
	void copyLocalAbundance(){
		if(hasHits!=-1){
			localOptAbundanceSt.genomeAbundance=optimumAbundanceSt.genomeAbundance;
			localOptAbundanceSt.expectedChunkAbudance=optimumAbundanceSt.expectedChunkAbudance;
			localOptAbundanceSt.chunksAbundanceVector.insert(localOptAbundanceSt.chunksAbundanceVector.end(),optimumAbundanceSt.chunksAbundanceVector.begin(),optimumAbundanceSt.chunksAbundanceVector.end());
		}
	}
	/*
	 * Extracting the chunk number from the chunk label . For example 10_3 returns 3.
	 */
	static int getChunkNum(string chunkName){
		vector<string> genomeChunk;
		split(genomeChunk, chunkName,"_");
		return atoi(genomeChunk[1].c_str());
	}
	/*
	 * @param chunkName: The chunk name in the lemon graph object. For example 10_3.
	 * @param removeSingleHitReads: based on this boolean value, the expected chunk abundance is recalculated. For example if the single hits are removed from the network flow,
	 * 			the expected abundance will be the original expected abundance - the number of single hit reads removed.
	 * @return : The expected chunk abundance
	 */
	int getChunkExpecetedAbundance(string chunkName, bool removeSingleHitReads){
		int chunkNum=getChunkNum (chunkName);
		if(removeSingleHitReads){
			// If the number of single hit reads >= the expected abundance, then the chunk shouldn't demand for more reads.
			if(optimumAbundanceSt.chunksMinAbundanceVector[chunkNum] >= optimumAbundanceSt.expectedChunkAbudance)
				return 0;
			else
				// The chunk demands for the difference between the original expected and the covered by the single hit reads.
				return optimumAbundanceSt.expectedChunkAbudance - optimumAbundanceSt.chunksMinAbundanceVector[chunkNum];
		}else{
			return optimumAbundanceSt.expectedChunkAbudance;
		}
	}

	/*
	 *  Initializing the genome abundance to zero
	 */
	void initGenome(){
		this->optimumAbundanceSt.genomeAbundance=0;
		this->optimumAbundanceSt.genomeMinAbundance=0;
		this->optimumAbundanceSt.relativeAbundance=0;
		this->optimumAbundanceSt.numOfHitsForChunksSharedReads=0;
		int numOfChunks=this->optimumAbundanceSt.chunksAbundanceVector.size();
		if(this->genomeId!=-1){
			for(int i=0; i<numOfChunks;i++){
				this->optimumAbundanceSt.chunksAbundanceVector[i]=0;
				this->optimumAbundanceSt.chunksMinAbundanceVector[i]=0;
			}
		}
	}

	/*
	 * calculates the expected abundance for the chunk based on the total number of reads assigned to the genome and their distribution between the chunks.
	 * For uniform expected abundance for chunks, the trimmed mean is used for calculating the abundance.
	 */
	void calcExpectedAbundanceInChunk(){
		//20% Trimmed Mean
		vector<int> chunkAbundanceCopy(optimumAbundanceSt.chunksAbundanceVector);
		sort(chunkAbundanceCopy.begin(),chunkAbundanceCopy.end());
		int trimmingPercent=chunkAbundanceCopy.size()*configSt.getDoubleProperty(TRIMMING_PERCENTAGE); // Remember to change the code and make it a configuration parameter.
		int totalSum=0;
		int i=0;
		int end=chunkAbundanceCopy.size()-trimmingPercent;
		for(i=trimmingPercent;i<end;i++){
			totalSum+=chunkAbundanceCopy[i];
		}
		optimumAbundanceSt.expectedChunkAbudance=round ((double)totalSum/(double)(i-trimmingPercent));
		mcfLogger.log("Expected Chunk Abundance for Genome "+genomeName+"="+getStringValue(optimumAbundanceSt.expectedChunkAbudance));
		//	Mean
		//	return round((double)gInfo.optimumAbundance/ (double)gInfo.numOfChunks);
	}

	/*
	 *
	 */
	double calcPercentageOfCoveredChunks(){
		return (double)optimumAbundanceSt.genomeAbundance/(double)numOfChunks;
	}
	double calcGenomeAbundance(double& averageReadLength){
		return (double)(optimumAbundanceSt.genomeAbundance*averageReadLength)/(double)genomeLength;
	}
	double calcPercentageOfLowCoveredChunks(MCFConfig& configSt){
		int numOfEmptyChunks=0;
		int minNumOfReadsPerChunk=configSt.getIntProperty(REQUIRED_AVERAGE_CHUNKS_COVERAGE);
		for(vector<int>::iterator itr=optimumAbundanceSt.chunksAbundanceVector.begin();itr!=optimumAbundanceSt.chunksAbundanceVector.end();itr++){
			if(*itr<minNumOfReadsPerChunk)
				numOfEmptyChunks++;
		}
		return (double)numOfEmptyChunks/(double)numOfChunks;
	}
	bool isLowCovered(double averageReadLength,MCFConfig& configSt){
		mcfLogger.log("********************checkCoverage Start****************");
		bool lowCovereged=false;
		double perOfCoveredChunks=calcPercentageOfCoveredChunks();
		double genomeAbundance=calcGenomeAbundance(averageReadLength);
		double requiredAvgChunksCoverage=configSt.getDoubleProperty(REQUIRED_AVERAGE_CHUNKS_COVERAGE);
		double requiredMinAbundance=configSt.getDoubleProperty(REQUIRED_MIN_ABUNDANCE);
		mcfLogger.log("Genome_name="+genomeName+"\tAbsolute Abundance="+getStringValue(genomeAbundance)+"\tChunks Average Coverage="+getStringValue(perOfCoveredChunks));
		if(perOfCoveredChunks< requiredAvgChunksCoverage){
			lowCovereged=true;
			mcfLogger.log("Delete genome: GenomeId="+getStringValue(genomeId)+"\tGenome_name="+genomeName+"\tChunks Average Coverage="+getStringValue(perOfCoveredChunks)+"\tRequired Chunks Average Coverage="+getStringValue(requiredAvgChunksCoverage));
		}
		else if(genomeAbundance< requiredMinAbundance){
			lowCovereged=true;
			mcfLogger.log("Delete genome: GenomeId="+getStringValue(genomeId)+"\tGenome_name="+genomeName+ "\tAbundance="+getStringValue(genomeAbundance)+"\tRequired abundance="+getStringValue(requiredMinAbundance));
		}
		else{
			double emptyChunks=calcPercentageOfLowCoveredChunks(configSt);
			mcfLogger.log("% of empty chunks="+getStringValue(emptyChunks));
			double allowedEmptyChunks=configSt.getDoubleProperty(REQUIRED_MAX_PER_OF_EMPTY_CHUNKS);
			if(emptyChunks >allowedEmptyChunks){
				lowCovereged=true;
				mcfLogger.log("Delete genome: GenomeId="+getStringValue(genomeId)+"\tGenome_name="+genomeName+ "\t% of empty chunks="+getStringValue(emptyChunks)+"\tAllowed Empty Chunks="+getStringValue(allowedEmptyChunks));
			}
		}
		mcfLogger.log("********************checkCoverage End****************");
		mcfLogger.log("*****************************************************");
		return lowCovereged;
	}
};


#endif /* MCFGENOMEINFO_H_ */
