# MetaFlow
MetaFlow is a program for community profiling of a metagenomic sample.

** Compiling the source code:
[This part will be updated when we finish and try compiling it]

** Running the tool:
./MCF_Abundance_Estimation {LGF} {Genome File}

** INPUT:
The tool takes two input files:
1- LGF(Lemon Graph File):-
The aligner output file should be transformed into LGF format to be processed by our tool. We provide a python script file which converts the tabular Blast output (with format=6) 
into our mapping file. If the user wishes to use another aligner, he should write his script for transforming this aligner’s output into our mapping file format.
The file formate is explained below.

2- Genome File:-
Genome file contains as in NCBI bacterial. We provide a file for NCBI Bacterial database taken on 10/06/2015. The file location is Tool_Home\NCBI\NCBI_Ref_Genome.txt.
If you are using different or updated database, you need to update or change the genome file to incorporate all the reference genomes. 
The file format is explained below.    


** OUTPUT:
The tool generates different output files. All the files are in CSV format to make any further analysis easy. All the output files will be generated in the same
folder where the LGF file is located, expect for the log file which will be generated in the Log folder.

|-------------------------------|-----------------------------------------------------------------------------------------------|
|	File 						|	Description																					|
|-------------------------------|-----------------------------------------------------------------------------------------------|
| abundance.csv 				| The main output file. It contains the final estimation of the species richness and abundance.	|
|-------------------------------|-----------------------------------------------------------------------------------------------|
| dist.csv 						| It contains the final distribution of the reads over all chunks in all genomes.				|
|-------------------------------|-----------------------------------------------------------------------------------------------|
| stage0.abundance.csv			| Intermediary internal files that contain the estimation of the species richness and abundance	|
| stage1.abundance.csv			| before starting and after the stages 1,2, and 3, respectively.								|
| stage2.abundance.csv			|																								|
| stage3.abundance.csv			|																								|
|-------------------------------|-----------------------------------------------------------------------------------------------|
| stage0.dist.csv				| Intermediary internal files that contains the distribution of the reads over the chunks before|
| stage1.dist.csv				| starting and after the stages 1,2, and 3, respectively.										|
| stage2.dist.csv				|																								|
| stage3.dist.csv				|																								|
|-------------------------------|-----------------------------------------------------------------------------------------------|
| .log 							| The running log file.																			|
|-------------------------------|-----------------------------------------------------------------------------------------------|


**** Genome File (if you are writing your own script to create the file):-
The genome file contains a list of the bacterial genomes and their length. Each line is GenomeName\tGenomeLength, with the following rules:
1- GenomeName is GeneraName_SpeciesName,
2- GenomeLength is the length of the genome. If one species has different strains with different lengths, select the shortest one.
3- Ignore Plasmid.

A simple example for the file:

Corynebacterium_diphtheriae	2395441
Methylomonas_methanica	5051681
Paenibacillus_sp._JDR-2	7184930
Psychroflexus_torquis	4321832
Shewanella_woodyi	5935403
Erwinia_amylovora	3805573

**** LGF Format (if you are writing your own script to create the file):-
The file formate is:
1- The file starts with two lines marks the beginning of:
@nodes
label\tgenome

2- A list of [GenomeId_ChunkNumber\tGenomeId]. For example if a read r1 maps to Brucella_ovis. The genomeId of Brucella_ovis in the genome database file is 402.
The genome length in the file is 1164220. If the chunk size we use is 2000, then we will have 583 chunks (numbered from 0 to 582) for this genome. These 583 chunks
should be added as following:
402_0	402
402_1	402
.....	...
.....	...
402_582	402

3- A list of all reads in the form [ReadName\t-1]. For example if we have 10 reads named in the fasta file(r1, r2, ....r10), we add the following lines
r1	-1
r2	-1
..	-1
..	-1
r10	-1

4- Two lines mark the start of mapping from read to chunks.
@arcs 
\t\tlabel\tweight

5- Mappings from read to chunks. if there is a read "r1" that maps to chunk "3" in genome "1", add a line [r1\t1_3\tCounter\tCost], Where Counter is a regular counter
that starts from 0, and Cost is the transformation of the alignment score.

6-One line marks the 
@attributes

7-Summary. If the total number of reads is 100, with average read length 50, and they map to 10 genomes, and maximum cost is 150, minimum cost is 100, add the following lines:
number_of_genomes	10
number_of_mapping_reads	100
avg_read_length	50
max_score	150
min_score	100


An simple example of how the LGF file looks like:-

@nodes
label	genome
0_0	0
0_1	0			
0_2	0
0_3	0
0_4	0
1_0	1
1_1	1
1_2	1
1_3	1
1_4	1
2_0	2
2_1	2
2_2	2
r1	-1
r2	-1
r3	-1
r4	-1
@arcs
		label	weight
r1	1_0	0	100
r1	2_1	1	110
r2	1_1	2	118
r2	1_0	3	111
r2	1_0	4	100
r3	0_2	5	105
r3	0_3	6	100
r3	0_4	7	100
r3	2_2	8	113
r4	1_3	9	100
@attributes
number_of_genomes	10
number_of_mapping_reads	100
avg_read_length	50
max_score	150
min_score	100


*Note: \t is replaced with actual tab in the file.