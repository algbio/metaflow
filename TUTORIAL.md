# MetaFlow Tutorial

Here is a step-by-step tutorial on using *MetaFlow*. 

We have the file **Tutorial/mock_sample.fna** of 454 Pyrosequencing reads that we want to analyze. This file looks like this:

	>r1.1 
	GCGGAAACGGCTTTAAAACAAGCAAAAGAAAACAAATACGAAAATTACTTCGCGTTGACTAAGTC
	AGACAAAACCGACACGGGCAGAAGTTTGGCGCTAAAAGCAGATTTAAAACGTGCATTAGCACAAA
	ATGAACTCGAGCTTTACTATCAACCAAAAGTCGATCTGACAACGCTAGAAGTAGTGGGGGCTGAG
	TGCTTACTTCGCTGGAATCACCCTCTTGATGGCGTCTTATTTCCAGGTCC
	>r2.1 
	GCGACTGCCTCGTAGTCTACATGTTCTACGTTACGTGCATGTTCAAGCGTGGCTTTAAACTCATC
	ACTGCGCACTACCGCTTGTACTGATGCATCGTCATACCCGTCAATGGCAGTCACATCGATATACA
	AATAATTAAGCCAGCGACGTGAGCTTGGTCCATAAGGCGAGCACGCGTTGGGATTTGCGGGGTAA
	AGCGCATGAATAGGATTGAGTCCGATAAAGTCTGCGCCCAC
	...

*MetaFlow* takes as input an LGF file containing the alignments of these reads inside a database of reference bacterial genomes. In our experiments we aligned the reads with BLAST; we describe this workflow onwards. However, you can use any other aligner, provided that you convert its output into the right input format for *MetaFlow* (see the [complete manual](https://github.com/alexandrutomescu/metaflow/blob/master/MANUAL.md) for how to do this).

## Step 0

We install *MetaFlow* from github by running

	git clone https://github.com/alexandrutomescu/metaflow.git
	
This creates the directory **metaflow**. We go to it and compile:

	cd metaflow
	make

This creates the executable called **metaflow** in the same directory.


## Step 1 

First, we need to install BLAST, as described e.g. [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). Let's assume that we installed BLAST in the directory

	~/ncbi-blast-2.3.0+/

## Step 2

Then, we need to create a BLAST database of complete reference bacterial genome. We use the Python script *Create_Blast_DB.py* of **MetaFlow**:

	python Create_Blast_DB.py
	
This downloads the complete reference bacterial genomes from NCBI ([ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.fna.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.fna.tar.gz)), uncompresses the file, and extracts the genomes (concatentates different chromosomes, in case there are multiple strains it selects the shortest one, and excludes plasmids and transposons). It creates the directory **NCBI_DB** containing the file **BLAST_DB.fasta**, which looks like this:

	>Campylobacter_jejuni
	GTGCCAAAAGCGATATCAGTGAAGTCAGTGAACTGCATACACATATACATAAAGATGGAAAAATGATGAT
	GCAAAAAATTCCTGAAATTATTATAAAAGCTCATTCAAGCACAGAGCTTAAATCTGGAGGTTACCATATA
	...
	>Agrobacterium_radiobacter
	CCTAGTTGCCCCGAATCATAGAATCGGAGTACTTCTTATCTACGGAAAATGCCGGTAGTACCGACTTTCT
	GCGTAAATAAAGAAACCAGTCACAAAAACCCAACCATAAAGGCGCTTAGTAACGTGGCCAAAAATGGGAG
	...
It also contains the file **NCBI_Ref_Genome.txt** containing the lengths of these genomes

	Campylobacter_jejuni    1718980
	Agrobacterium_radiobacter       6656043
	...

### Note

If you want to add other reference genomes to this database, e.g., partial contigs, scaffolds, etc, you have to edit **both** of these files yourself. 

## Step 3

We construct a *BLAST* database for this file:

	~/ncbi-blast-2.3.0+/makeblastdb -in NCBI_DB/BLAST_DB.fasta -out NCBI_DB/BLAST_DB.fasta -dbtype nucl

## Step 4

We align our file **Tutorial/mock_sample.fna** with BLAST against this database. Remember to add the parameter **-outfmt 6** (tabular format). Change **num_threads 8* *to the number of threads you want to use.

	~/ncbi-blast-2.3.0+/bin/blastn -query Tutorial/mock_sample.fna -out Tutorial/mock_sample.blast -outfmt 6 -db ./MetaFlow_Blast/MetaFlow_BLAST_DB -num_threads 8

We should now have a file **Tutorial/mock_sample.blast** looking like this

	r1.1    Alteromonas_macleodii   99.592  245     1       0       1       245     1658204 1658448 1.63e-123       448
	r1.1    Alteromonas_sp._SN2     78.261  230     44      6       13      239     2469491 2469265 8.60e-32        143
	r2.1    Alteromonas_macleodii   96.957  230     7       0       1       230     2686617 2686846 3.45e-105       387
	...

## Step 5

From these BLAST alignments stored in **Tutorial/mock_sample.blast**, we construct the input LGF file for *MetaFlow*:

	python BLAST_TO_LGF.py Tutorial/mock_sample.blast NCBI_DB/NCBI_Ref_Genome.txt 250 1

where **250** is the average read length in our sample, and **1** is the sequencing machine (**0** for Illumina, **1** for 454 Pyrosequencing)

## Step 6

We now run *MetaFlow* 

	./metaflow -m Tutorial/mock_sample.blast.lgf -g NCBI_DB/NCBI_Ref_Genome.txt -c metaflow.config
	
The main output is **Tutorial/mock_sample.blast.lgf.abundance.csv**, which looks like this:

	Genome_Id       Genome_Name     Genome_Length   Num_Of_Chunks   Num_Of_Mapped_Reads     Absolute_Abundance      Relative_Abundance
	72      Streptococcus_mutans    2013587 1007    29483   3.61658 0.0118378
	161     Alteromonas_macleodii   4448980 2225    1716830 95.3156 0.311986
	198     Thermobispora_bispora   4189976 2095    65148   3.84049 0.0125707
	219     Lactobacillus_crispatus 2043161 1022    100361  12.1328 0.0397129
	...

where the last column contains the relative abundances.

## Step 7 (Drawing a circular cladogram)

We install [Nicola Segata](http://cibiocm.bitbucket.org)'s [GraPhlAn](https://bitbucket.org/nsegata/graphlan/src) in

	~/graphlan/

Then, we run

	python Drawing/draw.py -i Tutorial/mock_sample.blast.lgf.abundance.csv -g ~/graphlan/
	
and Violà!

![Example tree image](Drawing/tree_stool_sample.png)