import sys
import os

# if len(sys.argv)<2:
# 	print "Params: BLAST bin directory path."
# 	sys.exit()
# blastDir=sys.argv[1]


fnaDir = os.path.join(os.getcwd(), "NCBI_Bacteria_Fna")
dbDir  = os.path.join(os.getcwd(), "NCBI_DB")
logsDir = os.path.join(os.getcwd(), "Logs") 
tarFile= "all.fna.tar.gz"
genomeFile = "NCBI_Ref_Genome.txt"
fastaFile = "BLAST_DB.fasta"


if not os.path.isdir(fnaDir):
	os.system("mkdir " + fnaDir)
if not os.path.isdir(dbDir):
	os.system("mkdir " + dbDir)

print "Downloading NCBI database..."
if not os.path.isfile(os.path.join(fnaDir,tarFile)):
	os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/"+tarFile+" -P "+fnaDir)
	print "Extracting genomes..."
 	os.system("tar -C " + fnaDir + " -xvf " + os.path.join(fnaDir,tarFile))


print "Preparing BLAST database..."
print "Parsing the genomes..."

genomes_concatenate_chromosomes = dict()
genomes_minimum_length_strain = dict()
genome_length_concatenate_chromosomes = dict()
genome_length_minimum_length_strain = dict()

for dirname in os.listdir(fnaDir):
	dirname= os.path.join(fnaDir, dirname)
	if not os.path.isdir(dirname):
		continue
	for filename in os.listdir(dirname):
		if filename.endswith(".fna"):				
			f= open(os.path.join(dirname,filename), "r")
			header = f.readline()
			# print header
			if ("plasmid" not in header.lower()) and ("transposon" not in header.lower()):

				# normalizing the species name
				splitLine=header.strip().split("|")[4]
				splitLine2=splitLine.strip().split(" ")
				speciesName = splitLine2[0] + "_" + splitLine2[1]
				if splitLine2[1] == "sp." or splitLine2[1] == "genomosp." or splitLine2[1] == "cf.":
					speciesName += "_" + splitLine2[2].strip(",")
				speciesName = speciesName.strip(".")
				speciesName = speciesName.strip(",")
				speciesName = speciesName.replace("(","")
				speciesName = speciesName.replace("'","")
				# speciesName = speciesName.replace(".","")
				# speciesName = speciesName.replace("-","")
				if speciesName.count("_") >= 3:
					n = 3
					groups = speciesName.split('_')
					speciesName = '_'.join(groups[:n]) + "-" + '-'.join(groups[n:])

				if "complete genome" in header.lower():
					sequence = ""
					length = 0
					for line in f:
						sequence += line
						length += len(line) - 1  # -1 for the \n character
					if (speciesName not in genomes_minimum_length_strain) or (len(sequence) < genomes_minimum_length_strain[speciesName]):
						genomes_minimum_length_strain[speciesName] = sequence
						genome_length_minimum_length_strain[speciesName] = length

				if "complete sequence" in header.lower():
					if (speciesName not in genomes_concatenate_chromosomes):
						genomes_concatenate_chromosomes[speciesName] = ""
						genome_length_concatenate_chromosomes[speciesName] = 0
					sequence = ""
					length = 0
					for line in f:
						sequence += line
						length += len(line) - 1  # -1 for the \n character
					genomes_concatenate_chromosomes[speciesName] += sequence
					genome_length_concatenate_chromosomes[speciesName] += length
			f.close()	

print "Writing the genomes to fasta file..."

blastDB= open(os.path.join(dbDir,fastaFile), "w")
genomeList= open(os.path.join(dbDir,genomeFile), "w")

for key in genomes_concatenate_chromosomes:
	if key in genomes_minimum_length_strain:
		if genome_length_minimum_length_strain[key] <= genome_length_concatenate_chromosomes[key]:
			# print "Skipping duplicate genome " + key
			continue
	genomeList.write(key + "\t" + str(genome_length_concatenate_chromosomes[key]) + "\n")
	blastDB.write(">" + key + "\n")
	blastDB.write(genomes_concatenate_chromosomes[key])

for key in genomes_minimum_length_strain:
	if key in genomes_concatenate_chromosomes:
		if genome_length_minimum_length_strain[key] > genome_length_concatenate_chromosomes[key]:
			# print "Skipping duplicate genome " + key
			continue
	genomeList.write(key + "\t" + str(genome_length_minimum_length_strain[key]) + "\n")
	blastDB.write(">" + key + "\n")
	blastDB.write(genomes_minimum_length_strain[key])

blastDB.close()
genomeList.close()

# if not os.path.isdir(dbDir):
# 	os.system("mkdir "+dbDir)

# if not os.path.isdir(logsDir):
# 	os.system("mkdir "+logsDir)

# print "Creating BLAST database"
# errLog=os.path.join(logsDir,"blastDB.err")
# blastLog=os.path.join(logsDir,"blastDB.out")
# os.system(os.path.join(blastDir,"makeblastdb") + " -in BLAST_DB.fna -title BLAST_DB -out "+dbDir+"/NCBI_DB -dbtype nucl > "+blastLog+" 2> "+errLog)

os.system("rm -r "+fnaDir)

print "Done"

# if os.stat(errLog).st_size < 20:
# 	os.system("rm -r "+fnaDir)
# 	os.system("rm -r "+logsDir)
# 	# os.system("rm BLAST_DB.fna")
# 	print("BLAST database successfully created")

# else:
# 	print os.stat(errLog).st_size
