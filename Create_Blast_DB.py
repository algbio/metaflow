import sys
import os

if len(sys.argv)<2:
	print "Params: Blast bin directory path."
	sys.exit()
blastDir=sys.argv[1]


fnaDir = os.path.join(os.getcwd(), "NCBI_Bacteria_Fna")
dbDir  = os.path.join(os.getcwd(), "MetaFlow_Blast")
logsDir = os.path.join(os.getcwd(), "Logs") 
tarFile= "all.fna.tar.gz"


if not os.path.isdir(fnaDir):
	os.system("mkdir "+fnaDir)

if not os.path.isfile(os.path.join(fnaDir,tarFile)):
	os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/"+tarFile+" -P "+fnaDir)
	os.system("tar -C " + fnaDir + " -xvf " + os.path.join(fnaDir,tarFile))


blastDB= open("Blast_DB.fna", "w")

for dirname in os.listdir(fnaDir):
	dirname= os.path.join(fnaDir, dirname)
	if not os.path.isdir(dirname):
		continue
	for filename in os.listdir(dirname):
		if filename.endswith(".fna"):				
			f= open(os.path.join(dirname,filename), "r")
			header = f.readline()
			if ("plasmid" not in header.lower()) and ("transposon" not in header.lower()):
				splitLine=header.strip().split("|")[4]
				speciesName=splitLine.split(",")[0].strip().replace(" ","_")
				blastDB.write(">"+speciesName+" | "+splitLine+"\n")
				for line in f:
					blastDB.write(line)
			else:
				f.close()
blastDB.close()


if not os.path.isdir(dbDir):
	os.system("mkdir "+dbDir)

if not os.path.isdir(logsDir):
	os.system("mkdir "+logsDir)

errLog=os.path.join(logsDir,"blastDB.err")
blastLog=os.path.join(logsDir,"blastDB.out")
os.system(blastDir+"/makeblastdb -in Blast_DB.fna -title MetaFlow_BLAST_DB -out "+dbDir+"/MetaFlow_BLAST_DB -dbtype nucl > "+blastLog+" 2> "+errLog)

if os.stat(errLog).st_size < 50:
	os.system("rm "+os.path.join(fnaDir,tarFile))
	os.system("rm -r "+fnaDir)
	os.system("rm Blast_DB.fna")
	print("Blast database was successfully created.")

else:
	print os.stat(errLog).st_size
