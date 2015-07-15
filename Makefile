CC = g++
FILES = MCFConfig.h MCFFlowSolver.cpp MCFGenetic.cpp MCFMapper.cpp MCFUtils.cpp MCFUtils_Temp.cpp OptionParser.cpp
OUT_EXE = metaflow
PATH_TO_LEMON = ./lemon_binaries_linux/include/

CPPFLAGS = -O3 -I $(PATH_TO_LEMON)

all:
	$(CC) $(CPPFLAGS) -o $(OUT_EXE) $(FILES)

example: 
	mkdir -p Example; cd Example; wget http://cs.helsinki.fi/u/tomescu/metaflow/example/MCF_Sample_100.blast.gz; gunzip MCF_Sample_100.blast.gz; 

clean:
	rm -rf $(OUT_EXE)