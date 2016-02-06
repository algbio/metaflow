CC = g++
FILES = MCFConfig.h MCFFlowSolver.cpp MCFGenetic.cpp MCFMapper.cpp MCFUtils.cpp MCFUtils_Temp.cpp OptionParser.cpp
OUT_EXE = ../metaflow
PATH_TO_LEMON = ./lemon_binaries_linux/include/

CPPFLAGS = -O3 -I $(PATH_TO_LEMON)

all:
	cd Src; $(CC) $(CPPFLAGS) -o $(OUT_EXE) $(FILES)

clean:
	rm -rf $(OUT_EXE)