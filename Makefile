CC = g++
FILES = MCFConfig.h MCFFlowSolver.cpp MCFGenetic.cpp MCFMapper.cpp MCFUtils.cpp MCFUtils_Temp.cpp
OUT_EXE = metaflow
PATH_TO_LEMON = ./lemon_binaries_linux/include/

CPPFLAGS = -O3 -I $(PATH_TO_LEMON)

all:
	$(CC) $(CPPFLAGS) -o $(OUT_EXE) $(FILES)

example: 
	

clean:
	rm -rf $(OUT_EXE)