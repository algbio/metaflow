#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <sstream>

using namespace std;

class Logger {

public:
	ofstream writer;
	Logger(){
		time_t theTime = time(NULL);
		struct tm *aTime = localtime(&theTime);
		stringstream fileName;
		fileName << "Logs/MCF_Log_"<<aTime->tm_year+1900<<aTime->tm_mon+1<<aTime->tm_mday<<aTime->tm_hour<<aTime->tm_min<<".log";
		writer.open(fileName.str().c_str());
	}
	void log(string log){
		writer << log<<endl;
		writer.flush();
	}
	void close(){
		writer.close();
	}
};
extern Logger logger;

