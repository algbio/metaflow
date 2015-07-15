/*
 * MCFLogger.h
 *
 *  Created on: 21.1.2014
 *      Author: sobih
 */

#ifndef MCFLOGGER_H_
#define MCFLOGGER_H_

#include <string>
#include <fstream>
#include "time.h"

class MCFLogger{
	private:
		static string getTime(){
			struct tm *aTime ;
			time_t theTime = time(NULL);
			aTime = localtime(&theTime);
			char timeBuffer[80];
			strftime(timeBuffer,80,"%G%m%d%H%M%S",aTime);
			return timeBuffer;
		}
		ofstream logFile;
	public:
		MCFLogger() {};
		// void init(string log_filename)
		// {
		// 	logFile.open((char*)(log_filename.c_str());
		// }
		void log(string msg){
			cout << msg << endl;
			// logFile << msg << endl;
		}
		void close(){
			// logFile.close();
		}
};
//MCFLogger mcfLogger;
#endif /* MCFLOGGER_H_ */
