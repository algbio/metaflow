/*
 * MCFUtil_Temp.cpp
 *
 *  Created on: Jan 26, 2014
 *      Author: ahmedsobih
 */

#include "MCFHeaders.h"

string getStringValue(int number){
	return static_cast<ostringstream*>( &(ostringstream() << number) )->str();
}
string getStringValue(double number){
	return static_cast<ostringstream*>( &(ostringstream() << number) )->str();
}
void split(vector<string>& splittedLine, string s, string separator){
	splittedLine.clear();
	int pos=s.find(separator);
	int start=0;
	while(pos!=string::npos){
		splittedLine.push_back(s.substr(start, pos-start));
		start=pos+separator.size();
		pos=s.find(separator, start);
	}
	splittedLine.push_back(s.substr(start, pos));
}
int64_t absolute(int64_t val){
	if (val<0)
		return -val;
	return val;
}

