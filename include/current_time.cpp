#include "current_time.h"
#include <string>

// ostringstream constructor
#include <iostream>     // std::cout, std::ios
#include <sstream>      // std::ostringstream

using namespace std;

string current_time(void){
	time_t rawtime;
	time(&rawtime); //get current time
	ostringstream sTime;
	sTime << ctime(&rawtime);
	string sTimeS = sTime.str();
	sTimeS = sTimeS.substr(0,sTimeS.length()-1);

return sTimeS;
}

