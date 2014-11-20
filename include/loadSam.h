#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <sstream>      // std::istringstream
#include <string>
#include <map>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;

/********************************************************************
 * Load SAM File *
 ********************************************************************/
void loadSam(string &inFile, map<string, map<string, map<int, string>>> &sam, bool flagMAQC);


