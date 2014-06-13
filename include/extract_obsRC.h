#include <vector>
#include <string>
#include <map>

using namespace std;

/********************************************************************
 * Pre-processing Step: Extract OBS Read Classes from Bowtie output. *
 ********************************************************************/
void extract_obsRC(string &rBowtie_file, map<vector<string>, int> &ReadsClasses);