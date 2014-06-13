#include <vector>
#include <string>
#include <map>

using namespace std;

/********************************************************************
 * Pre-processing Step: Extract Observed Read Classes from Bowtie output. *
 ********************************************************************/
int extract_observed_rc_cc(string &rBowtie_file, map<vector<string>, int> &ReadsClasses);