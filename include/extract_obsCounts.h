#include <vector>
#include <string>
#include <map>

using namespace std;

/********************************************************************
 * Pre-processing Step: Extract Observed Read Classes from Bowtie output. *
 ********************************************************************/
int extract_obsCounts(string &rcCounts_file, map<vector<string>, int> &ReadsClasses);