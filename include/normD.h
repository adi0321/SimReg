#include <vector>
#include <string>
#include <map>

using namespace std;

/********************************************************************
 * Normalize D matrix by column (by transcript) *
 ********************************************************************/
void normD(map< string, map< vector<string>, double> > &d_values, map<string, int> &sumTrCounts);


