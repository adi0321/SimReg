#include <vector>
#include <string>
#include <map>

using namespace std;

/********************************************************************
 * Pre-processing Step: Extract MC Read Classes from Bowtie output. *
 ********************************************************************/
void prepro_extract_read_classes_cc(
							string &rBowtie_file, 
							map<vector<string>, int> &ReadsClasses,
							map<string, string> &readsRef,
							map<vector<string>, map<string, int> > &readsClassesRef);


/**********************************************************
 * Extract Reads References from Grinder input. 
 **********************************************************/
void extract_reads_references(string &rGrinder_file, map<string, string> &readRef);

