#include <vector>
#include <string>
#include <map>

using namespace std;

/**********************************************************
 * Normalize values. 
 **********************************************************/

void norm(vector<float> &a);
void norm(vector<double> &a);


/**********************************************************
 * 	Used for normalizing: genes_frequency             *
 **********************************************************/
void norm(map<string, double> &a);


/***************************************************************************************************
 *      Used for normalizing: ADJUSTED TRUE transcript frequency: adjusted_true_tr_freq            *
 ***************************************************************************************************/
void norm(map<string, map<string, double> > &a);

