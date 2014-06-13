#include "normD.h"
#include <vector>
#include <string>
#include <map>

#define DEBUG 1

/********************************************************************
 * Normalize D matrix by column (by transcript) *
 ********************************************************************/
 
void normD(map< string, map< vector<string>, double> > &d_values, map<string, int> &sumTrCounts)
{

//for each transcript (column) in D matrix
for(const auto &tr: d_values){		
			
	//compute total sum
	double total_sum=0.0;
			
	//for each read class in the current transcript
	for (const auto &it : tr.second)
		total_sum+=it.second;
			
	sumTrCounts[tr.first]=total_sum;
	//Normalize
	//again for each read class in the current transcript
	for (const auto &it : tr.second)
	{
		if(total_sum==0)
			d_values[tr.first][it.first]=0;
		else
			d_values[tr.first][it.first]=it.second/total_sum;				
	}
			
}//end: for each transcript (column)
}//end normD


