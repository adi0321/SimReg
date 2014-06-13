#include <vector>
#include <string>
#include <map>

using namespace std;

/********************************************************************************
*	Prints vector of strings													*
*		Example: used for printing read classes mcreads_class in MCReg_v2.cpp	*
*********************************************************************************/
void print(const vector<string> &myVector);

void print(const map<vector<string>, int> &myMap);//this prints reads classes and sizes
void print(const map<vector<string>, double> &myMap);//this prints normalized observed class frequency

/********************************************************************************
*		Example: used for printing d_values in MCReg_vi.cpp	*
*********************************************************************************/
void print(const map<string, map<vector<string>,double> > &myMap);

/********************************************************************
 *      Prints d values	- new - with the component id	           	*
 ********************************************************************/
void print(const map<int, map<string, map<vector<string>,double> > > &myMap);

/********************************************************************************
 *      Prints ALL Normalized adjusted transcripts frequency a line           	*
 * 		-from position_transcript					*
 *******************************************************************************/
void print(const map<string, map<string, pair<double, double> > > &myMap);


/****************************************************
 *      Prints ALL READS CLASSES		           	*
 ****************************************************/
void print_reads_classes(const map< map<string, bool> , int> &myMap);


/****************************************
 *      Prints d values		           	*
 ****************************************/
void print(const map<string, map< map<string, bool> , double> > &myMap);

/****************************************
 *      Prints current read            	*
 ****************************************/
void print(const map<string, bool > & myMap);


/****************************************
 *      Prints reads references        	*
 ****************************************/
void print(const map<string, string> & myMap);

/****************************************
 *      Print from where each read came from each class       	*
 ****************************************/
void print(const map< map<string, bool> , map<string,int> > &myMap);

//The new function:
void print(const map< vector<string> , map<string,int> > &myMap);


/****************************************
 *      Print observed/virtual reads Classes     	*
 ****************************************/
void print(const map< map<string, bool> , int> &myMap);


