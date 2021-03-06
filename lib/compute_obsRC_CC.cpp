//============================================================================
// Project     : SimReg
// Author      : Adrian Caciula
// Description : calculate d for the bipartite graph(transcripts and reads)
// Created     : 06/01/2013
// Last Update: 09/01/2014
// Compile command:
// g++ -std=c++0x ./compute_obsRC_CC.cpp ../include/current_time.cpp ../include/norm.cpp ../include/print.cpp ../include/extract_obsRC.cpp -o ./compute_obsRC_CC -I /usr/local/include/ -L /usr/local/lib/
//============================================================================

//This is an extension of version of ~/code/MCReg/annotation_preprocessing/compute_read_classes_v2.cpp (and compute_RC_CC_v2.cpp)
//This version computes Connected Components in addition to MC read classes and d values
//2013.12.16 - Compute Observed RC and CC
	
//Note (8/28/2014): There is a problem with the read classes. We should only make a decision after all reads have been parsed...but this means to store in memory all read names :( ... we loose the quality ... of the read

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//TO DO: 9/3/2014
//start with the id of the last component
//if it doesn't belong to any other class then don't delete it if this transcript belong to other classes
//A transcript is a leaf if it belongs to a single class. Delete delete but remember which is a leaf (add letter L)
//t1t5t8t10 (t8 and t10 are leafs...return the names...and keep the)
//kernel and leaves.
//Keep Sub class 	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "../include/current_time.h"
#include "../include/norm.h"
#include "../include/print.h"
#include "../include/extract_obsRC.h"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>

#include <map>
#include <iterator>
#include <random>
#include <iomanip>      // std::setprecision

//Packages for boost:
#include <boost/config.hpp>
#include <iostream>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/labeled_graph.hpp>


//DEBUGGING CONSTANTS:
#define DEBUG 0
#define DEBUG_READ 0 //debug at read level
#define DEBUG_GF 0 //debug for genes frequency
#define DEBUG_T 0 //debug Transcripts
#define DEBUG_GEN 0 //Debug the read generator
#define DEBUG_GEN_S 0 //debug for simulated reads
#define DEBUG_GEN_O 0 //debug for observed reads
#define DEBUG_T_FREQ 0 //debug transcript frequency

			
#define DEBUG_O 0

using namespace std;
using namespace boost;

#define HELPMESSAGE "too few parameters: \n\
Input:\n\
argv[1] - File with mapped Observed reads (from Bowtie using options --best -v 3 -k 60) \n\
argv[2] - GTF File with exons coordinates \n\
argv[3] - FA File with transcript sequences \n\
"
		
void print(const vector<pair<int,int> > &myVector);
void print(const map<string, int> &myMap);
void print(const map<string, double> &myMap);
void print(const map<string, map<string, int> > &myMap); //Print Transcripts Length
void print(const map<string, map<string, double> > &myMap); //Print Normalized Adjusted Transcripts Length
void print(const map<string, vector<pair<int, int> > > &myMap);
void print(const map<string, map<string, vector<pair<int, int> > > > &myMap); //print all exons
void print(const map<string, vector<int> > &myMap);
void print(const map<string, map<string, vector<int> > >&myMap);
void print(const map<vector<int>, int > &myMap);
void print(const map<vector<int> , vector<pair<int, int> > > &myMap);
void print(const map< map<string, bool> , double> &myMap);
void print(const map<string, map< map<string, bool> , int> > &myMap);

void splitCmp (int cmpID, vector<string> trVec, map<int, vector<string>> &rcCounts, map<string, vector<string>> &rcID_reverse, map<vector<string>, int> &obsReadsClasses, map<vector<string>, string> &rcID, map<string, map<string, int>> &gtf, map <string, string> &fa);

void loadGTF(string &gtf_file, map<string, map<string, int>> &gtf);
void loadFA(string &fa_file, map<string, string> &fa);

void write2Files(int cmpID, vector<string> transcripts, map<string, map<string, int>> &gtf, map <string, string> &fa);

int readCountThresh=1; //initial read count threshold (counts used to collapse read classes)
int cmpSizeThresh=100; //initial number of transcripts that is allowed in one component
int subCmpId=100000;

int main(int argc,char *argv[]){

//cout<<"\nRunning "<<argv[0]<<endl;

	if(argc<4){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

double EPS = 0.01;
int iterations=0;

	string obsRBowtie_file=argv[1]; //aligned Observed reads
	string gtf_file=argv[2];
	string fa_file=argv[3];
	
clock_t begin = clock(); //used for measuring entire elapsed time for this function
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Load GTF "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		map<string, map<string, int>> gtf;
		loadGTF(gtf_file, gtf);
		
		//Write to file transcript lengths //This file is used to replace the use of grep for computing transcript lenght (grep is slow)
		ofstream outGlobalTrLength;
		string globalTrLengthFile="../trLen.txt";
		outGlobalTrLength.open(globalTrLengthFile.c_str());
	
		if(!outGlobalTrLength){
			cout<<"Unable to open " <<globalTrLengthFile<<endl;
			exit(1);
		}
		
		for (const auto &tr : gtf){
			outGlobalTrLength<<tr.first<<"\t"<<tr.second.begin()->second<<endl;
		}
		outGlobalTrLength.close();
		
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;	

cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Load FA File"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		map<string, string> fa;
		loadFA(fa_file, fa);
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;	
		
		
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Compute Observed Reads Classes "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		
//Extract OBS reads classes from Bowtie output
map<vector<string>, int> obsReadsClasses; 

cout<<"\nMain: Parsing Bowtie file: "<<endl;

extract_obsRC(obsRBowtie_file, obsReadsClasses);

//~~~~~~~~~~~~~~~~~~~~~~~
//08/26/2014~~~~~~~~~~~~~
//delete classes with count 10

//for each read class
//for (const auto &r_class : obsReadsClasses){
//	if(r_class.second<10)
//		obsReadsClasses.erase(r_class.first);
//}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~


//Create an ID for each Read Class --- this is just a current solution --- needs improvement --- not smart because will double the memory
map<vector<string>, string> rcID;
map<string, vector<string>> rcID_reverse; //in this way we have a bidirectional map (which is stupid anyway...but what to do?)
int total_obs_reads=0;
int id=0;
//for each readClass in obsReadsClasses add an id and Compute Total number of Observed reads
for (const auto &r_class : obsReadsClasses){
	stringstream cID_ss;
	cID_ss<<"["<<id<<"]";//these square brackets were inserted in order to avoid modifying the code down when it checks for "["

	rcID[r_class.first]=cID_ss.str();
	rcID_reverse[cID_ss.str()]=r_class.first;
	id++;
	total_obs_reads+=r_class.second;
}

cout<<"\n~~~ Total number of Observed reads = "<<total_obs_reads<<endl;
cout<<"Total number of Observed Reads classes is: "<<obsReadsClasses.size()<<endl;		

	#if DEBUG
		cout<<"\nOBS Reads Classes \t Size:"<<endl;
			//print(obsReadsClasses);
			cout<<endl;
			//exit(7);
	

	

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//This might be only for debugging
		//07/30/2014: Print all observed classes to file 
		ofstream allRCout; // all observed read classes
		string allRCFile="../allObsReadClasses.txt";
		allRCout.open(allRCFile.c_str(),ios::app);//open file in append mode (in order to avoid overwritting)
			
		if(!allRCout){
			cout<<"Unable to open " <<allRCFile<<endl;
			exit(1);
		}

			for (const auto &read_class: obsReadsClasses){ //for each read class
					allRCout<<"[ ";
					for (const auto &tr : read_class.first)//for each transcript in this class
						allRCout<<tr<<" ";
					
					allRCout<<"] "<<read_class.second<<endl;
				
				
			}


		allRCout.close();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		

	#endif

		
		
	
//cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
//cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;	

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cout<<"\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Computing adjacency matrix between reads (classes) and transcripts... "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

//We'll use d_value (p_value) to hold the adjacency matrix
map< string, map< vector<string>, double> > p_value_new2; //this is d_value


		//for each OBS Read class
		for (const auto &read_class: obsReadsClasses){
			
			//for each transcript in the current class
			for (const auto &tr: read_class.first){
				
				//if this value exists --> then we have an error
				if(p_value_new2[tr][read_class.first]){
					cout<<"Error: p value new2 already exists!!!"<<endl;
					exit(1);
				}
				
				p_value_new2[tr][read_class.first]=read_class.second;//enter the counts for each class
			}
		}
		
				#if DEBUG
					cout<<"\nPrint p_value_new2"<<endl;
					//print(p_value_new2);
					
					//exit(7);
				#endif


//cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//cout<<"["<<current_time()<<"] ~~~~~ Done: adjacency matrix has been computed"<<endl;
//cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;				

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/		

cout<<"\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Computing Connected Components ... "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

//Both DFS and BFS complete on O(n) time, using O(n) memory, where n is matrix size. But BFS it doesn't suffer from stack overflow problem, and it doesn't spend time on recursive calls.
//Anyway, let's first use boost library

//*********************************************For debug only*******************************//////////
//08/24/2014: Print all components and edges to file 
/*
	ofstream allCout; // all components
	string allCFile="../allComponents.txt";
	allCout.open(allCFile.c_str(),ios::app);//open file in append mode (in order to avoid overwritting)
	
	if(!allCout){
		cout<<"Unable to open " <<allCFile<<endl;
		exit(1);
	}
*/
//*****************************************////////////////////////////////////////////////////////////


	struct Vertex {
		std::string vertexName;
	};
	
	//typedef adjacency_list <vecS, vecS, undirectedS, Vertex> MyGraph;
	typedef boost::labeled_graph<adjacency_list<vecS, vecS, undirectedS, Vertex>, std::string> MyGraph;
	//The side effect of this is that one need to use graph() member function to make some algorithms work:
	
	typedef boost::graph_traits<MyGraph>::vertex_descriptor VertexID;
	
    MyGraph G;
	
	//for each transcript --  
	for (const auto &tr : p_value_new2){
		// Create vertices in that graph
		VertexID u = boost::add_vertex(tr.first,G);
		G[tr.first].vertexName = tr.first;
		
		//for each class in the current transcript
		for (const auto &read_class : tr.second ){
			
			//string rClassName = "["; //add [ in order to differentiate the classes that contains only one transcripts from the transcripts itself
			string rClassName=rcID[read_class.first]; //get just the read class ID
			//---Lines below have been replaced by line above
			//for each transcript in this class -- create a string with the name of the class
			//for (const auto &tr_class : read_class.first)
			//	rClassName+=tr_class;
			
			//cout<<"Class name is: "<<rClassName<<endl;
			
			VertexID v = boost::add_vertex(rClassName,G);
			G[rClassName].vertexName = rClassName; //in case vertex already exists then this should overwrite
			
			add_edge(u, v, G);
			//allCout<<"[Print 1] Edge btwn vertex "<<tr.first<<" and vertex "<<rClassName<<endl;
		}
	}
	
	
    std::vector<int> component(num_vertices(G));
	
    int num = connected_components(G, &component[0]);
    cout << "Total number of components: " << num << endl;
	
	/*

	std::vector<int>::size_type j;
	for (j = 0; j != component.size(); ++j)
      allCout << "Vertex " << j <<" is in component " << component[j]<<endl;
	  
    allCout << endl;
	
	

	cout<<"\nPrint components: "<<endl;
	for (const auto &c : component) {
		cout<<"Component: "<<c<<endl;

	}
	*/
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//08/22/2014: Print all components to file 
/*	
	MyGraph::vertex_iterator vertexIt, vertexEnd;
	boost::tie(vertexIt, vertexEnd) = vertices(G);
	for (; vertexIt != vertexEnd; ++vertexIt){
		VertexID vertexID = *vertexIt; // dereference vertexIt, get the ID
		Vertex & vertex = G.graph()[vertexID];
		//The side effect of boost::labeled_graph is that one need to use graph() member function to make some algorithms work
	
		allCout<<"[Print 2] Vertex name is : "<<vertex.vertexName<<endl;
    }

	allCout.close();
*/
//cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
//cout<<"["<<current_time()<<"] ~~~~~ Done: Connected Components had been computed"<<endl;
//cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
	
		
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~ Write Values to Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cout<<endl<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"Prepare Values for SimReg"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"Please wait, this may take few minutes"<<endl;
	cout<<"writing to files . . ."<<endl;
/*********************************
*********************************/
/*
//Write OBS Read Classes Names and Frequencies to File 
string read_classes_file2="obsCC.txt";

cout<<"\nWrite read classes to: \n\t"<<read_classes_file2<<endl;

ofstream read_classes_stream;
read_classes_stream.open(read_classes_file2.c_str());
		
if(!read_classes_stream){
	cout<<"Unable to open" <<read_classes_file2<<endl;
		exit(1);
}

//for each read class
for (const auto &read_class : obsReadsClasses) {

		//print to file the class name
		read_classes_stream<<"[ ";
		for (const auto &cluster : read_class.first){
				read_classes_stream<<cluster<<" ";
		}
		read_classes_stream<<"]\t"<<new_o[read_class.first]<<endl;
}
read_classes_stream.close();
/*********************************************************************************
**********************************************************************************/



/****************************
*	Write transcripts from each component to file	*
*****************************/

	
	//go thorugh components
	//cout<<"component size is: "<<component.size()<<endl;
	//CREATE MAP THAT WILL HOLD ALL THE COMPONENTS (transcript names):
	map<int, vector<string>> cc;
	
	//CREATE MAP THAT WILL HOLD ALL READ CLASSES AND COUNTS //   ### 9/2/2014: Not sure about the counts? The integer looks like the component ID?
	
	map<int, vector<string>> rcCounts;
	
	vector<int>::size_type i=0;
	Vertex & vertex = G.graph()[i];
	//cout<<"Vertex "<<i<<" has name: "<<vertex.vertexName<<endl;
	
	//Notes for each component (# of components == # of vertices but the component id changes only when we have a new component)
	for (i = 0; i < component.size(); i++){
		//cout<<"Component iterations: "<<i<<endl;
	
		Vertex & vertex = G.graph()[i];
		//cout<<"i="<<i<<" Element from component: "<<component[i]<<" is "<<vertex.vertexName<<endl;
		
		if(vertex.vertexName[0] != '[')//if the first character is not [
			cc[component[i]].push_back(vertex.vertexName);
		else
		{ //if the first character it is a "[" then this a read class
			rcCounts[component[i]].push_back(vertex.vertexName);
		}
	}//end: for each component
	
		#if DEBUG
			cout<<"Print component map:"<<endl;
			for(const auto &cmp : cc)
			{
				cout<<"component "<<cmp.first<<" - size: "<<cmp.second.size()<<" containts: "<<endl;
					for(const auto &elem : cmp.second)
						cout<<"\t\t"<<elem<<endl;
			}
		#endif


//Now print to file
cout<<"\nWrite observed transcripts names\n\t"<<endl;

ofstream resultStream;
string resultsFile="../singleTrGenes.txt";
resultStream.open(resultsFile.c_str(),ios::app);//open file in append mode (in order to avoid overwritting)

	if(!resultStream){
		cout<<"Unable to open " <<resultsFile<<endl;
		exit(1);
	}

//For each component
for(const auto &cmp : cc)	
{//for each component there is a vector with transcript names


//Check the size of the component and if it is 1 then just write the results to file
if( cmp.second.size()< 2 )
{
	resultStream<<cmp.first<<"\t"<<cmp.second.size()<<"\t"<<obsReadsClasses[cmp.second];
	resultStream<<"\t"<<cmp.second.front()<<"\t1"<<endl;
	//cmpID---#tr.---#ObsReads---Transcript names

	continue; //continue to the next component
}

//~~~~~~~
//8/28/2014
if(cmp.second.size() > cmpSizeThresh)
{

	//Prepare the prefix for sub-components
	//stringstream ss;
	//ss << cmp.first;
	//string dirPrefix=ss.str();
		
	//Recursive Function
	cout<<"Call the recursive function"<<endl;
	
	//Reset the thereshold value everytime this function is called
	readCountThresh=1;
	
	splitCmp(cmp.first, cmp.second, rcCounts, rcID_reverse, obsReadsClasses, rcID, gtf, fa);
	continue;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//#########################################################################

write2Files(cmp.first, cmp.second, gtf, fa);

//#########################################################################
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ofstream d_stream_new2;
//d_stream_new2.clear(); //reuse the same stream (just clear the state flags)

stringstream outFile;
//outFile.str(""); //clear the stringstream
outFile<<"./"<<cmp.first<<"/obsRCcounts.txt"; //directory already exists from previous 
string d_file_new2=outFile.str();

d_stream_new2.open(d_file_new2.c_str());
	if(!d_stream_new2){
		cout<<"Unable to open" <<d_file_new2<<endl;
		exit(1);
	}
	
	//for each read class ID
	for (const auto &class_ID : rcCounts[cmp.first]){
		d_stream_new2<<"[\t";
		
		//for each transcript in this read class rcID_reverse[class_ID]-->print the class and the counts
		for (const auto &tr : rcID_reverse[class_ID])
		{
			d_stream_new2<<tr<<"\t";
		}
		d_stream_new2<<"]"<<"\t"<<obsReadsClasses[rcID_reverse[class_ID]]<<endl;
	}//end: for each read class ID	

d_stream_new2.close();
	

//Write to file the number of transcripts in each readClass
d_stream_new2.clear(); //reuse the same stream (just clear the state flags)
outFile.str(""); //clear the stringstream
outFile<<"./"<<cmp.first<<"/obsRCsize.txt"; //directory already exists from previous 
d_file_new2=outFile.str();

d_stream_new2.open(d_file_new2.c_str());
	if(!d_stream_new2){
		cout<<"Unable to open" <<d_file_new2<<endl;
		exit(1);
	}
	
	
	//Matrix for size of read class
	vector<vector<double>> rcSize;
	//initialize matrix rcSize with zero
	
	for(int i=0;i<rcCounts[cmp.first].size();i++)
    {
		vector<double> row; // Create an empty row
         for(int j=0;j<rcCounts[cmp.first].size();j++)
         {
              row.push_back(0.0);   // Add an element (column) to the row
         }
		rcSize.push_back(row); // Add the row to the main vector		 
	}	
	
	/*
	cout<<"print(rcSize)"<<endl;
	for(int i=0;i<rcSize.size();i++)
    {
         for(int j=0;j<rcSize.size();j++)
         {
              cout<<rcSize[i][j];   
         }
		cout<<endl;
	}	
	exit(7);*/
	
	//fill the diagonal with read classes size.
	//for each read class ID
	int i=0;
	for (const auto &class_ID : rcCounts[cmp.first])
	{
		rcSize[i][i]=(double)1/rcID_reverse[class_ID].size();
		i++;
	}
	
	//print Matrix to file
	for(int i=0;i<rcCounts[cmp.first].size();i++)
    {
         for(int j=0;j<rcCounts[cmp.first].size();j++)
         {
				d_stream_new2<<rcSize[i][j]<<"\t";
         }   
         d_stream_new2<<endl;
	}	

d_stream_new2.close();
}//end:for each cc	


resultStream.close();

		
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
			
clock_t end = clock();
double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cout<<"\nDone - Elapsed time: "<<elapsed_secs<<endl;
}//end main





/*
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~ FUNCTIONS DEFINITION ~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*/


//~~~~~~~
//8/29/2014
//Recursive function
//Keep Spliting until component gets smaller than let's say 300
void splitCmp (int cmpID, vector<string> trVec, map<int, vector<string>> &rcCounts, map<string, vector<string>> &rcID_reverse, map<vector<string>, int> &obsReadsClasses, map<vector<string>, string> &rcID, map<string, map<string, int>> &gtf, map <string, string> &fa)
{

	//Note: Leave single transcript classes "as is"
	
	int tempCount=0;
	
		
		//cout<<"Component "<<cmpID<<" has size > 10"<<endl;
	
		//sort the current read class in another vector
		vector<string> sortedReadClass;
							

		//Adjust read classes in this component and resplit using DFS (or boost)
		//1. Get the read classes for this component and save them in a local variable
			//It's faster to search for matching read class in this small structure
			map<vector<string>, int> cmpObsReadsClasses; //All obs read class that belong to this component only
			//for each read class ID
			for (const auto &class_ID : rcCounts[cmpID]){
			
				sortedReadClass.clear();
				//for each transcript in this read class rcID_reverse[class_ID]-->print the class and the counts
				for (const auto &tr : rcID_reverse[class_ID])
				{
					sortedReadClass.push_back(tr);
					//cout<<tr<<"\t";
				}
				//cout<<"]"<<"\t"<<obsReadsClasses[rcID_reverse[class_ID]]<<endl;
				sort(sortedReadClass.begin(), sortedReadClass.end());
				cmpObsReadsClasses[sortedReadClass]=obsReadsClasses[rcID_reverse[class_ID]];
			}
			
						#if DEBUG
							cout<<"Print all initial read classes that belong to component"<<cmpID<<endl;
							print(cmpObsReadsClasses);
							//exit(7);
						#endif

						
		//2.Compress Read Classes
			//Delete the last transcript
			//for each read class ID
			for (const auto &class_ID : rcCounts[cmpID]){
				//cout<<"Processing Read Class: "<<endl;
			
				//cout<<"[\t";
		
				sortedReadClass.clear();
				//for each transcript in this read class rcID_reverse[class_ID]-->print the class and the counts
				for (const auto &tr : rcID_reverse[class_ID])
				{
					sortedReadClass.push_back(tr);
					//cout<<tr<<"\t";
				}
				//cout<<"]"<<"\t"<<obsReadsClasses[rcID_reverse[class_ID]]<<endl;
				
				sort(sortedReadClass.begin(), sortedReadClass.end());
				
				//if counts < 10 && readclass size is not 1
				if(obsReadsClasses[rcID_reverse[class_ID]] < readCountThresh && rcID_reverse[class_ID].size()>1)
				{
						
						//cout<<" Count of read class: "<<endl;
						//print(rcID_reverse[class_ID]);
						//cout<<" is : "<<obsReadsClasses[rcID_reverse[class_ID]]<<endl;
						
						//keep the initial class before removing elements from the end
						vector<string> initialReadClass;
							for (const auto &tr : rcID_reverse[class_ID])
								initialReadClass.push_back(tr);
							sort(initialReadClass.begin(), initialReadClass.end());
				
							
							bool foundClass=false;
							while(!foundClass && !rcID_reverse[class_ID].empty())
							{
									rcID_reverse[class_ID].pop_back();
								
								if(!rcID_reverse[class_ID].empty())
								{
										//if there is at least one more transcript in this class after removing 
									
										//sort again
										sortedReadClass.clear();
										for (const auto &tr : rcID_reverse[class_ID])
											sortedReadClass.push_back(tr);
									
										sort(sortedReadClass.begin(), sortedReadClass.end());
								
										//#if DEBUG
										//	cout<<"New Read Class "<<endl;
										//	cout<<" [";
										//	for (const auto tr : sortedReadClass)
										//		cout<<tr<<" ";
										//	cout<<" ]"<<endl;
										//#endif
								
										//no need to search inside all read class ... only those read classes that belong to this component
										//cmpObsReadsClasses have been created to address the issue above
										if ( cmpObsReadsClasses.find(sortedReadClass) != cmpObsReadsClasses.end() ){
											//if this read class is found among the given classes
												#if DEBUG
													cout<<" [";
													cout<<"Read class has been found"<<endl;
													for (const auto tr : sortedReadClass)
														cout<<tr<<" ";
													cout<<" ]"<<endl;
												#endif
											
											//increment the new founded class
											obsReadsClasses[sortedReadClass]+=obsReadsClasses[initialReadClass]; 	//global read classes
											cmpObsReadsClasses[sortedReadClass]+=obsReadsClasses[initialReadClass]; //local read classes
											
											//delete the old class
											obsReadsClasses.erase(initialReadClass);
											cmpObsReadsClasses.erase(initialReadClass);
											
											//Note: All structures that are related to obsReadsClasses MUST also be updated
												//This includes but not limited to: p_value_new2
												
											
											foundClass=true;
										}
								}
								else
								{ 
									//else after poping we got to the end of the stack-->this means no read class was found
									//drop the read and write this read to file the read to which was initially mapped that was initally mapped to initialReadClass
									
					/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
									//unused_reads<<read_name<<"\t mapped to class: [";
									cout<<"\t Initial class. (NOT Compatible with any other class): [";
									for (const auto tr : initialReadClass)
									cout<<tr<<" ";
									cout<<" ]"<<endl;
									//discarded_reads++;
					~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
					
									//delete the old class
									obsReadsClasses.erase(initialReadClass);
									cmpObsReadsClasses.erase(initialReadClass);
									
									#if DEBUG
										//cout<<"Read "<<read_name<<" was NOT mapped to any class"<<endl;
										cout<<"Initial Class: [";
										for (const auto tr : initialReadClass)
											cout<<tr<<" ";
										cout<<" ]"<<endl;
									#endif
									
								}
							} //end while class not found				
				}
				
			}//end: for each read class ID	
			
	//cout<<"Print adjusted read classes that belong to this component"<<endl;
	//print(cmpObsReadsClasses);	
		
		
		
	//2. Maybe create again the graph and run boost? (some recursion needs to be done here)
	//!!!the entire p_value_new must be updated?
		struct Vertex {
		std::string vertexName;
	};
	
	//typedef adjacency_list <vecS, vecS, undirectedS, Vertex> MyGraph;
	typedef boost::labeled_graph<adjacency_list<vecS, vecS, undirectedS, Vertex>, std::string> MyGraph;
	//The side effect of this is that one need to use graph() member function to make some algorithms work:
	
	typedef boost::graph_traits<MyGraph>::vertex_descriptor VertexID;
	
    MyGraph G;
	
	//We'll use d_value (p_value) to hold the adjacency matrix
	map< string, map< vector<string>, double> > p_value_new2; //this is d_value
	//local variable

		//for each OBS Read class
		for (const auto &read_class: cmpObsReadsClasses){
			
			//for each transcript in the current class
			for (const auto &tr: read_class.first){
				
				//if this value exists --> then we have an error
				if(p_value_new2[tr][read_class.first]){
					cout<<"Error: p value new2 already exists!!!"<<endl;
					exit(1);
				}
				
				p_value_new2[tr][read_class.first]=read_class.second;//enter the counts for each class
			}
		}
		
				#if DEBUG
					cout<<"\nPrint p_value_new2 (local):"<<endl;
					//print(p_value_new2);
					
					//exit(7);
				#endif

	//Build the graph

	//for each transcript --  
	for (const auto &tr : p_value_new2){
		// Create vertices in that graph
		VertexID u = boost::add_vertex(tr.first,G);
		G[tr.first].vertexName = tr.first;
		
		//for each class in the current transcript
		for (const auto &read_class : tr.second ){
			
			//string rClassName = "["; //add [ in order to differentiate the classes that contains only one transcripts from the transcripts itself
			string rClassName=rcID[read_class.first]; //get just the read class ID
			
			VertexID v = boost::add_vertex(rClassName,G);
			G[rClassName].vertexName = rClassName; //in case vertex already exists then this should overwrite
			
			add_edge(u, v, G);
			//allCout<<"[Print 1] Edge btwn vertex "<<tr.first<<" and vertex "<<rClassName<<endl;
		}
	}

	
    std::vector<int> component(num_vertices(G));
	
    int num = connected_components(G, &component[0]);
    //cout << "Total number of components: " << num << endl;
	
		//go thorugh components
	//cout<<"component size is: "<<component.size()<<endl;
	//CREATE MAP THAT WILL HOLD ALL THE COMPONENTS (transcript names):
	map<int, vector<string>> cc;
	
	//CREATE MAP THAT WILL HOLD ALL READ CLASSES AND COUNTS
	map<int, vector<string>> rcCountsLocal;
	
	vector<int>::size_type i=0;
	Vertex & vertex = G.graph()[i];
	//cout<<"Vertex "<<i<<" has name: "<<vertex.vertexName<<endl;
	
	//Notes for each component (# of components == # of vertices but the component id changes only when we have a new component)
	for (i = 0; i < component.size(); i++){
		//cout<<"Component iterations: "<<i<<endl;
	
		Vertex & vertex = G.graph()[i];
		//cout<<"i="<<i<<" Element from component: "<<component[i]<<" is "<<vertex.vertexName<<endl;
		
		if(vertex.vertexName[0] != '[')//if the first character is not [
			cc[component[i]].push_back(vertex.vertexName);
		else
		{ //if the first character it is a "[" then this a read class
			rcCountsLocal[component[i]].push_back(vertex.vertexName);
		}
	}//end: for each component
	
		#if DEBUG
			cout<<"Print component map:"<<endl;
			for(const auto &cmp : cc)
			{
				cout<<"component "<<cmp.first<<" - size: "<<cmp.second.size()<<" containts: "<<endl;
					for(const auto &elem : cmp.second)
						cout<<"\t\t"<<elem<<endl;
			}
		#endif
	

//Now print to file
//cout<<"\nWrite observed transcripts names\n\t"<<endl;

ofstream d_stream_new2;
ofstream resultStream;
string resultsFile="../singleTrGenes.txt";
resultStream.open(resultsFile.c_str(),ios::app);//open file in append mode (in order to avoid overwritting)
	
if(!resultStream){
	cout<<"Unable to open " <<resultsFile<<endl;
		exit(1);
}


//For each component
for(const auto &cmp : cc)	
{//for each component there is a vector with transcript names

//string glblCmpId=dirPrefix+"_"+cmp.first; //global component id
subCmpId++;

//Check the size of the component and if it is 1 then just write the results to file
if( cmp.second.size()< 2 )
{
	resultStream<<subCmpId<<"\t"<<cmp.second.size()<<"\t"<<obsReadsClasses[cmp.second];
	resultStream<<"\t"<<cmp.second.front()<<"\t1"<<endl;
	//cmpID---#tr.---#ObsReads---Transcript names
	
	continue; //continue to the next component
}

//~~~~~~~
//8/28/2014
if(cmp.second.size() > cmpSizeThresh)
{
	//Prepare directory prefix
	//stringstream ss;
	//ss << cmp.first;
	//dirPrefix+="_"+ss.str();
	
	//increment read count threshold (required for the next function call)
	readCountThresh*=2;
	//otherwise we'll get the same read classes (no new colapse will happen)
	
	//Recursive Function
	//cout<<"splitCmp: Call the recursive function"<<endl;
	splitCmp(cmp.first, cmp.second, rcCountsLocal, rcID_reverse, obsReadsClasses, rcID, gtf, fa);
	
	continue;
}

//~~~~~~~
//Else write to files

write2Files(subCmpId, cmp.second, gtf, fa);


d_stream_new2.clear(); //reuse the same stream (just clear the state flags)
stringstream outFile;
//outFile.str(""); //clear the stringstream
//outFile<<"./"<<dirPrefix<<"_"<<cmp.first<<"/obsRCcounts.txt"; //directory already exists from previous 
outFile<<"./"<<subCmpId<<"/obsRCcounts.txt";
string d_file_new2=outFile.str();

d_stream_new2.open(d_file_new2.c_str());
	if(!d_stream_new2){
		cout<<"Unable to open" <<d_file_new2<<endl;
		exit(1);
	}
	
	//for each read class ID
	for (const auto &class_ID : rcCountsLocal[cmp.first]){
		d_stream_new2<<"[\t";
		
		//for each transcript in this read class rcID_reverse[class_ID]-->print the class and the counts
		for (const auto &tr : rcID_reverse[class_ID])
		{
			d_stream_new2<<tr<<"\t";
		}
		d_stream_new2<<"]"<<"\t"<<obsReadsClasses[rcID_reverse[class_ID]]<<endl;
	}//end: for each read class ID	

d_stream_new2.close();
	

//Write to file the number of transcripts in each readClass
d_stream_new2.clear(); //reuse the same stream (just clear the state flags)
outFile.str(""); //clear the stringstream
//outFile<<"./"<<dirPrefix<<"_"<<cmp.first<<"/obsRCsize.txt"; //directory already exists from previous 
outFile<<"./"<<subCmpId<<"/obsRCsize.txt";
d_file_new2=outFile.str();

d_stream_new2.open(d_file_new2.c_str());
	if(!d_stream_new2){
		cout<<"Unable to open" <<d_file_new2<<endl;
		exit(1);
	}
	
	
	//Matrix for size of read class
	vector<vector<double>> rcSize;
	//initialize matrix rcSize with zero
	
	for(int i=0;i<rcCountsLocal[cmp.first].size();i++)
    {
		vector<double> row; // Create an empty row
         for(int j=0;j<rcCountsLocal[cmp.first].size();j++)
         {
              row.push_back(0.0);   // Add an element (column) to the row
         }
		rcSize.push_back(row); // Add the row to the main vector		 
	}	
	
	
	/*cout<<"print(rcSize)"<<endl;
	for(int i=0;i<rcSize.size();i++)
    {
         for(int j=0;j<rcSize.size();j++)
         {
              cout<<rcSize[i][j];   
         }
		cout<<endl;
	}*/	
	
	//fill the diagonal with read classes size.
	//for each read class ID
	int i=0;
	for (const auto &class_ID : rcCountsLocal[cmp.first])
	{
		rcSize[i][i]=(double)1/rcID_reverse[class_ID].size();
		i++;
	}
	
	//print Matrix to file
	for(int i=0;i<rcCountsLocal[cmp.first].size();i++)
    {
         for(int j=0;j<rcCountsLocal[cmp.first].size();j++)
         {
				d_stream_new2<<rcSize[i][j]<<"\t";
         }   
         d_stream_new2<<endl;
	}	

d_stream_new2.close();

}//end:for each cc	


resultStream.close();
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}//end splitCmp


void write2Files(int cmpID, vector<string> transcripts, map<string, map<string, int>> &gtf, map <string, string> &fa){

stringstream command1;
command1<<"mkdir "<<cmpID;
system(command1.str().c_str());

ofstream outStreamObsTr;
ofstream outStreamFa;
ofstream outStreamTrLength;

stringstream outObsTrFile;
stringstream outFaFile;
stringstream outTrLengthFile;

outObsTrFile<<"./"<<cmpID<<"/obsTr.txt";
outFaFile<<"./"<<cmpID<<"/"<<cmpID<<".fa";
outTrLengthFile<<"./"<<cmpID<<"/tr_length.txt";

string obsTrFile=outObsTrFile.str();
string faFile=outFaFile.str();
string trLengthFile=outTrLengthFile.str();

outStreamObsTr.open(obsTrFile.c_str());
outStreamFa.open(faFile.c_str());
outStreamTrLength.open(trLengthFile.c_str());

	if(!outStreamObsTr){
		cout<<"Unable to open " <<outObsTrFile<<endl;
		exit(1);
	}
	if(!outStreamFa){
		cout<<"Unable to open " <<faFile<<endl;
		exit(1);
	}
	if(!outStreamTrLength){
		cout<<"Unable to open " <<trLengthFile<<endl;
		exit(1);
	}

	//for each transcript -- print gene name -- transcript name -- transcript length
	for (const auto &tr : transcripts){
		
			outStreamObsTr<<gtf[tr].begin()->first<<"\t"<<tr<<endl;
			outStreamTrLength<<tr<<"\t"<<gtf[tr].begin()->second<<endl;
			outStreamFa<<">"<<tr<<endl<<fa[tr]<<endl;
	}//end: for each transcript	

outStreamObsTr.close();
outStreamFa.close();
outStreamTrLength.close();
}

/*
void write2FileObsTr(){
ofstream out;
stringstream outFile;
//outFile<<"./"<<dirPrefix<<"_"<<cmp.first<<"/obsTr.txt";
outFile<<"./"<<subCmpId<<"/obsTr.txt";
string d_file_new2=outFile.str();

}

void writeTrSeq2FaFile(int cmpID){
//Write to FASTA File the transcripts sequences

}


void write2FileTrLengths(){

}
*/
void print(const vector<pair<int,int> > &myVector){
	for (const auto& kv : myVector) {
		cout<<"\t\t\t\t"<<kv.first<<"	"<<kv.second<<endl;
	}
	cout<<endl;
}

void print(const map<string, int > & myMap) {

	for (const auto& kv : myMap) {
		std::cout << kv.first << " has value " << kv.second << std::endl;
	}
}

void print(const map<string, double > & myMap) {

	for (const auto& kv : myMap) {
		std::cout << kv.first << " = " << kv.second << std::endl;
	}
}


// ~~~~~ Print Transcripts Length ~~~~~
void print(const map<string, map<string, int > > &myMap) {
cout<<"~~~~~ Print Transcripts Length ~~~~~"<<endl;
	for (const auto& it: myMap)
	for (const auto& kv : it.second) {
		cout<<"["<<it.first<<"][";
		cout << kv.first<<"]"<< " has length " << kv.second <<endl;
	}
}

// ~~~~~ Print Normalized Adjusted Transcripts Length ~~~~~
void print(const map<string, map<string, double > > &myMap) {
cout<<"~~~~~ Print Normalized Adjusted Transcripts Length ~~~~~"<<endl;
	for (const auto& it: myMap)
	for (const auto& kv : it.second) {
		cout<<"["<<it.first<<"][";
		cout <<kv.first<<"]"<<" has length " << kv.second <<endl;
	}
}

void print(const map<string, vector<pair<int, int> > > &myMap) {
	for (const auto &kv : myMap) {
		cout <<"\t\t\tLocal Exons for transcript "<<kv.first << " :"<<endl; 
		vector<pair<int, int> > vect=kv.second;
			for (int it=0; it<vect.size(); it++){
				cout<<"\t\t\t\t"<< vect[it].first<<"\t"<<vect[it].second<<endl;
			}
	}

}

void print(const map<string, map<string, vector<pair<int, int> > > > &myMap) {
cout<<"\nPRINT ALL EXONS"<<endl<<endl;
	for (const auto &it : myMap){
		cout<<"Gene: "<<it.first<<" ("<<it.second.size()<<" transcript(s))"<<endl;
		for (const auto &kv : it.second) {
		cout <<"\tTranscript: "<<kv.first<<endl; 
		vector<pair<int, int> > vect=kv.second;
			for (int it=0; it<vect.size(); it++){
				cout<<"\t\t"<< vect[it].first<<"\t"<<vect[it].second<<endl;
			}
		}
	}

}

void print(const map<string, vector<int> > &myMap) {

	for (const auto &kv : myMap) {
		cout <<"Positions for transcript "<<kv.first << " :"<<endl; 
		vector<int> vect=kv.second;
			for (int it=0; it<vect.size(); it++){
				cout<<" "<< vect[it];
			}
		cout<<endl;
	}

}

void print(const map<string, map<string, vector<int> > > &myMap) {

for (const auto &it: myMap){
	for (const auto &kv : it.second) {
		cout <<"Positions for transcript "<<kv.first << " :"<<endl; 
		vector<int> vect=kv.second;
			for (int it=0; it<vect.size(); it++){
				cout<<" "<< vect[it];
			}
		cout<<endl;
	}

}
}


void print(const map<vector<int>, int > &myMap) {
cout<<"~~~~~ Print read counts ~~~~~"<<endl;
cout<<"Size of each class \t Read Classes"<<endl;
	for (const auto &kv : myMap) {
		cout <<kv.second<<"\t\t";
			for(int i=0; i<kv.first.size(); i++)
				cout<<"\t"<<kv.first[i];
			cout<<endl;
			
	}

}

void print(const map< map<string, bool> , double> &myMap) {
cout<<"~~~~~ Print observed/estimated reads frequency ~~~~~"<<endl;
cout<<"Read Class \t Size"<<endl;

	for (const auto &kv : myMap) {
		for (const auto &cluster : kv.first){
				cout<<cluster.second<<" ";
		}
		
		cout <<"\t\t"<< kv.second<<endl;		
	}

}



void print(const map<string, map< map<string, bool> , int> > &myMap) {
cout<<"~~~~~ Print local virtual reads ~~~~~"<<endl;
cout<<"Read Class \t Size"<<endl;

	for (const auto &tr : myMap){
	cout<<tr.first<<endl;
		for (const auto &kv : tr.second) {
			for (const auto &cluster : kv.first){
				cout<<cluster.second<<" ";
			}
		
			cout <<"\t\t"<< kv.second<<endl;		
		}
	cout<<endl;
	}

}

void print(const map<vector<int> , vector<pair<int, int> > > &myMap) {
cout<<"~~~~~ Print origins ~~~~~"<<endl;
cout<<"Read Cluster \t\t Origins"<<endl;

	for (const auto &kv : myMap) {
		for (const auto &cluster : kv.first){
			cout << cluster << " ";
		}
		cout<<endl;
		
		for (const auto &segment : kv.second){
			cout <<"\t\t"<< segment.first<<"\t"<<segment.second<<endl;
		}
		cout<<endl;			
	}

}

void loadGTF(string &gtf_file, map<string, map<string, int>> &gtf){

	ifstream gtfIfStream;
	gtfIfStream.open(gtf_file.c_str());
	if(!gtfIfStream){
		cout<<"Unable to open "<<gtf_file.c_str()<<endl;
		exit(1);
	}
		
	string line;
	string field;
	string type;
	int startExon=0;
	int endExon=0;
	int trLength=0;
	string geneName;
	string trName;
	
	while(gtfIfStream.good()){
		getline(gtfIfStream,line);	
		
		if((!line.empty())){

			istringstream iss(line);
			
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			type=field.c_str();
			
			if(type=="exon"){
								
				getline(iss,field,'\t');
				startExon=atoi(field.c_str());
				getline(iss,field,'\t');
				endExon=atoi(field.c_str());
				getline(iss,field,'\t');
				getline(iss,field,'\t');
				getline(iss,field,'\t');
				getline(iss,field,' '); //there is a space between gene_id and gene name
				geneName=field.c_str();
				
				if(geneName!="gene_id")
					getline(iss,field,' '); //skip space - Ensembl gtf files contains a space before gene_id
				
							
				getline(iss,field,' '); //there is a space between gene name and transcript_id
				geneName=field.c_str();
				geneName.erase(0,1); //remove the first quote mark
				geneName.erase(geneName.end()-2,geneName.end()); //delete the last quote mark and semicolon
				
				getline(iss,field,' '); //there is a space between gene name and transcript_id
				getline(iss,field,' '); //get transcript name
				trName=field.c_str();
				trName.erase(0,1); //remove the first quote mark
				trName.erase(trName.end()-2,trName.end()); //delete the last quote mark and semicolon
				
				//this eliminates the need of sorting gtf
				gtf[trName][geneName]+=endExon-startExon+1;
				//cout<<trName<<" "<<geneName<<" "<<gtf[trName][geneName];
			}
		}
	}
}

void loadFA(string &fa_file, map<string, string> &fa){
	
	ifstream inIfStream;
	inIfStream.open(fa_file.c_str());
	if(!inIfStream){
		cout<<"Unable to open "<<fa_file.c_str()<<endl;
		exit(1);
	}
		
	string line;
	string field;
	string trName;
	int tempCount=0;
	
	while(inIfStream.good()){
		getline(inIfStream,line);	
		
		if((!line.empty())){

			//cout<<"line = "<<line<<endl;
			if(line.substr(0,1)==">"){
				istringstream iss(line);
				getline(iss,field,'|'); //some fasta files have a | after transcript name
				trName=field.c_str();
				trName.erase(0,1);		
			}
			else{
				fa[trName]+=line;
			}
			
			
			/*tempCount++;
			if(tempCount==3){
				cout<<trName<<" "<<fa[trName]<<endl;
				exit(7);
			}*/
			
		}
	}
}