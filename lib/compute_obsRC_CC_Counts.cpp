//============================================================================
// Project     : Life Technologies Collaborative Research Grant
// Author      : Serghei Mangul / Adrian Caciula
// Description : calculate d for the bipartite graph(transcripts and reads)
// Created     : 06/01/2013
// Last Update: 12/16/2013
// Compile command:
// g++ -std=c++0x ./compute_obsRC_CC_Counts.cpp ../include/current_time.cpp ../include/norm.cpp ../include/normD.cpp ../include/print.cpp ../include/extract_obsRC.cpp ../include/prepro_extract_read_classes_cc.cpp -o ./compute_obsRC_CC_Counts -I /usr/local/include/ -L /usr/local/lib/
//============================================================================

//!!!!WARNING: This program is only when SimReg reads are used as observed reads (perfect case)

//This is an extension of version of ~/code/MCReg/annotation_preprocessing/compute_read_classes_v2.cpp (and compute_RC_CC_v2.cpp)
//This version computes Connected Components in addition to MC read classes and d values
//2013.12.16 - Compute Observed RC and CC
	
#include "../include/current_time.h"
#include "../include/norm.h"
#include "../include/normD.h"
#include "../include/print.h"
#include "../include/extract_obsRC.h"
#include "../include/prepro_extract_read_classes_cc.h"

#include <iostream>
#include <stdlib.h> // atof
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
argv[1] - File with mapped Observed reads (from Bowtie) \n\
argv[2] - File with True Transcript Frequencies \n\
argv[3] - Reads File to extract frequency (for SimReg reads either pair 1 or 2) \n\
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



int main(int argc,char *argv[]){

cout<<"\nRunning "<<argv[0]<<endl;

	if(argc<5){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

double EPS = 0.01;
int iterations=0;

	string readsBowtie_file=argv[1]; //aligned Observed reads
	string trueTrFreq_file=argv[2]; // file with true transcripts frequencies
	string mcReadsFile=argv[3]; //raw reads to extract the reference (SimReg reads either pair1 or pair2)
	int totalNReads=atoi(argv[4]); //Total number of generated reads
	
clock_t begin = clock(); //used for measuring entire elapsed time for this function
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Compute Observed Reads Classes "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		
//Extract OBS reads classes from Bowtie output
map<vector<string>, int> obsReadsClasses; 

cout<<"\nMain: Parsing Bowtie file ..."<<endl;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//!!!!!We also need the matrix with the references
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cout<<"\n["<<current_time()<<"] Extracting reads References..\n"<<endl;
	map<string, string> mcReadsRef;

	//This function is located in prepro_extract_read_classes_cc.cpp
	extract_reads_references(mcReadsFile, mcReadsRef);

cout<<"\n\t["<<current_time()<<"] Reads References had been collected"<<endl;
		
		#if DEBUG_READ
			cout<<"Print MC Reads references"<<endl;
			print(mcReadsRef);
			
			//exit(7);
		#endif

map<vector<string>, map<string, int> > mcReadsClassesRef; //references for OBS reads classes.
//references for reads classes. Reports how many reads comes from each transcript.

prepro_extract_read_classes_cc(readsBowtie_file, obsReadsClasses, mcReadsRef, mcReadsClassesRef);
//extract_obsRC(obsRBowtie_file, obsReadsClasses);



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
			print(obsReadsClasses);
			cout<<endl;
	#endif


	
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;	

//cout<<"Adjust the counts based on the file with true frequencies"<<endl;
ifstream trFreqStream;
trFreqStream.open(trueTrFreq_file.c_str());
	if(!trFreqStream){
		cout<<"Unable to open "<<trueTrFreq_file.c_str()<<endl;
		exit(1);
	}
	
string line;
string field;
string trName;
double trFreq;
map<string, double> mapTrFreq;
while(trFreqStream.good()){
	getline(trFreqStream,line);	
		
	if(!line.empty()){

		istringstream iss(line);
			
		getline(iss,field,'\t');
			trName=field.c_str();
		
		getline(iss,field,'\t');
			trFreq=atof(field.c_str());
		
		mapTrFreq[trName]=trFreq;
		
			#if DEBUG
				cout<<"Transcript name = "<<trName<<endl;
				cout<<"Transcript frequency = "<<trFreq<<endl;
					
				//exit(7);
			#endif
	}
}

trFreqStream.close();

//cout<<"Print transcript frequencies:"<<endl;
//print(mapTrFreq);


//Compute how many transcripts where maps on each transcript
map< string, map< vector<string>, double> > d_value; //this is d_value

		//for each read class
		for (const auto &read_class: obsReadsClasses){
			
			//for each transcript in the current class
			for (const auto &tr: read_class.first){
				
				//if this value exists --> then we have an error
				if(d_value[tr][read_class.first]){
					cout<<"Error: d value already exists!!!"<<endl;
					exit(1);
				}
				
				//As for the first step we only insert the references
				d_value[tr][read_class.first]=mcReadsClassesRef[read_class.first][tr];
				
				#if DEBUG
				if(mcReadsClassesRef[read_class.first][tr]==0)
				{
					cout<<"Reference = "<<mcReadsClassesRef[read_class.first][tr]<<endl;
					cout<<"Transcript "<<tr<<" does not contribute to class "<<endl;
					print(read_class.first);
					
					//check the other transcripts from this class
					cout<<endl;
					for (const auto &temp_tr : read_class.first){
						cout<<"Transcript "<<temp_tr<<" contributes with "<<mcReadsClassesRef[read_class.first][temp_tr]<<" to class: "<<endl;
						print(read_class.first);
					}
				}
				
				//cout<<"mcReadsClassesRef[read_class.first]["<<tr.first<<"] = "<<mcReadsClassesRef[read_class.first][tr.first]<<endl;
				#endif
			}
		}
		
				#if DEBUG
					cout<<"\nPrint d_value BEFORE Normalization"<<endl;
					print(d_value);
					
					//exit(7);
				#endif


cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done: d values had been computed"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;

//Sum of all counts per transcript
map<string, int> sumTrCounts;

//Normalize d per comlumn
normD(d_value, sumTrCounts);

		#if DEBUG
			cout<<"\nPrint d_value AFTER Normalization"<<endl;
			print(d_value);
					
			cout<<"Total sum counts for each transcript:"<<endl;
			print(sumTrCounts);
			//exit(7);
		#endif

//Compute Crude Frequency:
map<string, double> crudeFreq;
double sumCrudeFreq=0.0;
for(const auto &tr : mapTrFreq)
{
	crudeFreq[tr.first]=tr.second*sumTrCounts[tr.first];
	sumCrudeFreq+=crudeFreq[tr.first];
}

//cout<<"Print Crude Freq step1:"<<endl;
//print(crudeFreq);

//Update Crude Freq
for(const auto &tr : crudeFreq)
	crudeFreq[tr.first]=crudeFreq[tr.first]/sumCrudeFreq;

	/*
		cout<<"Print Crude Freq step2:"<<endl;
		print(crudeFreq);

		cout<<"Current Observed Counts:"<<endl;
		print(obsReadsClasses);
	*/
	
//Update observed counts:
for(const auto &rc : obsReadsClasses)
{
	double newObsCounts=0.0;
	//loop through all transcripts involved in this class
	for(const auto &trName : rc.first)
	{	
		//new counts = TrCounts * crudeFreq
		newObsCounts+=(double)d_value[trName][rc.first]*crudeFreq[trName];
	}
	
	//cout<<"new Obs Counts = "<<newObsCounts<<endl;
	obsReadsClasses[rc.first]=newObsCounts*totalNReads;
}


//cout<<"Observed Counts after UPDATE:"<<endl;
//print(obsReadsClasses);
//cout<<"where N = "<<totalNReads<<endl;
//exit(7);


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
					print(p_value_new2);
					
					//exit(7);
				#endif


cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done: adjacency matrix has been computed"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;				

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/		

cout<<"\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Computing Connected Components ... "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

//Both DFS and BFS complete on O(n) time, using O(n) memory, where n is matrix size. But BFS it doesn't suffer from stack overflow problem, and it doesn't spend time on recursive calls.
//Anyway, let's first use boost library


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
		}
	}
	
	
    std::vector<int> component(num_vertices(G));
	
    int num = connected_components(G, &component[0]);
    cout << "Total number of components: " << num << endl;
	
	/*
	std::vector<int>::size_type i;
	for (i = 0; i != component.size(); ++i)
      cout << "Vertex " << i <<" is in component " << component[i]<<endl;
	  
    cout << endl;
	
	cout<<"\nPrint components: "<<endl;
	for (const auto &c : component) {
		cout<<"Component: "<<c<<endl;

	}
	
	MyGraph::vertex_iterator vertexIt, vertexEnd;
	boost::tie(vertexIt, vertexEnd) = vertices(G);
	for (; vertexIt != vertexEnd; ++vertexIt){
		VertexID vertexID = *vertexIt; // dereference vertexIt, get the ID
		Vertex & vertex = G.graph()[vertexID];
		//The side effect of boost::labeled_graph is that one need to use graph() member function to make some algorithms work
	
		cout<<"Vertex name is : "<<vertex.vertexName<<endl;
    }
	
	*/
	
	
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done: Connected Components had been computed"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
	
		
		
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~ Write Values to Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cout<<endl<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"Write Values to Files"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;

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
	
	//CREATE MAP THAT WILL HOLD ALL READ CLASSES AND COUNTS
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
				cout<<"component "<<cmp.first<<"containts: "<<endl;
					for(const auto &elem : cmp.second)
						cout<<"\t\t"<<elem<<endl;
			}
		#endif
	
//Now print to file
cout<<"\nWrite observed transcripts names\n\t"<<endl;

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


//Check the size of the component and if it is 1 then just write the results to file
if( cmp.second.size()< 2 )
{
	resultStream<<cmp.first<<"\t"<<cmp.second.size()<<"\t"<<obsReadsClasses[cmp.second];
	resultStream<<"\t"<<cmp.second.front()<<"\t1"<<endl;
	//cmpID---#tr.---#ObsReads---Transcript names

	continue; //continue to the next component
}

stringstream command1;
command1<<"mkdir "<<cmp.first;
system(command1.str().c_str());

d_stream_new2.clear(); //reuse the same stream (just clear the state flags)
stringstream outFile;
outFile<<"./"<<cmp.first<<"/obsTr.txt";
string d_file_new2=outFile.str();

d_stream_new2.open(d_file_new2.c_str());
	if(!d_stream_new2){
		cout<<"Unable to open " <<d_file_new2<<endl;
		exit(1);
	}
//for each transcript -- print transcript -- print class -- print d values
	for (const auto &tr : cmp.second){
		d_stream_new2<<tr<<endl;
	}//end: for each transcript	

d_stream_new2.close();
	
d_stream_new2.clear(); //reuse the same stream (just clear the state flags)
outFile.str(""); //clear the stringstream
outFile<<"./"<<cmp.first<<"/obsRCcounts.txt"; //directory already exists from previous 
d_file_new2=outFile.str();

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
