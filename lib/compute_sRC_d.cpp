//============================================================================
// Author      : Adrian Caciula
// Description : calculate d for the bipartite graph(transcripts and reads)
// Created     : 2013/06/01
// Last Update : 2014/05/29
// Compile command:
// g++ -std=c++0x ./compute_sRC_d.cpp ../include/current_time.cpp ../include/norm.cpp ../include/print.cpp ../include/prepro_extract_read_classes_cc.cpp ../include/extract_obsCounts.cpp -o ./compute_sRC_d -I /usr/local/include/ -L /usr/local/lib/
//============================================================================

//This is an extension of version of ~/code/MCReg/annotation_preprocessing/compute_read_classes_v2.cpp (and compute_RC_CC_v2.cpp)


#include "../include/current_time.h"
#include "../include/norm.h"
#include "../include/print.h"
#include "../include/prepro_extract_read_classes_cc.h"
#include "../include/extract_obsCounts.h"

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
#define DEBUG_G2J 0 //debug gtf2junction method
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
coverage \n\
Input:\n\
argv[1] - gtf file \n\
argv[2] - File with reads for getting the references: (e.g., for Grinder: grinder-reads.fa and for Marius simulator: true.sam) \n\
argv[3] - File with mapped Monte Carlo reads from Bowtie (default output file) \n\
argv[4] - File with observed reads classes counts \n\
"

double compare(vector<double> a, vector<double> b);
			
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
//void print(const map< map<string, bool> , int> &myMap); /delete this (already declared in print.h)
void print(const map< map<string, bool> , double> &myMap);
void print(const map<string, map< map<string, bool> , int> > &myMap);

int find_index( const vector<int> &where, int searchParameter );

double A(double a[], double b[], int NUM_TRANSCRIPTS);

int main(int argc,char *argv[]){
	if(argc<5){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

double EPS = 0.01;
int iterations=0;

	string gtf=argv[1]; //gtf
	string mcReadsFile=argv[2]; //reference file: raw Monte Carlo reads (grinder1-reads.fa for Grinder or true.sam for Marius simulator, respectively)
								//SimReg Simulator has same format as grinder
	string mcRBowtie_file=argv[3]; //aligned Monte Carlo reads
	string obsRCcounts_file=argv[4]; 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cout<<"\n["<<current_time()<<"] Extracting Monte Carlo reads References..\n"<<endl;
	map<string, string> mcReadsRef;

	//This function is located in prepro_extract_read_classes_cc.cpp
	extract_reads_references(mcReadsFile, mcReadsRef);
		
cout<<"\n\t["<<current_time()<<"] Monte Carlo Reads References had been collected"<<endl;
		
		#if DEBUG_READ
			cout<<"Print MC Reads references"<<endl;
			print(mcReadsRef);
			
			//exit(7);
		#endif


		
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Compute Monte Carlo Reads Classes "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		
//Extract MC reads classes from Bowtie output
map<vector<string>, int> mcReadsClasses; 

map<vector<string>, map<string, int> > mcReadsClassesRef; //references for MC reads classes.
//references for reads classes. Reports how many reads comes from each transcript.

cout<<"\nMain: Parsing Bowtie files ..."<<endl;

prepro_extract_read_classes_cc(mcRBowtie_file, mcReadsClasses, mcReadsRef, mcReadsClassesRef);

//Compute Total number of Monte Carlo reads
int total_mc_reads=0;
for (const auto &r_class : mcReadsClasses){
	total_mc_reads+=r_class.second;
}

cout<<"\n~~~ Total number of Monte Carlo reads = "<<total_mc_reads<<endl;
cout<<"Total number of Monte Carlo Reads classes is: "<<mcReadsClasses.size()<<endl;		

	#if DEBUG
		cout<<"\nMONTE CARLO Reads Classes \t Size:"<<endl;
			//print(mcReadsClasses);
			cout<<endl;
			
		cout<<"MC Reads Classes References"<<endl;
			//print(mcReadsClassesRef);
				//exit(7);
	#endif

cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;	
	
	
////2014/5/1 - We need to do the union of read classes
//For this we'll need to load observed read classes
map<vector<string>, int> obsReadsClasses;

int discarded_reads=extract_obsCounts(obsRCcounts_file, obsReadsClasses);

		#if DEBUG
			cout<<"Size Simulated Read Classes: "<<mcReadsClasses.size()<<endl;
			cout<<"Size Observed Read Classes: "<<obsReadsClasses.size()<<endl;
			//exit(7);
		#endif
//for each observed read class check if it is in monte carlo
for (const auto &obsRC : obsReadsClasses)
{
	//if observed read class cannot be found in MC read class than add it with a count of zero
	if ( mcReadsClasses.find(obsRC.first) == mcReadsClasses.end() )
	{
		mcReadsClasses[obsRC.first]=0; // if not found then add it to Monte Carlo classes
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


cout<<"\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Computing d values... "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

				
/////Compute new d values: d_value_new2
map< string, map< vector<string>, double> > p_value_new2; //this is d_value


		//for each MC class
		for (const auto &read_class: mcReadsClasses){
			
			//for each transcript in the current class
			for (const auto &tr: read_class.first){
				
				//if this value exists --> then we have an error
				if(p_value_new2[tr][read_class.first]){
					cout<<"Error: p value new2 already exists!!!"<<endl;
					exit(1);
				}
				
				p_value_new2[tr][read_class.first]=mcReadsClassesRef[read_class.first][tr];
				
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
					cout<<"\nPrint p_value_new2 BEFORE Normalization"<<endl;
					print(p_value_new2);
					
					//exit(7);
				#endif


cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done: d values had been computed"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;				

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/		
//No Need For Connected Components anymore

	
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~ Write Values to Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cout<<endl<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"Write Values to Files"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;

//Write MC Read Classes (and each transcript name) to File 
string read_classes_file=gtf.substr(0,gtf.find_last_of("."))+"_read_classes.txt";

cout<<"\nWrite read classes to: \n\t"<<read_classes_file<<endl;

ofstream read_classes_stream;

read_classes_stream.open(read_classes_file.c_str());
		
if(!read_classes_stream){
	cout<<"Unable to open" <<read_classes_file<<endl;
		exit(1);
}

//for each read class
for (const auto &read_class : mcReadsClasses) {
	
	//for each transcript in the read class
	for (const auto &tr : read_class.first){
		//print to file the class name and the transcript
		read_classes_stream<<"[ ";
		for (const auto &cluster : read_class.first){
				read_classes_stream<<cluster<<" ";
		}
		read_classes_stream<<" ] \t "<< tr<<endl;	
	}
}


read_classes_stream.close();
/*********************************
*********************************/

//Write MC Read Classes (only) to File 
string read_classes_file2=gtf.substr(0,gtf.find_last_of("."))+"_read_classes2.txt";

cout<<"\nWrite read classes to: \n\t"<<read_classes_file2<<endl;

read_classes_stream.clear(); //reuse the same stream (just clear the state flags - it should be ok)
read_classes_stream.open(read_classes_file2.c_str());
		
if(!read_classes_stream){
	cout<<"Unable to open" <<read_classes_file2<<endl;
		exit(1);
}

//for each read class
for (const auto &read_class : mcReadsClasses) {

		//print to file the class name
		read_classes_stream<<"[ ";
		for (const auto &cluster : read_class.first){
				read_classes_stream<<cluster<<" ";
		}
		read_classes_stream<<"]"<<endl;
}
read_classes_stream.close();
/*********************************************************************************
**********************************************************************************/

/****************************
*	Write d values to file	*
*****************************/
string d_file_new2=gtf.substr(0,gtf.find_last_of("."))+"_d_values.txt";

cout<<"\nWrite d_t,r values to: \n\t"<<d_file_new2<<endl;

ofstream d_stream_new2;

		d_stream_new2.open(d_file_new2.c_str());
		
	if(!d_stream_new2){
		cout<<"Unable to open" <<d_file_new2<<endl;
		exit(1);
	}
	

//Now print to file

//for each transcript -- print transcript -- print class -- print d values
	for (const auto &tr : p_value_new2){

		//for each class in the current transcript
		for (const auto &read_class : tr.second ){		
				d_stream_new2<<tr.first<<"\t[\t";

			//for each transcript in this class
			for (const auto &tr_class : read_class.first)
				d_stream_new2<<tr_class<<"\t";
				
				
			d_stream_new2<<"]\t"<<read_class.second<<"\n";
			//No need for scientific since now there are only integers
			//d_stream_new2<<"]\t"<<scientific<<read_class.second<<"\n";
			//d_stream_new2<<" ] \t "<<setprecision(50)<<read_class.second<<"\n";
			//d_stream_new2<<p_value_new2[tr.first][read_class.first]<<"\t";
			
		}//end: for each class
	}//end: for each transcript	


d_stream_new2.close();
		
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
			
cout<<"\nmain: Done"<<endl;
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


//----------------------------------------------------------------------------------------------------------
//calculate min(a_j-b_j)
//----------------------------------------------------------------------------------------------------------
double compare(vector<double> a, vector<double> b){
	cout<<"a[0]="<<a[0]<<" , b[0]="<<b[0]<<endl;
	double min=abs(a[0]-b[0]);
		for (int i=0; i<a.size(); i++){
			if (abs(a[i]-b[i]) < min)
				min=abs(a[i]-b[i]);
				
		}
			cout<<"min abs(a[0]-b[0])="<<min<<endl;
		return min;
	}


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

int find_index( const vector<int> &where, int searchParameter )
{
    for( int i = 0; i < where.size(); i++ ) {
       if(where[i]==searchParameter) {
           return i;
       }
    }
    return -1;
}



//----------------------------------------------------------------------------------------------------------
//calculate min(a_j-b_j)
//----------------------------------------------------------------------------------------------------------
double A(double a[], double b[], int NUM_TRANSCRIPTS){
	float min=abs(a[0]-b[0]);
		for (int i=0; i<NUM_TRANSCRIPTS; i++){
			if (abs(a[i]-b[i]) < min)
				min=abs(a[i]-b[i]);
				
		}			
		return min;
	}
//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------