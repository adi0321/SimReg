//============================================================================
// Project     : Life Technologies Collaborative Research Grant
// Author      : Serghei Mangul / Adrian Caciula
// Description : calculate d for the bipartite graph(transcripts and reads)
// Created     : 12/29/2013
//
// Compile command:
// g++ -std=c++0x ./computeTrFreq.cpp ../include/current_time.cpp ../include/norm.cpp ../include/print.cpp ../include/prepro_extract_read_classes_cc.cpp -o ./computeTrFreq -I /usr/local/include/ -L /usr/local/lib/
//============================================================================

//This is an extension of version of ~/code/MCReg/annotation_preprocessing/compute_read_classes_v2.cpp (and compute_RC_CC_v2.cpp)
//This version computes Connected Components in addition to MC read classes and d values
	
#include "../include/current_time.h"
#include "../include/norm.h"
#include "../include/print.h"
#include "../include/prepro_extract_read_classes_cc.h"

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

#define HELPMESSAGE "Too few parameters: \n\
Input:\n\
argv[1] - gtf file \n\
argv[2] - File with local results \n\
argv[3] - total_sum_reads_portion \n\
"

double compare(vector<double> a, vector<double> b);
void getTrNames(string gtf,	map<string, double> &tr_freq);
			
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

double A(double a[], double b[], int NUM_TRANSCRIPTS);

int main(int argc,char *argv[]){

cout<<"\nRunning "<<argv[0]<<endl;

	if(argc<4){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

double EPS = 0.01;
int iterations=0;

	string gtf=argv[1]; //gtf
	string resultsFile=argv[2]; //raw Monte Carlo reads (grinder1-reads.fa for Grinder or true.sam for Marius simulator, respectively)
	double total_sum_reads_portion=atof(argv[3]);
	
	map<string, double> tr_freq;
	
	cout<<"\nRun "<<argv[0]<<endl;
	cout<<"\n["<<current_time()<<"] Collecting Transcript Names ... "<<endl;

		
		//~~~~~ get transcript names and load them in tr_freq
		getTrNames(gtf, tr_freq); 
		//see function definition at the end

		cout<<"\n\tInput Values:"<<endl;
		cout<<"\t\tTotal number of Transcripts: "<<tr_freq.size()<<endl;
		//print(tr_freq);
		
	cout<<"\n["<<current_time()<<"] Done!\n"<<endl;


//Compute MCReg Frequencies
//Open the file and compute /load the frequencies into tr_freq structure
//frequency=$(echo "${solver_estimates[$i]} * ( ${no_reads_c} / $local_avg_tr_len ) / $total_sum_reads_portion" | bc -l)


	string line;
	string field;
	
	ifstream inResults;
	
	inResults.open(resultsFile.c_str());
	if(!inResults){
		cout<<"Unable to open "<<resultsFile.c_str()<<endl;
		exit(1);
	}
	
	getline(inResults,line);//skip the headers line
	
	while(inResults.good())
	{
		getline(inResults,line);
		
		if(!line.empty())
		{
			istringstream iss(line);
			getline(iss,field,'\t'); //component ID
			getline(iss,field,'\t'); //number of transcripts
			getline(iss,field,'\t'); //no_reads_c , number of reads in current component
				//cout<<"Field = "<<field<<endl;
				int no_reads_c=atoi(field.c_str());
			getline(iss,field,'\t'); //transcript name
				string tr_name=field.c_str();
			getline(iss,field,'\t'); //local transcript frequency from solver
				double solver_estimates=atof(field.c_str());
			getline(iss,field,'\t'); //average transcript length in current component
				double local_avg_tr_len=atof(field.c_str());
			
			double frequency = solver_estimates * (no_reads_c / local_avg_tr_len) / total_sum_reads_portion;
			//$(echo "${solver_estimates[$i]} * ( ${no_reads_c} / $local_avg_tr_len ) / $total_sum_reads_portion" | bc -l)
				
				#if DEBUG
					cout<<"no_reads_c = "<<no_reads_c<<endl;
					cout<<"tr_name = "<<tr_name<<endl;
					cout<<"solver_estimates = "<<solver_estimates<<endl;
					cout<<"local_avg_tr_len = "<<local_avg_tr_len<<endl;
					cout<<"total_sum_reads_portion = "<<total_sum_reads_portion<<endl;
					cout<<"Transcript Frequency = "<<frequency<<endl;
				#endif
				
			tr_freq[tr_name]=frequency;
			
		}
	}
	
	inResults.close();
	
	
	
	//cout<<"Print transcript Frequencies"<<endl;
	//print(tr_freq);
	
ofstream resultStream;
string finalResultsFile="../MCReg.iso.estimates";
resultStream.open(finalResultsFile.c_str(),ios::app);//open file in append mode (in order to avoid overwritting)
	
if(!resultStream){
	cout<<"Unable to open " <<finalResultsFile<<endl;
		exit(1);
}

for (const auto &tr : tr_freq)
{
	resultStream<<tr.first<<"\t"<<tr.second<<endl;
}

resultStream.close();

	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

//----------------------------------------------------------------------------------------------------------
//---- Get ALL Transcripts Names ----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
//open gtf and get all information
void getTrNames(string gtf,	map<string, double> &tr_freq){
	//Warning: if input gtf is not sorted by the gene name then the same gene may appear several times
	
	string line;
	string field;
	string gene_name;
	string isoform;
	string isoform_prev;
	ifstream ingtf;
	int flag_first=0;
	int x_prev=0;
	int y_prev=0;

	ingtf.open(gtf.c_str());
	if(!ingtf){
		cout<<"Unable to open "<<gtf.c_str()<<endl;
		exit(1);
	}

	int current_length=0;
	
	while(ingtf.good()){
		getline(ingtf,line);
		
			//skip the headers
	 		while(line.substr(0,1)=="#")
	 			getline(ingtf,line);
				
		if(!line.empty()){
			istringstream iss(line);
			//cout<<line<<endl;
			getline(iss,field,'\t');
			//cout<<field<<endl;
			getline(iss,field,'\t');
			getline(iss,field,'\t');
				string lineType=field.c_str();
				if(lineType.compare("gene")==0)
					continue;
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'"');
			getline(iss,field,'"');
			getline(iss,field,'"');
			getline(iss,field,'"');
				isoform=field.c_str();

			tr_freq[isoform]=0.0;
		}
	}
	
ingtf.close();
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