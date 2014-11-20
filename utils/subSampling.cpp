//============================================================================
// Project     	: SimReg Simulated Regression Algorithm for Transcriptome Quantification from RNA-Seq Data
// Author      	: Adrian Caciula, Department of Computer Science, Georgia State University, Atlanta, GA 30303, USA 
// Contact		: Adrian Caciula at adrian.caciula@gmail.com   
// Description 	: SubSampling the SAM File
// 					- Input: SAM File:
//					- Output: New SamFile 

// Created     	: 23/10/2014
// @Copyright 2014, All Rights Reserved.

// Compile command:
// g++ -std=c++0x subSampling.cpp ../include/current_time.cpp ../include/loadSam.cpp -o subSampling 
//============================================================================

#include "../include/current_time.h"
#include "../include/loadSam.h"

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <sstream>      // std::istringstream
#include <string>
#include <vector>
#include <algorithm>    // std::reverse, copy
#include <iterator>     // std::ostream_iterator
#include <ctime>
#include <map>

#define DEBUG 0

using namespace std;

#define HELPMESSAGE "Usage: \n\
argv[1] - SAM File with Transcriptome alignments. \n\
		- Align reads with the folloing command: bowtie --best -v 3 -k 60 -p 16 --chunkmbs 128 $REF_File ${input_option} -1 $obs_pair1_file -2 $obs_pair2_file -I $fil_min -X $fil_max -S $PWD/bowtie_OBS_k60v3best.sam \n\
argv[2] - true or false Specify if the given SAM file represents MAQC mapped reads \n\
"

int main(int argc,char *argv[]){
	
	cout << endl <<" SubSampling SAM File"<<endl; 
    cout << "==================================================="<<endl<<endl;
	
	if(argc<2){
		cout<<"Error: Too few parameters!"<<endl;
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

string samFile=argv[1];
bool flagMAQC=argv[2];

cout<<samFile<<endl;
cout<<"subSampling: flagMAQC = "<<flagMAQC<<endl;

	cout<<"\n["<<current_time()<<"] Run: "<<argv[0]<<endl;
	
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Load SAM File"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

	//map<readName, map<tr_name, pair<flag, wholeLine>>>
	map<string, map<string, map<int, string>>> sam;
	//load and redirect
	loadSam(samFile, sam, flagMAQC);
				
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;	


	ofstream outSam;
	string resultsFile="subSample.sam";
	outSam.open(resultsFile.c_str(),ios::app);//open file in append mode (in order to avoid overwritting the headers added by loadSam)

	if(!outSam){
		cout<<"Unable to open " <<resultsFile<<endl;
		exit(1);
	}

//1. Subsampling Simulated Data //No need anymore
//For Each read generate, flip the coin, and keep only heads
//Redirect the read to another sam file


//2. Subsampling Real Data
//add index to sam structure
//map<index, map<readName, map<tr_name, pair<flag, wholeLine>>>>
map<int, map<string, map<string, map<int, string>>>> samWithIndex;
int index=1;
for(const auto &read : sam){
	samWithIndex[index][read.first]=read.second;
	index++;
}

			#if DEBUG
			for(const auto &index : samWithIndex){
				cout<<"Index = "<<index.first<<endl;
				for(const auto &read : index.second){
					//cout<<read.first<<endl;
					for(const auto &t : read.second){
						//cout<<t.first<<endl;
						for(const auto &f : t.second){
							//cout<<f.first<<endl;
							cout<<f.second<<endl;
						}
					}
				}
			}
			#endif

//N = total number of reads 
const int N = samWithIndex.size();
int readID=1;

//unique reads
cout<<"Total number of reads: "<<N<<endl;
cout<<"Subsampling! Please wait..."<<endl;

//Randomly pick a read from the original SAM File and redirect it to a new SAM file (assign a new id (readID)
//then readID++;

/* initialize random seed: */
srand (time(NULL));
//loop from 1 to N
for(int i=0; i<N; i++){
	int randIndex=rand() % N + 1;
	//cout<<"Random Index = "<<randIndex<<endl;
	
	//send read with index randIndex to file
				for(const auto &read : samWithIndex[randIndex]){
					//cout<<read.first<<endl;
					for(const auto &t : read.second){
						//cout<<"t.first"<<endl<<t.first<<endl;
						for(const auto &f : t.second){
							//cout<<f.first<<endl;
							//cout<<f.second<<endl;
							string line = f.second;
							string part = line.substr(line.find_first_of('\t')+1);
							//1 was added to exclude the delimiter
							
							//cout<<"Part"<<endl<<part<<endl;
							//cout<<readID<<"\t"<<part<<endl;
							outSam<<readID<<"\t"<<part<<endl;
						}
					}
					readID++;
				}
}

cout<<"Done SubSampling"<<endl;
outSam.close();
}