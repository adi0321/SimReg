//============================================================================
// Project     : SimReg
// Author      : Adrian Caciula
// Description : calculate d for the bipartite graph(transcripts and reads)
// Created     : 10/20/2014
// Last Update: 
// Compile command:
// g++ -std=c++0x ./combineFiles.cpp ../include/current_time.cpp ../include/norm.cpp ../include/print.cpp -o ./combineFiles -I /usr/local/include/ -L /usr/local/lib/
//============================================================================


#include "../include/current_time.h"
#include "../include/norm.h"
#include "../include/print.h"

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
#include <iostream>
#include <utility>


//DEBUGGING CONSTANTS:
#define DEBUG 0

using namespace std;

#define HELPMESSAGE "too few parameters: \n\
Input:\n\
argv[1] - singleTrGenes.txt \n\
argv[2] - trLen.txt \n\
argv[3] - results.txt \n\
"
		
void loadTrLen(string &inFile, map<string, int> &tr);

int main(int argc,char *argv[]){

//cout<<"\nRunning "<<argv[0]<<endl;

	if(argc<4){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}


	string singleTrGenes_file=argv[1]; //aligned Observed reads
	string trLen_file=argv[2];
	string resultsFile=argv[3];
	
clock_t begin = clock(); //used for measuring entire elapsed time for this function
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Load Transcripts lengths"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		map<string, int> tr;
		loadTrLen(trLen_file, tr);
				
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;	

//Append singleTrGenes.txt and trLen.txt to results
ofstream resultStream;
ifstream inSingleTrGenesStream;

inSingleTrGenesStream.open(singleTrGenes_file.c_str());
resultStream.open(resultsFile.c_str(),ios::app);//open file in append mode (in order to avoid overwritting)

	if(!resultStream){
		cout<<"Unable to open " <<resultsFile<<endl;
		exit(1);
	}
	
	if(!inSingleTrGenesStream){
		cout<<"Unable to open " <<singleTrGenes_file<<endl;
		exit(1);
	}

	string line;
	string field;
	string trName;
	
	while(inSingleTrGenesStream.good()){
		getline(inSingleTrGenesStream,line);	
		
		if((!line.empty())){
			
			istringstream iss(line);
						
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'\t');
	
			getline(iss,field,'\t');
			trName=field.c_str();
	
			resultStream<<line<<"\t"<<tr[trName]<<endl;
		}
	}

resultStream.close();
inSingleTrGenesStream.close();
		
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


void loadTrLen(string &inFile, map<string, int> &tr){

	ifstream inIfStream;
	inIfStream.open(inFile.c_str());
	if(!inIfStream){
		cout<<"Unable to open "<<inFile.c_str()<<endl;
		exit(1);
	}
		
	string line;
	string field;
	string type;
	int trLength=0;
	string trName;
	
	while(inIfStream.good()){
		getline(inIfStream,line);	
		
		if((!line.empty())){

			istringstream iss(line);
			
			getline(iss,field,'\t');
			trName=field.c_str();
			
			getline(iss,field,'\t');
			trLength=atoi(field.c_str());
			tr[trName]=trLength;
		}
	}
}
