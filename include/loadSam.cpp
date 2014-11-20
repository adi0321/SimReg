#include "loadSam.h"
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <sstream>      // std::istringstream
#include <string>
#include <map>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#define DEBUG 0

/********************************************************************
 * Load SAM File *
 ********************************************************************/
 
//For MAQC Read Names we need to extract the second half from the name of the reads
//Here is an example of one read mapped to multiple transcripts
//SRR039628.8201139       99      ENST00000461638.2       621     255     50M     =       798     227     GCCCACCCCCCTGGACCCTGAAGAGACAGCCTACCCTAGCCTGAGTGGGG      <<<<<<<<<<<9<<<<<<39.9<<<<<:<<<:;<<<662<:888884448     XA:i:0  MD:Z:50 NM:i:0
//SRR039629.8201139       147     ENST00000461638.2       798     255     50M     =       621     -227    CCTAAGCAATCCCAACCTCCAGGCTTCCCTGAGCAGTCCTCAGCCCCAGC      *79<<9;;<;477<;;:::;<<<<5<:<<<<<<<<<<<<<<<<<<<<<<<     XA:i:0  MD:Z:50 NM:i:0
//SRR039628.8201139       99      ENST00000368633.1       989     255     50M     =       1166    227     GCCCACCCCCCTGGACCCTGAAGAGACAGCCTACCCTAGCCTGAGTGGGG      <<<<<<<<<<<9<<<<<<39.9<<<<<:<<<:;<<<662<:888884448     XA:i:0  MD:Z:50 NM:i:0
//SRR039629.8201139       147     ENST00000368633.1       1166    255     50M     =       989     -227    CCTAAGCAATCCCAACCTCCAGGCTTCCCTGAGCAGTCCTCAGCCCCAGC      *79<<9;;<;477<;;:::;<<<<5<:<<<<<<<<<<<<<<<<<<<<<<<     XA:i:0  MD:Z:50 NM:i:0
//SRR039628.8201139       99      ENST00000303569.6       629     255     50M     =       806     227     GCCCACCCCCCTGGACCCTGAAGAGACAGCCTACCCTAGCCTGAGTGGGG      <<<<<<<<<<<9<<<<<<39.9<<<<<:<<<:;<<<662<:888884448     XA:i:0  MD:Z:50 NM:i:0
//SRR039629.8201139       147     ENST00000303569.6       806     255     50M     =       629     -227    CCTAAGCAATCCCAACCTCCAGGCTTCCCTGAGCAGTCCTCAGCCCCAGC      *79<<9;;<;477<;;:::;<<<<5<:<<<<<<<<<<<<<<<<<<<<<<<     XA:i:0  MD:Z:50 NM:i:0

 
void loadSam(string &inFile, map<string, map<string, map<int, string>>> &sam, bool flagMAQC){
	
	/* initialize random seed: */
	//srand (time(NULL));
	
	ifstream inIfStream;
	inIfStream.open(inFile.c_str());
	if(!inIfStream){
		cout<<"Unable to open "<<inFile.c_str()<<endl;
		exit(1);
	}
	
	
	ofstream outSam;
	string resultsFile="subSample.sam";
	outSam.open(resultsFile.c_str(),ios::app);//open file in append mode (in order to avoid overwritting)

	if(!outSam){
		cout<<"Unable to open " <<resultsFile<<endl;
		exit(1);
	}
	
		
	string line;
	string field;
	string readName;
	int readFlag;
	string trName;
	string readToSkip="";
	bool firstRead=true;
	
	while(inIfStream.good()){
		getline(inIfStream,line);	
		
		//skip the headers
	 	while(line.substr(0,1)=="@"){
			outSam<<line<<endl;
			getline(inIfStream,line);//get next line
		}
		
		if((!line.empty())){

			istringstream iss(line);
			
			//Get read name
			getline(iss,field,'\t');
				readName=field.c_str();
			
			//For MAQC we only need to keep the name after the extension
			if(flagMAQC) 
				readName = readName.substr(readName.find_first_of('.')+1);
			
			/*
			if(readName.compare(readToSkip) == 0)
				continue; 
				
			//Else Flip a coin
			int flipResult;
			flipResult = rand() % 2 + 1;
					
			if(flipResult == 1){
				readToSkip=readName;
				continue;
			}
			*/
			//Else keep the read
			getline(iss,field,'\t');
				readFlag=atoi(field.c_str());
			getline(iss,field,'\t');
				trName=field.c_str();

			if(trName.compare("*")==0)
				continue;
			
			//This ovoids overwriting since we alway add in addition that what we already have
			sam[readName][trName].insert(pair<int,string>(readFlag, line));
			
			#if DEBUG		
				//Redirect the read to another file
				//resultStream<<line<<endl;
			
			
				cout<<"readName = "<<readName<<endl;
				cout<<"readFlag = "<<readFlag<<endl;
				cout<<"trName = "<<trName<<endl;
				cout<<"Entire line: "<<line<<endl;
				
				//exit(7);
			#endif
		}
	}
	
		#if DEBUG
				for(const auto &read : sam){
					//cout<<read.first<<endl;
					for(const auto &t : read.second){
						//cout<<t.first<<endl;
						for(const auto &f : t.second){
							//cout<<f.first<<endl;
							cout<<f.second<<endl;
						}
					}
				}
		#endif
	
	inIfStream.close();
	outSam.close();

}//end loadSam


