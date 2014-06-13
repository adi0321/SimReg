#include "extract_obsRC.h"
#include "print.h"
#include <iostream>		// std::cout
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <algorithm>    // std::sort
#include <vector>       // std::vector

#define DEBUG 0

/********************************************************************************
*	Extract initial Observed Read Classes --- Required for computing the initial CC		*
*********************************************************************************/

using namespace std;

//Extract read classes from Bowtie output
//This version is for observed reads only (reduced version of extract_read_classes2.cpp)

void extract_obsRC(string &rBowtie_file, map<vector<string>, int> &ReadsClasses)
{
	//This version uses vector<string> in the map for ReadsClasses instead of map<string,bool>
	
	//cout<<"\n\textract_read_classes: "<<endl;
	//cout<<"\t\tParameters: "<<endl;
	//cout<<"\t\t\t	Bowtie file: "<<rBowtie_file<<endl;
	cout<<rBowtie_file<<endl;
	//cout<<"\t\t\t2. Reads Classes map "<<endl;
	cout<<"Please wait, this may take few minutes"<<endl;
	cout<<"Processing . . ."<<endl;
	/****************************************************
	*					Error Log file					*
	*****************************************************/
	//Write unused reads to file
	string unused_reads_log_file="unusedObsReads.txt";

	//cout<<"\nWrite unused reads to: \t"<<unused_reads_log_file<<endl;

	ofstream unused_reads;
	unused_reads.open(unused_reads_log_file.c_str());

	if(!unused_reads){
		cout<<"Unable to open" <<unused_reads_log_file<<endl;
		exit(1);
	}
	/****************************************************
	*					END Log file					*
	*****************************************************/
	
	ifstream readsBowtie;
	readsBowtie.open(rBowtie_file.c_str());
	if(!readsBowtie){
		cout<<"Unable to open "<<rBowtie_file.c_str()<<endl;
		exit(1);
	}
	
	
	int temp_counter=0;
	vector<string> current_read; //current read class
	string transcript_name; //name of the transcript where the read was mapped
	
	string read_name;
	string read_name_prev;
	bool first_read=true;
	
	string line;
	string field;
	while(readsBowtie.good()){
		getline(readsBowtie,line);
		
			//skip the headers
	 		while(line.substr(0,1)=="@")
	 			getline(readsBowtie,line);
		
		
		if((!line.empty()) && (line.substr(0,1)!="@")){

			istringstream iss(line);
			
			getline(iss,field,'\t');
				read_name=field.c_str();
			
			getline(iss,field,'\t');
			getline(iss,field,'\t');
				transcript_name=field.c_str();
				
					#if DEBUG
						cout<<"Read_name = "<<read_name<<endl;
						cout<<"transcript_name = "<<transcript_name<<endl;
					
						//exit(7);
					#endif
			
			if(transcript_name.compare("*")==0)
			{
					#if DEBUG
						cout<<"Line = "<<line<<endl;
					#endif
					
				//unused_reads<<line<<endl;
				continue;//skip to the next line
			}
			
			//first read
			if(first_read){
						#if DEBUG
							cout<<"\t1.The Very First Read: ";
						#endif
						
				first_read=false;

				current_read.push_back(transcript_name);
				//we assume that there are more than 1 reads in the file
			}else{
				
				if(read_name==read_name_prev){
					
						#if DEBUG
							cout<<"\t\t2.1. Same read"<<endl;
						#endif
					
					//Check if the transcript already exists in the vector
					if (find(current_read.begin(), current_read.end(),transcript_name)==current_read.end())
						current_read.push_back(transcript_name);//current_read is the current read class
						
				}//End Same Read
				else{
					#if DEBUG
						cout<<"\t\t2.2. New Read: "<<read_name<<endl;
						cout<<"Print current read classes"<<endl;
							print(current_read);
							
							//exit(7);
					#endif
						
						//first sort the current read class
						sort(current_read.begin(), current_read.end());
						
					//ADD PREVIOUS set of transcripts to the class and reset the current class (current_read)
						//also for previous class increment the correspnding transcript from where the current read came from
						//if value exists then do ++ otherwise assign 1
						
						if ( ReadsClasses.find(current_read) != ReadsClasses.end() ){
							ReadsClasses[current_read]++;
						}
						else{
							ReadsClasses[current_read]=1;
						}
						
							#if DEBUG
								cout<<"\t @@@ read "<<read_name_prev<<" was mapped to class [ ";
								for(const auto tr : current_read){
									cout<<tr<<" ";
								}
								cout<<"] @@@"<<endl;
								
								/*
								if (temp_counter == 2)
									exit(7);
								else
									temp_counter++;
								*/

								cout<<"\nPrint read_name_prev: "<<read_name_prev<<endl;
								cout<<"\nPrint current read:"<<endl;
									print(current_read);

								exit (7);
							#endif			
					
					//Reset the vector (this is the current reads class)
					current_read.clear();
					
					//add again the current transcript
					current_read.push_back(transcript_name);
										
				}//end: else: New Read: i.e. read_name!=read_name_prev
			}//end: else: not the first read in the file
			
			read_name_prev=read_name;			
			
		}//end: if((!line.empty()) && (line.substr(0,1)!="@"))
		
		if(line.empty()){
			//this shows that we got to the end of the file and we have to add the last class
							
							#if DEBUG
								cout<<"\t @@@ LAST read "<<read_name_prev<<" was mapped to class [";
								for(const auto tr : current_read){
									cout<<tr<<" ";
								}
								cout<<" ] @@@"<<endl;
								
							#endif
			
			//first sort the current read class (just to make sure it is sorted)
			sort(current_read.begin(), current_read.end());
			
			//if value exists then do ++ otherwise assign 1
			if ( ReadsClasses.find(current_read) != ReadsClasses.end() ){
				ReadsClasses[current_read]++;
			}
			else{
				ReadsClasses[current_read]=1;
			}
			
		}//end: if(line.empty())
	}//end: while(readsBowtie.good())
	
			#if DEBUG
				cout<<"\nextract_read_classes2: Print read classes"<<endl;
					print(ReadsClasses);
								
				//exit (7);
			#endif	

	
readsBowtie.close();

}