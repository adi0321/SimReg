#include "prepro_extract_read_classes_cc.h"
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
// this version uses vector<string> in the map for ReadsClasses instead of map<string,bool> same as:
//~/code/wt/d/include/extract_read_classes3.cpp (it is actually the next version of this code)

using namespace std;

//Extract read classes from Bowtie output
//This version is for observed reads only (reduced version of extract_read_classes2.cpp)

void prepro_extract_read_classes_cc(
						string &rBowtie_file, 
						map<vector<string>, int> &ReadsClasses,
						map<string, string> &readsRef,
						map<vector<string>, map<string, int> > &readsClassesRef)
{
	//This version uses vector<string> in the map for ReadsClasses instead of map<string,bool>
	
	cout<<"\n\textract_read_classes: "<<endl;
	cout<<"\t\tParameters: "<<endl;
	cout<<"\t\t\t1. Bowtie file: "<<rBowtie_file<<endl;
	cout<<"\t\t\t2. Reads Classes map "<<endl;
	
	/****************************************************
	*					Error Log file					*
	*****************************************************/
	//Write unused reads to file
	string unused_reads_log_file=rBowtie_file.substr(0,rBowtie_file.find_last_of("/"))+"/unused_reads.txt";

	cout<<"\nWrite unused reads to: \t"<<unused_reads_log_file<<endl;

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
					
				unused_reads<<line<<endl;
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
							#endif

						
						//Find current read class and increment its count
						if ( readsClassesRef.find(current_read) != readsClassesRef.end() ){
							readsClassesRef[current_read][readsRef[read_name_prev]]++;
						}
						else{
							readsClassesRef[current_read][readsRef[read_name_prev]]=1;
						}
							
							#if DEBUG
								cout<<"\nPrint read_name_prev: "<<read_name_prev<<endl;
								//cout<<"\nreadsRef = "<<endl;
									//print(readsRef);
								cout<<"\nPrint current read:"<<endl;
									print(current_read);
								cout<<"\nPrint read classes references"<<endl;
								print(readsClassesRef);
								
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
			
			
			if ( readsClassesRef.find(current_read) != readsClassesRef.end() ){
				readsClassesRef[current_read][readsRef[read_name_prev]]++;
			}
			else{
				readsClassesRef[current_read][readsRef[read_name_prev]]=1;
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



/********************************************************************
*			Extract reads references from Grinder output			*
*********************************************************************/
void extract_reads_references(string &rGrinder_file, map<string, string> &readRef)
{

#define DEBUG_REF 0
//TO DO: Adjust the program to read the read references from marius simulator program
	cout<<"\n\textract_reads_references: "<<endl;
	cout<<"\t\tInput Parameters: "<<endl;
	cout<<"\t\t\t1. Grinder file: "<<rGrinder_file<<endl;
	cout<<"\t\t\t2. Reads Classes map "<<endl;
	
	bool simReg=true; //flag used if SimReg simulator is used
	
	ifstream readsGrinder;
	readsGrinder.open(rGrinder_file.c_str());
	if(!readsGrinder){
		cout<<"Unable to open "<<rGrinder_file.c_str()<<endl;
		exit(1);
	}
	
	
	string line;
	string field;
	string current_read;
	string previous_read;
	string refTr; //reference Transcript
	bool first_read=true;
	bool first_read_m2=true;
	
	string read_type="none";
	
	while(readsGrinder.good()){
		getline(readsGrinder,line);
		if(!line.empty()){
			istringstream iss(line);
			
						#if DEBUG_REF
							cout<<"Parse line: "<<line<<endl;
						#endif
			
			if(first_read)
			{//only for the first line ... see if those are Marius simulator's reads
				//first_read will become false few lines below
				
				if (line.compare("@HD	VN:1.0	SO:unsorted") == 0)
				{
					cout<<"Those are reads from Marius Simulator"<<endl;
				
					getline(readsGrinder,line);//get the next line
					istringstream iss(line);
					
						#if DEBUG_REF
							cout<<"\tNew line is: "<<line<<endl;
						#endif
					
					read_type="Marius";
				}
				else
					read_type="Grinder";
			}
			
			
			if(read_type.compare("Grinder") == 0)
			{
				if(simReg)
					getline(iss,field,' '); //SimReg Simulator
				else
					getline(iss,field,'/'); //Grinder Simulator
				
			}
			else{
				istringstream iss(line);
				getline(iss,field,'\t');//Marius simulator's read
			}				

			char first_char=field.at(0);
			
								#if DEBUG_REF
									cout<<"read_type = "<<read_type<<endl;
									cout<<"line = "<<line<<endl;
									cout<<"field = "<<field<<endl;
									cout<<"First char  = "<<first_char<<endl;
								#endif
			
			if (first_char =='>' || first_char == 'r'){ //new read: '>' for grinder and 'r' for Marius simulator reads
				
				current_read=field.c_str();
				
				if (first_char == 'r'){ //new read from Marius
					
					string ref;
					istringstream iss(field);
					getline(iss,refTr,'_');
					getline(iss,refTr);//no need for delimiter since this field is composed of only two words
				}
				else
				{//new read from Grinder
					
					current_read.erase(0,1); //from position 0 erase one character (the '>')
				
					if(!simReg)
						getline(iss,field,' ');
					
					getline(iss,field,'=');
					getline(iss,field,' ');
						refTr=field.c_str();
				}
						#if DEBUG_REF
							cout<<"\tread "<<current_read<<" from reference Transcript = "<<refTr<<endl;
						#endif
				
				if(first_read){ //enter here only the very first time
					readRef[current_read]=refTr;
					first_read=false;
					
						#if DEBUG_REF
							cout<<"\tThe Very First Read: "<<current_read<<" -> add to readRef"<<endl;
							print(readRef);
								//exit(7);
						#endif
				}
				else
				{
							#if DEBUG_REF
								cout<<"\tRead: "<<current_read<<endl;
								cout<<"\t\tfirst_char: "<<first_char<<endl;
								cout<<"\t\tPrevious Read: "<<previous_read<<endl;
								cout<<"\t\tvs. Current Read: "<<current_read;
								cout<<" = "<<current_read.compare(previous_read)<<endl;
							#endif
							
					if( (current_read.compare(previous_read)) == 0 ){//if equal zero then this is m2
						//so we can ignor it
							
						first_read_m2=false;//set to false in order to enter next time
						
					}else{//else this is m1
						readRef[current_read]=refTr;

							#if DEBUG_REF
								cout<<"\tNew Read: "<<current_read<<" -> add to readRef"<<endl;
								print(readRef);
									//exit(7);
							#endif
					}
				}
				
				if(read_type.compare("Marius") == 0)
				{//for Marius simulator's each read is on a single line --> i.e., it will never go in the else condition below
								#if DEBUG_REF
									cout<<"Read type = "<<read_type<<endl;
								#endif
								
					previous_read = current_read;
				}
				
			} //end new read
			else
			{//same read --> i.e., we can skip this part
				
				previous_read = current_read;
				
						#if DEBUG_REF
							cout<<"\t\tPrevious read: "<<previous_read<<endl;
						#endif
			}

		}//end if line not empty
		
	}
readsGrinder.close();
}