#include "load_obsCounts.h"
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
/*
Features:
	- this version uses vector<string> in the map for ReadsClasses instead of map<string,bool> same as:
	~/code/wt/d/include/extract_read_classes3.cpp (it is actually the next version of this code)

TO DO:
	-also save the discarded classes
using namespace std;
*/


//Extract read classes from Bowtie output
//This version is for observed reads only (reduced version of extract_read_classes2.cpp)

void load_obsCounts(string &rcCounts_file, map<vector<string>, int> &ReadsClasses)
{
	//This version uses vector<string> in the map for ReadsClasses instead of map<string,bool>
	
	cout<<"\n\textract_read_classes: "<<endl;
	cout<<"\t\tParameters: "<<endl;
	cout<<"\t\t\t1. Read Classes Counts File: "<<rcCounts_file<<endl;
	cout<<"\t\t\t2. Reads Classes map "<<endl;
	
		
	ifstream readsBowtie;
	readsBowtie.open(rcCounts_file.c_str());
	if(!readsBowtie){
		cout<<"Unable to open "<<rcCounts_file.c_str()<<endl;
		exit(1);
	}
	
	
	int temp_counter=0; //used for debugging
	
	vector<string> current_read; //current read class
	string transcript_name; //name of the transcript where the read was mapped
	
	string line;
	string field;
	while(readsBowtie.good()){
		getline(readsBowtie,line);

		
		if(!line.empty()){

			int counts=0;//reset counts reading
			
			istringstream iss(line);
			getline(iss,field,'\t'); //skip the first "["
				
				//if(temp_counter==2)
				//	exit(7);
				//temp_counter++;
				
			while(getline(iss,field,'\t'))
			{
				transcript_name=field.c_str();
				//cout<<"Transcript_name: "<<transcript_name<<endl;
				
				if(transcript_name.compare("]")==0)
				{
					//cout<<"End of list--no more transcripts"<<endl;
					getline(iss,field,'\t');//skip it and get the counts
					counts=atoi(field.c_str());
					//cout<<"Counts="<<counts<<endl;
				}
				else
				{
					current_read.push_back(transcript_name);
				}
				//cout<<"Field="<<field.c_str()<<endl;
			}
			
			//cout<<"Counts="<<counts<<endl;
			//exit(7);
			
			//fsort the current read class
			sort(current_read.begin(), current_read.end());
			
			//cout<<"Print Current Read Class"<<endl;
			//print(current_read);

			
			//Add the counts
			ReadsClasses[current_read]=counts;
			
			//Reset the vector (this is the current reads class)
			current_read.clear();	
			
		}//end: if(!line.empty())
	}//end: while(readsBowtie.good())
	
	
			#if DEBUG
				cout<<"\nextract_read_classes2: Print read classes"<<endl;
					print(ReadsClasses);
								
				//exit (7);
			#endif	

	
readsBowtie.close();

}

