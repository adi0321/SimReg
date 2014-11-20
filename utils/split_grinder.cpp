//============================================================================
// Project     : Life Technologies Collaborative Research Grant
// Author      : Adrian Caciula / Serghei Mangul /
// Description : Split reads from grinder into 2 files: pair_1.fa and pair_2.fa
// Created     : 06/01/2013
//
// Compile command:
// g++ -std=c++0x -o split_grinder split_grinder.cpp ../include/norm.cpp
//============================================================================
//


#include "/home/adrian/MCReg_v15/include/norm.h"

#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

#include <map>
#include <iterator>
#include <random>

#define DEBUG 0
#define DEBUG_GF 0 //debug for genes frequency

using namespace std;

#define HELPMESSAGE "too few parameters: \n\
coverage \n\
Input:\n\
argv[1] - grinder_reads.fa \n\
"

double compare(vector<double> a, vector<double> b);
void get_names(string gtf, map<string, vector<string> > &gt_names);
			
void print(const vector<pair<int,int> > &myVector);
void print(const map<string, int> &myMap);
void print(const map<string, double> &myMap);
void print(const map<string, bool > & myMap);
void print(const map<string, map<string, int> > &myMap); //Print Transcripts Length
void print(const map<string, map<string, double> > &myMap); //Print Normalized Adjusted Transcripts Length
void print(const map<string, vector<pair<int, int> > > &myMap);
void print(const map<string, map<string, vector<pair<int, int> > > > &myMap); //print all exons
void print(const map<string, vector<int> > &myMap);
void print(const map<string, map<string, vector<int> > >&myMap);
void print(const map<vector<int>, int > &myMap);
void print(const map<vector<int> , vector<pair<int, int> > > &myMap);
void print(const map< map<string, bool> , int> &myMap);
void print(const map< map<string, bool> , double> &myMap);
void print(const map<string, map< map<string, bool> , int> > &myMap);
void print(const map<string, map< map<string, bool> , double> > &myMap);

//TO DO: Write read m1 only if read m2 exists: Done (11/8/2013)
//if read m1 do not have a corresponding read m2 then bowtie returns an error during alignment

int main(int argc,char *argv[]){
	if(argc<2){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

	string grinder_fa = argv[1];
	
	cout<<"\nRun "<<argv[0]<<endl;
	cout<<"\nInput File: "<<argv[1]<<endl;
	cout<<"\nSplitting File ... "<<endl;
	
	
	
	// ~~~~~~ Open reads file ~~~~~
	ifstream input_stream;
	input_stream.open(grinder_fa.c_str());
	if(!input_stream){
		cout<<"Unable to open "<<grinder_fa.c_str()<<endl;
		exit(1);
	}
	
	string m1_file=grinder_fa.substr(0,grinder_fa.find_last_of("."))+"_pair_1.fa";
	string m2_file=grinder_fa.substr(0,grinder_fa.find_last_of("."))+"_pair_2.fa";
	
		#if DEBUG
			cout<<"m1_file = "<<m1_file<<endl;
			cout<<"m2_file = "<<m2_file<<endl;
				//exit(7);
		#endif
	
	ofstream out_m1;
	ofstream out_m2;
	out_m1.open(m1_file.c_str());
	out_m2.open(m2_file.c_str());
	
	
	string line;
	string field;
	string current_read;
	string previous_read;
	bool read_p1=true;
	bool first_read=true;
	bool first_read_m2=true;
	
	string read_m1; //save read m1 until we make sure we also have read m2
	//if read m1 do not have a corresponding read m2 then bowtie returns an error during alignment
	
	//Print also the count (How many reads were generated from each transcript)
	map<string, int> reads_count; //for each transcript
	string ref;
	
	while(input_stream.good()){
		getline(input_stream,line);
		if(!line.empty()){
			istringstream iss(line);
			
			getline(iss,field,'/');
			char first_char=field.at(0);
			
						#if DEBUG
							cout<<"Process Field: "<<field<<endl;
						#endif
			
			if (first_char =='>'){ //new read
				
				current_read=field.c_str();
									
				if(first_read){ //enter here only the very first time
					//out_m1<<current_read<<"\n";
					getline(iss,field,' '); //get the rest of the line after read name/id
					getline(iss,field); //rest of the line is in field
					read_m1=current_read+" "+field+"\n";
					first_read=false;
					
					//get the reference for counting
					istringstream iss(field);
					getline(iss,field,'='); //remove the reference string in front of =	
					getline(iss,field,' ');
						//cout<<field<<endl;
					ref=field;
					reads_count[ref]=1; //we assign 1 because this is the first count
					
						#if DEBUG
							cout<<"\tFirst Read: "<<current_read<<" -> write to m1 file"<<endl;
							cout<<"Read m1 is: "<<read_m1<<endl;
							cout<<"Reference: "<<ref<<endl;
							
							//exit(7);
						#endif
				}
				else
				{
							#if DEBUG
								cout<<"\tRead: "<<current_read<<endl;
								cout<<"\t\tfirst_char: "<<first_char<<endl;
								cout<<"\t\tPrevious Read: "<<previous_read<<endl;
								//previous read is assigned first read at the end: previous_read = current_read;
								cout<<"\t\tvs. Current Read: "<<current_read;
								cout<<" = "<<current_read.compare(previous_read)<<endl;
								//exit(7);
							#endif
							
					if( (current_read.compare(previous_read)) == 0 ){//if equal zero then this is m2
						
							#if DEBUG
								cout<<"\tPair Read! (first_read_m2 = "<<first_read_m2<<" )"<<endl;
								cout<<"\t\tRead: "<<current_read<<" -> write to m2 file"<<endl;
							#endif
					
						if(!first_read_m2){//add new line for all new reads in m2 (except the first)
								#if DEBUG
									cout<<"\t\tAdd new line in m2"<<endl;
								#endif
							out_m2<<"\n";
						}
							first_read_m2=false;//set to false in order to enter next time
						
						getline(iss,field,' '); //get the rest of the line after read name/id
						getline(iss,field); //rest of the line is in field
						out_m2<<current_read<<" "+field+"\n";
						
						read_p1=false;
						
						//Now we can also write m1 to file
						out_m1<<read_m1;
						
						//get the reference for counting
						istringstream iss(field);
						getline(iss,field,'='); //remove the reference string in front of =	
						getline(iss,field,' ');
						ref=field;
						reads_count[ref]++; //increment the counter for this particular transcript
						
					}else{//else this is m1
							#if DEBUG
								cout<<"\tNew Read: "<<current_read<<" -> write to m1 file"<<endl;
							#endif
					
						//add new line in m1 before each new read
						//out_m1<<"\n"<<current_read<<"\n";
						
						getline(iss,field,' '); //get the rest of the line after read name/id
						getline(iss,field); //rest of the line is in field
						read_m1="\n"+current_read+" "+field+"\n";
						read_p1=true;
						
						//get the reference for counting
						istringstream iss(field);
						getline(iss,field,'='); //remove the reference string in front of =	
						getline(iss,field,' ');
						ref=field;
						reads_count[ref]++; //increment the counter for this particular transcript
					}
				}
			}
			else
			{//same read
						#if DEBUG
							cout<<"\t\tSame Read: "<<current_read<<endl;
							cout<<"\t\t\tfirst_char = "<<first_char<<endl;
						#endif
			
				if(read_p1){
					//out_m1<<field;
					read_m1+=field;
							#if DEBUG
								cout<<"\t\tWrite Read: "<<current_read<<" -> to m1 file"<<endl;
							#endif
				}
				else{
					out_m2<<field;
							
							#if DEBUG
								cout<<"\t\tWrite Read: "<<current_read<<" -> to m2 file"<<endl;
							#endif
				}
				
				previous_read = current_read;
				
						#if DEBUG
							cout<<"\t\tPrevious read: "<<previous_read<<endl;
						#endif
				

			}
		}
	}

out_m1<<"\n";
out_m2<<"\n";

cout<<"Files had been created"<<endl;
cout<<"\t"<<m1_file<<endl;
cout<<"\t"<<m2_file<<endl<<endl;

#if DEBUG
cout<<"Read Count Summary:"<<endl;
for(const auto &t : reads_count)
	cout<<t.first<<": "<<t.second<<endl;
#endif
	
cout<<endl<<"Done Spliding Grinder FA file!"<<endl<<endl;

//cout<<"NOTE: Check to see if both files have the same number of reads!"<<endl<<endl;
//This issue had been solved.

input_stream.close();
out_m1.close();
out_m2.close();

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
//---- GTF2JUNCTION ----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
//open gtf and get all information
void get_names(string gtf, map<string, vector<string> > &gt_names)
{

	#define DEBUG_GN 0
	
	string line;
	string field;
	string gene_name;
	string isoform;
	string isoform_prev;
	ifstream ingtf;
	int flag_first=0;


	ingtf.open(gtf.c_str());
	if(!ingtf){
		cout<<"Unable to open "<<gtf.c_str()<<endl;
		exit(1);
	}

	int current_length=0;
	
	while(ingtf.good()){
		getline(ingtf,line);
		if(!line.empty()){
			istringstream iss(line);
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'"');
			getline(iss,field,'"');
				gene_name=field.c_str();
				
			getline(iss,field,'"');
			getline(iss,field,'"');
				isoform=field.c_str();

			
			//first exon
			if(flag_first==0){
						#if DEBUG_GN
							cout<<"\t1.The Very First exon: ";
						#endif
						
				flag_first=1;
				
				gt_names[gene_name].push_back(isoform);
				
			}else{
				
				if(isoform!=isoform_prev){
					#if DEBUG_GN
							cout<<"\t\tNew transcript"<<endl;
							cout<<"\t\t\ti.e., isoform!=isoform_prev"<<endl;
						#endif

					gt_names[gene_name].push_back(isoform);
					
				}//End New Transcript
			}

			isoform_prev=isoform;
		}
	}

	#if DEBUG_GN
		cout<<"~~~~~~~~~~~~~~~~~"<<endl;
		cout<<"~~~~~ Done! ~~~~~"<<endl;
		cout<<"~~~~~~~~~~~~~~~~~"<<endl<<endl;
	#endif
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

void print(const map<string, bool > & myMap) {

	for (const auto& kv : myMap) {
		cout<<kv.first<<" = "<<kv.second<<endl;
		
		//exit(7);
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

void print(const map< map<string, bool> , int> &myMap) {
cout<<"~~~~~ Print observed/virtual reads ~~~~~"<<endl;
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

void print(const map<string, map< map<string, bool> , double> > &myMap) {
cout<<"~~~~~ Print p value ~~~~~"<<endl;
cout<<"Read Class \t Size"<<endl;

	for (const auto &tr : myMap){
		for (const auto &kv : tr.second) {
		cout<<"p["<<tr.first<<"][";
			for (const auto &cluster : kv.first){
				cout<<cluster.second<<" ";
			}
		cout<<"] = "<< kv.second<<endl;		
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
//----------------------------------------------------------------------------------------------------------