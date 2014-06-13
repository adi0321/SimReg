//============================================================================
// Project     : 
// Author      : Adrian Caciula
// Description : Compute read coverage from sam file
// Created     : 04/06/2014
//
// Compile command:
// g++ -std=c++0x compute-coverage.cpp -o compute-coverage 
//============================================================================
//




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

#define DEBUG 1
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

void print(const map<int, int> & myMap);
vector<int> cigar(string , int );


int main(int argc,char *argv[]){
	if(argc<2){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

	//string grinder_fa = argv[1];
	string input_sam=argv[1];
	
	cout<<"\nRun "<<argv[0]<<endl;
	cout<<"\nInput File: "<<argv[1]<<endl;
	
	ifstream readsBowtie;
	readsBowtie.open(input_sam.c_str());
	if(!readsBowtie){
		cout<<"Unable to open "<<input_sam.c_str()<<endl;
		exit(1);
	}
	
	string output_file="coverage_profile.txt";
	//string m2_file=grinder_fa.substr(0,grinder_fa.find_last_of("."))+"_pair_2.fa";
	
		#if DEBUG
			cout<<"output_file = "<<output_file<<endl;
			//cout<<"m2_file = "<<m2_file<<endl;
				//exit(7);
		#endif
	
	ofstream out_m1;
	//ofstream out_m2;
	out_m1.open(output_file.c_str());
	//out_m2.open(m2_file.c_str());
	
	string line;
	string field;
	string current_read;
	string previous_read;
	bool read_p1=true;
	bool first_read=true;
	bool first_read_m2=true;
	
	map<int, int> positions; //this maps keeps the coverage for each position
	string transcript_name; //name of the transcript where the read was mapped
	
	while(readsBowtie.good()){
	
		vector<int> p1;
		p1.clear();
		int pos;
		
		getline(readsBowtie,line);
		
		//skip the headers
	 		while(line.substr(0,1)=="@")
	 			getline(readsBowtie,line);
		
		if((!line.empty()) && (line.substr(0,1)!="@")){
		
			istringstream iss(line);
			pos=0;
			
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			iss>>pos;
			
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			
			p1=cigar(field,pos);
			
				for (int i=0; i<p1.size(); i++)
				{
					positions[p1[i]]++;
					//cout<<p1[i]<<endl;
				}
			
		}//end if line not empty
	}

			//cout<<"Print map positions"<<endl;
			//print(positions);
			//cout<<"Done"<<endl;
			
			for (const auto& p : positions) {
				out_m1<<"[ "<<p.first<<" ] = " << p.second << endl;
			}
			
			
//out_m1<<"\n";
//out_m2<<"\n";

cout<<"Files had been created"<<endl;
cout<<"\t"<<output_file<<endl;
//cout<<"\t"<<m2_file<<endl<<endl;

//cout<<"NOTE: Check to see if both files have the same number of reads!"<<endl<<endl;
//This issue had been solved.

readsBowtie.close();
out_m1.close();
//out_m2.close();

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




//-----------------------------------------------------------
vector<int> cigar(string cigar, int start){
	 vector<int> c;//after spiting CIGAR
	 vector<int> coord;


	 stringstream stream(cigar);
	 string w1;
	 string w2;
	 while( getline(stream, w1, 'M') ){
		 stringstream ss(w1);
		 while( getline(ss, w2, 'N') )
			 c.push_back(atoi(w2.c_str()));
	 }

	 int previous=start-1;
	 for(int i=0;i<c.size();i++){
		 int start=previous+1;
		 for(int j=start;j<start+c[i];j++){
			 //cout<<j<<"	";
			 previous=j;
				if(i%2==0)
					coord.push_back(j);
		  }
	 }

	 return coord;
}

//--------------------------------------------------








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

void print(const map<int, int> & myMap) {

	for (const auto& kv : myMap) {
		cout <<"[ "<<kv.first<<" ] = " << kv.second << endl;
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