//============================================================================
// Project     : Life Technologies Collaborative Research Grant
// Author      : Serghei Mangul / Adrian Caciula
// Description : calculate d for the bipartite graph(transcripts and reads)
// Created     : 06/01/2013
//
// Compile command:
// g++ -std=c++0x ./compute_RC_CC_v2.cpp ../include/current_time.cpp ../include/norm.cpp ../include/print.cpp ../include/prepro_extract_read_classes_cc.cpp -o ./compute_RC_CC_v2 -I /usr/local/include/ -L /usr/local/lib/
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
"

double compare(vector<double> a, vector<double> b);
void gtf2junction(string gtf,
			vector<string> &tr_names, 
			map<string, map<string, int> > &transcript_len,
			map<string, map<string, vector<pair<int, int> > > > &exons_t,
			map<string, map<string,vector<int> > > &transcripts);
			
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
	if(argc<4){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

double EPS = 0.01;
int iterations=0;

	string gtf=argv[1]; //gtf
	string mcReadsFile=argv[2]; //raw Monte Carlo reads (grinder1-reads.fa for Grinder or true.sam for Marius simulator, respectively)
	string mcRBowtie_file=argv[3]; //aligned Monte Carlo reads
	
	vector<string> tr_names;
	map<string, map<string, int> > transcript_len;
	map<string, map<string, vector<pair<int, int> > > > exons_t;
	map<string, double> genes_frequency;
	
	cout<<"\nRun "<<argv[0]<<endl;
	cout<<"\n["<<current_time()<<"] Collecting Input Information ... "<<endl;

		//convert gtf to consecutive integers
		map<string, map<string, vector<int> > > transcripts; //transcript positions
	
		//~~~~~ get transcript names, and compute transcript lengths and transcript positions
		gtf2junction(gtf, tr_names, transcript_len, exons_t, transcripts); 

			#if DEBUG_T				
				cout<<"Transcripts lengths"<<endl;
					print(transcript_len);
			#endif
	
	

		cout<<"\n\tInput Values:"<<endl;
		cout<<"\t\tTotal number of Genes: "<<transcripts.size()<<endl;
		cout<<"\t\tTotal number of Transcripts: "<<tr_names.size()<<endl;
	cout<<"\n["<<current_time()<<"] Done!\n"<<endl;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cout<<"\n["<<current_time()<<"] Extracting Monte Carlo reads References..\n"<<endl;
	map<string, string> mcReadsRef;

	extract_reads_references(mcReadsFile, mcReadsRef);

cout<<"\n\t["<<current_time()<<"] Monte Carlo Reads References had been collected"<<endl;
		
		#if DEBUG_READ
			cout<<"Print MC Reads references"<<endl;
			print(mcReadsRef);
			
			exit(7);
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
			print(mcReadsClasses);
			cout<<endl;
			
		cout<<"MC Reads Classes References"<<endl;
			print(mcReadsClassesRef);
				//exit(7);
	#endif

cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;	
	
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
				
				//add 0.5 to get rid of zero values -->No need!!!
				//Many d_values have value zero (because reads had been mapped to those transcripts but the transcripts itself didn’t contribute to the class (i.e., no reads came from that transcript)

				//p_value_new2[tr][read_class.first]=mcReadsClassesRef[read_class.first][tr]+0.5;
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
						
					//exit(7);
				}
				
	
				//cout<<"mcReadsClassesRef[read_class.first]["<<tr.first<<"] = "<<mcReadsClassesRef[read_class.first][tr.first]<<endl;
				
				//exit(7);
				
				#endif
			}
		}
		
				#if DEBUG
					cout<<"\nPrint p_value_new2 BEFORE Normalization"<<endl;
					print(p_value_new2);
					
					//exit(7);
				#endif

		/*
		NOTE: Normalization was removed from the preprocessing step because we need the initial values in the 2nd phase in order to compute the number of reads per component
		Normalization per column will be done in the 2nd Phase (when the observed reads are computed in MCReg_CC_v1.cpp )
		10/8/2013: Note: No need for this reason---
			-- Normalization must be in 2nd Phase because the normalization must be done per component.
		
		
		//@@@NEW:
		//Normalize p values new2 - per column
		//for each transcript (column)
		for(const auto &tr: p_value_new2){
			
			//compute total sum
			double total_sum=0.0;
			
			//for each read class in the current transcript
			for (const auto &it : tr.second){
				total_sum+=it.second;
			}
			
			//cout<<"Total sum for transcript "<<tr.first<<" is: "<<total_sum<<endl;
			
			//Normalize
			//again for each read class in the current transcript
			for (const auto &it : tr.second){
				p_value_new2[tr.first][it.first]=it.second/total_sum;
			}

		}
		
				#if DEBUG
					cout<<"\nPrint Normalized p values (new2) - normalized by column"<<endl;
					print(p_value_new2);
				#endif
		*/
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done: d values had been computed"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;				

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/		

cout<<"\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Computing Connected Components ... "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

//Both DFS and BFS complete on O(n) time, using O(n) memory, where n is matrix size. But BFS it doesn't suffer from stack overflow problem, and it doesn't spend time on recursive calls.
//Anyway, let's first use boost library



   

	struct Vertex {
		std::string vertexName;
	};
	
	//typedef adjacency_list <vecS, vecS, undirectedS, Vertex> MyGraph;
	typedef boost::labeled_graph<adjacency_list<vecS, vecS, undirectedS, Vertex>, std::string> MyGraph;
	//The side effect of this is that one need to use graph() member function to make some algorithms work:
	
	typedef boost::graph_traits<MyGraph>::vertex_descriptor VertexID;
	
    MyGraph G;
	
	//for each transcript --  
	for (const auto &tr : p_value_new2){
		// Create vertices in that graph
		VertexID u = boost::add_vertex(tr.first,G);
		G[tr.first].vertexName = tr.first;
		
		//for each class in the current transcript
		for (const auto &read_class : tr.second ){
			
			string rClassName = "["; //add [ in order to differentiate the classes that contains only one transcripts from the transcripts itself
			
			//for each transcript in this class -- create a string with the name of the class
			for (const auto &tr_class : read_class.first)
				rClassName+=tr_class;
			
			//cout<<"Class name is: "<<rClassName<<endl;
			
			VertexID v = boost::add_vertex(rClassName,G);
			G[rClassName].vertexName = rClassName; //in case vertex already exists then this should overwrite
			
			add_edge(u, v, G);
		}
	}
	
	
    std::vector<int> component(num_vertices(G));
	
    int num = connected_components(G, &component[0]);
    cout << "Total number of components: " << num << endl;
	
	/*
	std::vector<int>::size_type i;
	for (i = 0; i != component.size(); ++i)
      cout << "Vertex " << i <<" is in component " << component[i]<<endl;
	  
    cout << endl;
	
	cout<<"\nPrint components: "<<endl;
	for (const auto &c : component) {
		cout<<"Component: "<<c<<endl;

	}
	
	MyGraph::vertex_iterator vertexIt, vertexEnd;
	boost::tie(vertexIt, vertexEnd) = vertices(G);
	for (; vertexIt != vertexEnd; ++vertexIt){
		VertexID vertexID = *vertexIt; // dereference vertexIt, get the ID
		Vertex & vertex = G.graph()[vertexID];
		//The side effect of boost::labeled_graph is that one need to use graph() member function to make some algorithms work
	
		cout<<"Vertex name is : "<<vertex.vertexName<<endl;
    }
	
	*/
	
	
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done: Connected Components had been computed"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
	
		
		
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
	
	//go thorugh components
	//cout<<"component size is: "<<component.size()<<endl;
	//CREATE MAP THAT WILL HOLD ALL THE COMPONENTS:
	map<int, vector<string>> cc;
	
	vector<int>::size_type i=0;
	Vertex & vertex = G.graph()[i];
	//cout<<"Vertex "<<i<<" has name: "<<vertex.vertexName<<endl;
	
	//Notes for each component (# of components == # of vertices but the component id changes only when we have a new component)
	for (i = 0; i < component.size(); i++){
		//cout<<"Component iterations: "<<i<<endl;
	
		Vertex & vertex = G.graph()[i];
		//cout<<"i="<<i<<" Element from component: "<<component[i]<<" is "<<vertex.vertexName<<endl;
		
		if(vertex.vertexName[0] != '[')//if the first character is not [
			cc[component[i]].push_back(vertex.vertexName);
	}//end: for each component
	
		#if DEBUG
			cout<<"Print component map:"<<endl;
			for(const auto &cmp : cc)
			{
				cout<<"component "<<cmp.first<<"containts: "<<endl;
					for(const auto &elem : cmp.second)
						cout<<"\t\t"<<elem<<endl;
			}
		#endif
	
//Now print to file
//For each component
for(const auto &cmp : cc)	
{//for each component there is a vector with transcript names

//for each transcript -- print transcript -- print class -- print d values
	for (const auto &tr : cmp.second){

		//for each class in the current transcript
		for (const auto &read_class : p_value_new2[tr] ){		
				d_stream_new2<<tr<<"\t[\t";

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
	
d_stream_new2<<"---"<<endl;
}//end:for each cc	
	
	
	
	
	
	/*
	//for each read class
	for (const auto &read_class : mcReadsClasses) {
		
		//for each transcript in T (T contains all transcripts)
		for (const auto &tr : p_value_new2)	
			d_stream_new2<<p_value_new2[tr.first][read_class.first]<<"\t";

		d_stream_new2<<"\n";

	}
	*/

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

//----------------------------------------------------------------------------------------------------------
//---- GTF2JUNCTION ----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
//open gtf and get all information
void gtf2junction(string gtf, 
		vector<string> &tr_names, 
		map<string, map<string, int> > &transcript_len,
		map<string, map<string,vector<pair<int, int> > > > &exons_t,
		map<string, map<string,vector<int> > > &transcripts){
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
		if(!line.empty()){
			istringstream iss(line);
			//cout<<line<<endl;
			getline(iss,field,'\t');
			//cout<<field<<endl;
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			getline(iss,field,'\t');
			int x=atoi(field.c_str());
			getline(iss,field,'\t');
			int y=atoi(field.c_str());
			getline(iss,field,'"');
			getline(iss,field,'"');
				gene_name=field.c_str();
				
			getline(iss,field,'"');
			getline(iss,field,'"');
				isoform=field.c_str();

					
			for(int i=x;i<=y;i++)
				transcripts[gene_name][isoform].push_back(i);
				
						#if DEBUG_G2J
							cout<<"\nGene name = "<<gene_name<<endl;
							cout<<"\tIsoform ----- Exons Start ----- Exon End"<<endl;
							cout<<"\t"<<isoform<<"-----"<<x<<"-----"<<y<<endl;
							
							
								//cout<<"Print transcripts"<<endl;
								//print(transcripts);
							cout<<"\tCurrent transcript length is "<<transcripts[gene_name][isoform].size()<<endl;
						#endif
			
			//first exon
			if(flag_first==0){
						#if DEBUG_G2J
							cout<<"\t1.The Very First exon: ";
						#endif
						
				flag_first=1;
				
				tr_names.push_back(isoform);
				current_length+=y-x+1;

				exons_t[gene_name][isoform].push_back(make_pair(x,y));
				
						#if DEBUG_G2J
							cout<<"exons_t["<<gene_name<<"]["<<isoform<<"] = "<<exons_t[gene_name][isoform][0].first<<"\t"
							<<exons_t[gene_name][isoform][0].second<<endl;
						#endif
			}else{
				
				if(isoform==isoform_prev){
					current_length+=y-x+1;
					transcript_len[gene_name][isoform]=current_length;
					
						#if DEBUG_G2J
							cout<<"\t\t2.1. Same transcript"<<endl;
							cout<<"\t\t --> current_length = "<<current_length<<endl;
							cout<<"\t\tstart="<<x<<" end="<<y<<" (prev exon end = "<<y_prev<<")"
							<<endl;
						#endif
						

					
					if(x!=y_prev+1){ //solve back-to-back exon problem
							//if x!=y_prev shows that we have another exon
								#if DEBUG_G2J
									cout<<"\n\t\t\t2.1.1. if(x!=y_prev+1) --> i.e., Not a back-to-back exon"<<endl;
								#endif
							
						exons_t[gene_name][isoform].push_back(make_pair(x,y));
							
								#if DEBUG_G2J
									//if(gene_name.compare("117_at")){
										cout<<"\t\t\tPrint Current Exons for gene "<<gene_name<<" transcript "<<isoform<<" ("<<exons_t[gene_name][isoform].size()<<" exons):"<<endl;
										print(exons_t[gene_name][isoform]); //Print Current Exons
									//}
								#endif

					}
					else{
								#if DEBUG_G2J
									cout<<"\t2.1.2 Else: x==y_prev+1 --> Back-to-back exon"<<endl;
									
									cout<<"\nGene name = "<<gene_name<<endl;
									cout<<"\tIsoform ----- Exons Start ----- Exon End"<<endl;
									cout<<"\t"<<isoform<<"-----"<<x<<"-----"<<y<<endl<<endl;
									
									
									cout<<"\t\t\tPrint Current Exons for transcript: "<<isoform<<" ("<<exons_t[gene_name][isoform].size()<<" exons):"<<endl;
										print(exons_t[gene_name][isoform]); //Print Current Exons
									
									//exit(7);
								#endif
						
						//update the last value in the vector
						int last_value=exons_t[gene_name][isoform].size()-1;
						exons_t[gene_name][isoform][last_value].second=y;
								#if DEBUG_G2J
										cout<<"\t\t\tPrint Current Exons for transcript: "<<isoform<<" ("<<exons_t[gene_name][isoform].size()<<" exons):"<<endl;
										print(exons_t[gene_name][isoform]); //Print Current Exons
									
										//exit(7);
								#endif
					}
				}//End Same transcript
				else{
					#if DEBUG_G2J
							cout<<"\t\t2.2. New transcript"<<endl;
							cout<<"\t\t\ti.e., isoform!=isoform_prev"<<endl;
							cout<<"\t\t\tcurrent_length = "<<current_length<<endl;
						#endif
					

					
					tr_names.push_back(isoform);
					current_length=0; //reset previous length
					current_length+=y-x+1;
					transcript_len[gene_name][isoform]=current_length;
										
					exons_t[gene_name][isoform].push_back(make_pair(x,y));
					
								#if DEBUG_G2J
									cout<<"\t\t\tPrint Current Exons for gene "<<gene_name<<" transcript "<<isoform<<" ("<<exons_t[gene][isoform].size()<<" exons):"<<endl;
										
									print(exons_t[gene_name]); //Print Current Exons
								#endif
					
				}//End New Transcript
			}

			x_prev=x;
			y_prev=y;
			isoform_prev=isoform;
		}
	}

	#if DEBUG_G2J
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