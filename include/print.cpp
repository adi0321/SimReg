#include "print.h"
#include <iostream> //cout
#include <vector>
#include <string>
#include <map>

using namespace std;

void print(const vector<string> &myVector){
	cout<<"\t[ ";
	for (const auto &kv : myVector) {
		cout<<kv<<" ";
	}
	cout<<"]"<<endl;
}

void print(const map<vector<string>, int> &myMap) {
	for (const auto &read_class: myMap){ //for each read class
			cout<<"[ ";
			for (const auto &tr : read_class.first)//for each transcript in this class
				cout<<tr<<" ";
			
			cout<<"] "<<read_class.second<<endl;
		
		
	}
}

void print(const map<vector<string>, double> &myMap) {
	for (const auto &read_class: myMap){ //for each read class
			cout<<"[ ";
			for (const auto &tr : read_class.first)//for each transcript in this class
				cout<<tr<<" ";
			
			cout<<"] "<<read_class.second<<endl;
		
		
	}
}

void print(const map<string, map<vector<string>,double> > &myMap) {
	for (const auto& it: myMap){ //for each transcript
		for (const auto &read_class : it.second) { //for each classes in this transcript
			cout<<it.first<<" "; //print current transcript name for each line
			cout<<"[ ";
			for (const auto &tr : read_class.first)
				cout<<tr<<" ";
			
			cout<<" ] "<<read_class.second<<endl;
		}
		
	}
}


/********************************************************************
 *      Prints d values	- new - with the component id	           	*
 ********************************************************************/
void print(const map<int, map<string, map<vector<string>,double> > > &myMap) {
for (const auto &comp : myMap){ //for each component
	for (const auto& it: comp.second){ //for each transcript
		for (const auto &read_class : it.second) { //for each classes in this transcript
			cout<<"Component "<<comp.first<<" - ";
			cout<<it.first<<" "; //print current transcript name for each line
			cout<<"[ ";
			for (const auto &tr : read_class.first)
				cout<<tr<<" ";
			
			cout<<" ] "<<read_class.second<<endl;
		}
		
	}
}
}


void print(const map<string, map<string, pair<double, double> > > &myMap) {
cout<<"\nPRINT ALL Normalized adjusted transcripts frequency on a line"<<endl<<endl;
	for (const auto &it : myMap){
	cout<<"Gene: "<<it.first<<" ("<<it.second.size()<<" transcript(s))"<<endl;
		for (const auto &kv : it.second) {
		cout <<"\tTranscript: "<<kv.first<<endl; 
				cout<<"\t\t"<<kv.second.first<<"\t"<<kv.second.second<<endl;
		}
	}
}

void print_reads_classes(const map< map<string, bool> , int> &myMap) {
cout<<"~~~~~ Print observed/virtual Reads Classes ~~~~~"<<endl;
cout<<"Read Class \t\t Size"<<endl;

	for (const auto &kv : myMap) {
		cout<<"[ ";
		for (const auto &cluster : kv.first){
				cout<<cluster.first<<" ";
		}
		cout<<" ] = "<< kv.second<<endl;	
	}

}

void print(const map<string, map< map<string, bool> , double> > &myMap) {
cout<<"~~~~~ Print p value ~~~~~"<<endl;
cout<<"Read Class \t Size"<<endl;

	for (const auto &tr : myMap){
		for (const auto &kv : tr.second) {
		cout<<"p["<<tr.first<<"][";
			for (const auto &cluster : kv.first){
				cout<<cluster.first<<" ";
			}
		cout<<"] = "<< kv.second<<endl;		
		}
	cout<<endl;
	}
}

void print(const map<string, bool> &myMap) {

	for (const auto& kv : myMap) {
		cout<<kv.first<<" = "<<kv.second<<endl;
	}
}

void print(const map<string, string> &myMap) {

	for (const auto& kv : myMap) {
		cout<<"\t\tRead "<<kv.first<<" comes from transcript "<<kv.second<<endl;
		
		//exit(7);
	}
}

void print(const map< map<string, bool> , map<string,int> > &myMap) {
cout<<"~~~~~ Print observed/virtual Reads Classes References ~~~~~"<<endl;
cout<<"Read Class \t\t References"<<endl;

	for (const auto &kv : myMap) {
		cout<<"[ ";
		for (const auto &cluster : kv.first){
				cout<<cluster.first<<" ";
		}
		cout<<"] = ";

		for (const auto &tr : kv.second){
				cout<<tr.first<<" = "<<tr.second<<"; ";
		}
		cout<<endl;
		
	}
}

//New print function for printing Reads Classes References
void print(const map< vector<string> , map<string,int> > &myMap) {
cout<<"~~~~~ Print MC Reads Classes References ~~~~~"<<endl;
cout<<"Read Class \t\t References"<<endl;

	for (const auto &kv : myMap) {
		cout<<"[ ";
		for (const auto &cluster : kv.first){
				cout<<cluster<<" ";
		}
		cout<<"] = ";

		for (const auto &tr : kv.second){
				cout<<tr.first<<" = "<<tr.second<<"; ";
		}
		cout<<endl;
		
	}
}


void print(const map< map<string, bool> , int> &myMap) {
cout<<"~~~~~ Print observed/virtual reads Classes ~~~~~"<<endl;
cout<<"Read Class \t\t Size"<<endl;

	for (const auto &kv : myMap) {
		cout<<"[ ";
		for (const auto &cluster : kv.first){
				cout<<cluster.first<<" ";
		}
		cout <<"] = "<< kv.second<<endl;	
	}

}