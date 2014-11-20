//============================================================================
// Author      : Adrian Caciula
// Description :-calculate observed reads classes for the bipartite graph(transcripts and reads)
//				-Also writes d values to proper format (required by the solver ~/code/wt/regression/qp.py)
//
// Created     : 08/20/2013
//
// Compile command:
// g++ -std=c++0x ./MCReg_CC_v2.cpp ../include/current_time.cpp ../include/print.cpp ../include/extract_obsCounts.cpp -o ./MCReg_CC_v2
//============================================================================

//This version is an extension of ~/code/wt/d/MCReg_CC_v1.cpp

/***
Features
	//10/15/2013: add option for connectied components 0 - without cc and 1 with cc
	//10/10/2013: add option for relative least squares min(d/o - 1)^2 or regula least squares min(d-o)^2 
	//Default Features:
		- read d file component by component
		- uses vector<string> in the map for ReadsClasses instead of map<string,bool>

	// The MC Reads Classes and d values had already by computed by:
		// ~/code/MCReg/annotation_preprocessing/compute_RC_CC_v2
		//	An example of pre-processed annotation for chr1 can be found in /data1/adrian/MCReg/preprocess_hg19_chr1
	- d values contain now a new parameter s.t. when it loads the components are preserved
***/

#include "../include/current_time.h"
#include "../include/print.h"
#include "../include/extract_obsCounts.h"

#include <iostream>     // std::cout
#include <fstream>		// Input/output stream class to operate on files.
#include <sstream>      // std::stringstream //Stream class to operate on strings.
#include <string>       // std::string
#include <vector>
#include <algorithm>
#include <map>

#define DEBUG 0
#define DEBUG_D 0
#define DEBUG_RC 1

using namespace std;

#define HELPMESSAGE "too few parameters: \n\
Input:\n\
argv[1] - d_values file (from the pre-processed directory)\n\
			An Example of the expected file is:\n\
			/data1/adrian/MCReg/preprocess_hg19_chr1/knownGeneGnfAtlas2_chr1_no_DUMMYCLUSTER_d_values.txt\n\
argv[2] - MC reads_classes file (from the same pre-processed directory as argv[1])\n\
argv[3] - File with observed reads classes counts \n\
argv[4] - Least Squares Type: 0 - relative least squares min(d/o - 1)^2 ; 1 - for regular least squares min(d-o)^2 \n\
argv[5] - Connected Components Option: 0 - without CC ; 1 - with CC \n\
"

int main(int argc,char *argv[]){
	if(argc<6){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}
	
	cout<<"\nRunning "<<argv[0]<<" ..."<<endl;
	
	string d_values_file=argv[1];
	string mcReadsClasses_file=argv[2];
	string obsRCcounts_file=argv[3]; 
	int least_square_type=atoi(argv[4]);
	//0 - for relative least squares min(d/o - 1)^2 or 1 - for regular least squares min(d-o)^2 
	
	int cc_option=atoi(argv[5]);
	
/********************************************************************	
Step1: Parse d_values file (in order to divide by o values 
(the correct format, required by the solver, will be done by a sh script)		*
//new edit: just ignore lines that start with "---"
*********************************************************************/
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Step1: Parse d_values file "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

	ifstream d_values_stream;
	d_values_stream.open(d_values_file.c_str());
	if(!d_values_stream){
		cout<<"Unable to open "<<d_values_file.c_str()<<endl;
		exit(1);
	}

	string line;
	string field;
	string tr;
	string tr_class; //transcript name from the class
	double value=0.0;
	map<int, map<string, map<vector<string>,double> > > d_values; 
	int component_id=0; 
	
	vector<string> mcreads_class;
	
	while(d_values_stream.good()){
		getline(d_values_stream,line);
		
			if(line.substr(0,1)=="-"){
	 			component_id++;
				getline(d_values_stream,line);//go to the next line
			}
		
		if(!line.empty()){
			istringstream iss(line);
			
			getline(iss,field,'\t');
			tr=field.c_str();
				//cout<<"Transcript name is: "<<tr<<endl;
			getline(iss,field,'\t'); //skip the '['
				//cout<<"2nd field is: "<<field.c_str()<<endl;

			while (getline(iss,field,'\t')){
					//cout<<"Next field is: "<<field.c_str()<<endl;
				tr_class=field.c_str();
				if (tr_class.compare("]") != 0)
					mcreads_class.push_back(field.c_str());
				else{
					getline(iss,field,'\t'); //skip the tab space after ']'
					//getline(iss,field,' '); //get the value
					value = atof(field.c_str());
						//cout<<"Value = "<<value<<endl;
				}
			}
			//They are already sorted -- check later to see if this step is still necessary
			//sort(mcreads_class.begin(), mcreads_class.end());
				
				//cout<<"Print transcripts from the vector class"<<endl;
				//print(mcreads_class);
				
				d_values[component_id][tr][mcreads_class]=value;
				mcreads_class.clear(); //clear vector
			//exit(7);
		}
	}
	
		#if DEBUG_D
			cout<<"Print d values: "<<endl;
			print(d_values);
			
			//exit(7);
		#endif
		
d_values_stream.close();

cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
/*******************************************************************************************************
*	END	Step 1	*
********************************************************************************************************/	



/********************************************************	
Step2: Load MC read classes from the pre-processed file	*
*********************************************************/
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Step2: Load MC read classes from the pre-processed file "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

	ifstream mcReadsClasses_stream;
	mcReadsClasses_stream.open(mcReadsClasses_file.c_str());
	if(!mcReadsClasses_stream){
		cout<<"Unable to open "<<mcReadsClasses_file.c_str()<<endl;
		exit(1);
	}

	line.clear();
	field.clear();
	tr.clear();
	mcreads_class.clear(); //clear vector //this is current read class
	map<vector<string>, int> obsReadsClasses;
	
	while(mcReadsClasses_stream.good()){
		getline(mcReadsClasses_stream,line);
		if(!line.empty()){
			istringstream iss(line);
			
			getline(iss,field,'['); //skip the '['
			getline(iss,field,' '); //skip the space after '['
				//cout<<"3rd field is: "<<field.c_str()<<endl;
			while (getline(iss,field,' ')){
					//cout<<"Next field is: "<<field.c_str()<<endl;
				tr=field.c_str();
				if (tr.compare("]") != 0)
					mcreads_class.push_back(field.c_str());
				else
					break;	
			}
			
			obsReadsClasses[mcreads_class]=0; //just initialize all with 0
			mcreads_class.clear(); //clear vector
		}
	}
	
	mcReadsClasses_stream.close();
	
		#if DEBUG_RC
			cout<<"Print MC Reads Classes: - just classes names "<<endl;
			print(obsReadsClasses);
		#endif
	cout<<"Monte Carlo Reads classes had been loaded into memory! - current size: "<<obsReadsClasses.size()<<endl;

cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
/*******************************************************************************************************
*	END	Step 2	*
********************************************************************************************************/	
	
	
	
/************************************************************	
*	Step3: Extract observed read classes from Bowtie output	*
*************************************************************/
cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] Step3: Load observed reads counts "<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	

//here nothing to extract ....just load the values from file obsRCcounts.txt ... compare and see if those classes exists in obsReadsClasses .... otherwise discard those counts
//2014/5/1 - We need to do the union of read classes
int discarded_reads=extract_obsCounts(obsRCcounts_file, obsReadsClasses);
//TO DO: Return a pair<int, int> where pair.first is total number of observed reads
//and pair.sedond is the number of discarded reads
cout<<"Done extract_obsCounts"<<endl;


//Compute Total number of observed reads
int total_obs_reads=0;
for (const auto &r_class : obsReadsClasses){
	total_obs_reads+=r_class.second;
}

cout<<"Total number of discarded reads is: "<<discarded_reads<<endl;
cout<<"~~~ Remining reads = "<<total_obs_reads<<" reads"<<endl;
cout<<"Total number of Observed Reads Classes is: "<<obsReadsClasses.size()<<endl;	

exit(7);
		#if DEBUG_RC
			cout<<"\nOBSERVED Reads Classes \t Size:"<<endl;
				print(obsReadsClasses);
				exit(7);
		#endif

cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"["<<current_time()<<"] ~~~~~ Done! ~~~~~"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
/*******************************************************************************************************
*	END	Step 3	*
********************************************************************************************************/

////2014/5/1 - We need to write again the file with obsRCSize
//Write to file the number of transcripts in each readClass
ofstream obsRCsize_stream;
string outFile="./obsRCsize.txt"; 

obsRCsize_stream.open(outFile.c_str());
	if(!obsRCsize_stream){
		cout<<"Unable to open" <<outFile<<endl;
		exit(1);
	}

	//Matrix for size of read class
	vector<vector<double>> rcSize;
	//initialize matrix rcSize with zero
	for(const auto &rc : obsReadsClasses)
    {
		vector<double> row; // Create an empty row
         for(int j=0;j<obsReadsClasses.size();j++)
         {
              row.push_back(0.0);   // Add an element (column) to the row
         }
		rcSize.push_back(row); // Add the row to the main vector		 
	}

	//fill the diagonal with read classes size.
	//for each read class ID
	int i=0;
	for(const auto &rc : obsReadsClasses)
	{
		if(rc.second != 0)
			rcSize[i][i]=(double)1/rc.first.size();
		else
			rcSize[i][i]=0;
			
		i++;
	}

	//print Matrix to file
	for(int i=0;i<obsReadsClasses.size();i++)
    {
         for(int j=0;j<obsReadsClasses.size();j++)
         {
				obsRCsize_stream<<rcSize[i][j]<<"\t";
         }   
         obsRCsize_stream<<endl;
	}	

	obsRCsize_stream.close();	



/************************************************************	
*	Step4: Compute Observed Reads Classes FREQUENCY	*
*************************************************************/

cout<<"\nCompute Observed Reads Classes FREQUENCY:"<<endl;
map <vector<string>, double> new_o;
	for (const auto &reads_it : obsReadsClasses){
		double temp_o=(double)1/total_obs_reads;
		//we normalize but all observed reads (including those that were deleted)

		new_o[reads_it.first]=temp_o*reads_it.second;
			
			#if DEBUG_RC
				cout<<"new_o = 1/total_obs_reads * class_size = 1/"
				<<total_obs_reads<<" * "<<reads_it.second<<" = "<<temp_o<<" * "<<reads_it.second<<" = "
				<<new_o[reads_it.first]<<endl;				
			#endif
	}

#if DEBUG_RC
	cout<<"Print new_o"<<endl;
	print(new_o);
	
	cout<<"\nOBSERVED Reads Classes \t Size:"<<endl;
		print(obsReadsClasses);
	//exit(7);
#endif

cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"Observed Classes size: "<<obsReadsClasses.size()<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;
/*******************************************************************************************************
*	END	Step 4	*
********************************************************************************************************/	
	
/************************************************************	
*	Step5: Compute d values	*
*************************************************************/
				#if DEBUG_D
					cout<<"\n\nPrint d values before normalization "<<endl;
					print(d_values);
					
					//exit(7);
				#endif

		//@@@NEW:
		map<int, map<string, map<vector<string>,double> > > d_value_new3;
		//Normalize p values new2 - per column
		
		//Normalization had been moved here from the 1st phase (preprocesing step)
		//because we needed the values here in order to compute the number of reads per component
		//--no need anymore for this scope
		//Update: Normalization it stays here cause it needs to be done per component ---
		
	for(const auto &comp: d_values){//for each component
	
		//for each transcript (column) in the current component
		for(const auto &tr: comp.second){		
			
			//compute total sum
			double total_sum=0.0;
			
			//for each read class in the current transcript
			for (const auto &it : tr.second)
			{
				total_sum+=it.second;
				
					#if DEBUG
						if(comp.first==1018)
						{
							cout<<"\nRead class in current transcript:"<<endl;
								print(tr.second);
							cout<<"New value it.second="<<it.second<<endl;
							cout<<"Current sum for transcript "<<tr.first<<" is: "<<total_sum<<endl;
						}
					#endif
			}
			
			#if DEBUG
			if(comp.first==1018)
			{
				cout<<"Total sum for transcript "<<tr.first<<" is: "<<total_sum<<endl;
			}
			#endif
			
			//Normalize
			//again for each read class in the current transcript
			for (const auto &it : tr.second)
			{
						
				if(total_sum==0)
					d_values[comp.first][tr.first][it.first]=0;
				else
					d_values[comp.first][tr.first][it.first]=it.second/total_sum;
						
						#if DEBUG
							if(comp.first==1018)
							{
								cout<<"d_values[comp.first][tr.first][it.first] = "<<d_values[comp.first][tr.first][it.first]<<endl;
								cout<<"it.second = "<<it.second<<endl;
								cout<<"total_sum = "<<total_sum<<endl;
							}
						#endif
			}
			
				#if DEBUG
					if(comp.first==1018)
					{
					cout<<"\nPrint Normalized d values for component: "<<comp.first<<endl;
					print(d_values[comp.first]);
					
					//exit(7);
					}
					
				#endif
			
			/*
			//Note: 10/8/2013 Skip this for now - Don't forget to modify in the printing to file section
				//Now you need to print d_values and not d_value_new3 -- done (10/10/2013)
			Update: 10/10/2013 - option for relative or not relative had been created
			*/
			
			if(least_square_type == 0)
			{
			//@@@NEW ~~~ Now divide by the corresponding observed value
			for (const auto &it : tr.second){
						#if DEBUG_D
							cout<<"new_o[it.first] = "<<new_o[it.first]<<endl;
							cout<<"print it.first"<<endl;
							print(it.first);
						#endif
				
				if(new_o[it.first]==0)
					d_value_new3[comp.first][tr.first][it.first]=0;
				else
					d_value_new3[comp.first][tr.first][it.first]=it.second/new_o[it.first];
			}
			}
			
		}//end: for each transcript (column) in the current component
	}//end: for each component
		
				#if DEBUG_D
					cout<<"\nPrint d values / o"<<endl;
					print(d_value_new3);
					
					//exit(7);
				#endif
		
/*******************************************************************************************************
*	END	Step 5	*
********************************************************************************************************/	

	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~ Write Values to Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cout<<endl<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;
cout<<"Write Values to Files"<<endl;
cout<<"~~~~~~~~~~~~~~~~~~~~~~"<<endl;

/*******************************************************************************************************
*	O VALUES	- Are written in the same time with d values	*
********************************************************************************************************/		

	
/******************************************************************
*	D values and stats info *
******************************************************************/	
if (cc_option == 1)
{
//Prepare the structure for writing to file and compute the number of observed reads in each component:	
map<int, map<vector<string>, map<string, double> > > d_values2file;
//map<comp, map<RC, map<tr,value> > >
	
map<int, map<string, bool> > ntr;// number of transcripts for each component

//# of reads for each component
map<int, int> ncReads;
map<vector<string>,bool> ncReadsClasses;
		

map<int, map<string, map<vector<string>,double> > > *dMatrix;
if(least_square_type == 0)
	dMatrix = &d_value_new3;
else
	dMatrix = &d_values;
	
		#if DEBUG
			//This was done to check if the pointer works correct
			cout<<"dMatrix size = "<<(*dMatrix).size()<<endl;
			cout<<"d_value size = "<<d_values.size()<<endl;
			//Done
			
			//exit(7);
		#endif

//for each component
for (const auto &comp : (*dMatrix))
{
	cout<<"For each component"<<endl;
	//d_values (or d_value_new3) definition : map<int, map<string, map<vector<string>,double> > >
	
	ncReads[comp.first]=0;//initilize ncReads
	
	//for each transcript
	for (const auto &tr : comp.second) {
	
		//for each read class
		for (const auto &read_class : tr.second) {
		
			d_values2file[comp.first][read_class.first][tr.first]=read_class.second;
			
			//Collect the Read Classes names in this component
			//this loop covers all classes
			//even if some classes will be called several times (because they appear in different transcripts)
			//the values will be overwritten
			ncReadsClasses[read_class.first]=true;
		
				#if DEBUG
						cout<<"For each transcript in component"<<endl;
						cout<<"component = "<<comp.first<<endl;
						cout<<"read_class.second = "<<read_class.second<<endl;
				#endif
		
		
		}//end: for each read class
	ntr[comp.first][tr.first]=true;
	}//end: for each transcript
	
			#if DEBUG
				cout<<"Number of Reads Classes in this Component: "<<ncReadsClasses.size()<<endl;
					if(comp.first==1018)
						//exit (7);
			#endif
	//Sum up the total # of reads in this component
	//Note the definition for obsReadsClasses: map<vector<string>, int>
	//cout<<"For each read class in ncReadsClasses"<<endl;
	for(const auto &read_class : ncReadsClasses){
		ncReads[comp.first]+=obsReadsClasses[read_class.first];
		//cout<<"obsReadsClasses[read_class.first] = "<<obsReadsClasses[read_class.first]<<endl;
	}
	//cout<<"Number of reads in component "<<comp.first<<" is "<<ncReads[comp.first]<<endl;
	
	//reset ncReadsClasses (just to save space and reuse the same structure in the next component)
	ncReadsClasses.clear();
	
	
	/*Step2:Assign zero to all the other classes
		-check all read classes and if this transcript does not exist than assign zero to it
	*/
	//for each transcript
	for (const auto &tr : comp.second) {
		for (const auto &read_class : d_values2file[comp.first]){
			//cout<<"Test: "<<d_values2file[comp.first][read_class.first][tr.first]<<endl;
			
			if(d_values2file[comp.first][read_class.first][tr.first]==0)
				d_values2file[comp.first][read_class.first][tr.first]=0;
		}
		//it seems that if value doesn't exist, only because we've called it will assign 0 (zero) to it
	}
	
}//end for each component	



//Write d values to file
string d_file_new2="d_values.txt";
cout<<"\nWrite d_t,r values to: \n\t"<<d_file_new2<<endl;

ofstream d_stream_new2;
d_stream_new2.open(d_file_new2.c_str());

	if(!d_stream_new2){
		cout<<"Unable to open" <<d_file_new2<<endl;
		exit(1);
	}

//for each component
for (const auto &comp : d_values2file)
{//Definition for d_values2file: //map<comp, map<RC, map<tr,value> > >

	//Step0: Compute sum of observed frequencies
	//To DO: Improve running time for this sumation(maybe move it together with another loop)
	
	double sum_obs_freq=0.0; //Sum of observed reads classes frequencies per component
	//for each Read Class in current component
	for (const auto &read_class : comp.second)
		sum_obs_freq+=new_o[read_class.first];
	
	
	//Step1: Print stats
	d_stream_new2<<"Component: "<<comp.first<<" - #Reads: "<<ncReads[comp.first]<<" - #ReadClasses: "<<comp.second.size()<<" - SumObsFreq: "<<sum_obs_freq;
	d_stream_new2<<" - #transcripts: "<<ntr[comp.first].size();
	
	cout<<"Number of reads in component "<<comp.first<<": "<<ncReads[comp.first]<<endl;
	//exit(7);
	
	for (const auto &tr : ntr[comp.first])
		d_stream_new2<<" "<<tr.first;
	d_stream_new2<<endl;
	
	//Step2: Print read classes and values:
	//for each Read Class
	for (const auto &read_class : comp.second)
	{
		//Print first the corresponding observer read class frequency
		if(least_square_type == 0)
		{
			//0 - Relative -- i.e., we have already divided d/o -- just add 1
			d_stream_new2<<"1 ";
		}
		else
		{
			if(fabs(sum_obs_freq) < 1.e-20)
				d_stream_new2<<"0 ";
			else
				d_stream_new2<<new_o[read_class.first]/sum_obs_freq<<" ";
		}
		//print read class
		d_stream_new2<<"[";
		//for each tr in the read class
		for (const auto &tr_in_rc : read_class.first)
			d_stream_new2<<tr_in_rc<<"-";
		d_stream_new2<<"] ";
		
		//for each transcript
		for (const auto &tr : read_class.second)
		{
			d_stream_new2<<tr.second<<" ";
			
				#if DEBUG
					if(comp.first==1018)
					{
						cout<<"tr.second = "<<tr.second<<endl;
						
					}
				#endif
		}
	d_stream_new2<<endl;
	}//end: for each RC
	
d_stream_new2<<"---"<<endl;

	#if DEBUG
		if(comp.first==1018)
		{
			cout<<"Comp id = "<<comp.first<<endl;
			exit(7);
		}
	#endif
}//end: for each component

d_stream_new2.close();		

}
else
{//else no CC

//Step1: Prepare the structure for d_values without components
map<string, map<vector<string>,double> > d_values_noCC;

map<int, map<string, map<vector<string>,double> > > *dMatrix;
if(least_square_type == 0)
	dMatrix = &d_value_new3;
else
	dMatrix = &d_values;
	
//for each component
for (const auto &comp : (*dMatrix))
{
	//d_values (or d_value_new3) definition : map<int, map<string, map<vector<string>,double> > >

	//for each transcript
	for (const auto &tr : comp.second) {
	
		//for each read class
		for (const auto &read_class : tr.second) {
			d_values_noCC[tr.first][read_class.first]=read_class.second;
		}
	}
}

//Step2: Sort alphabeticatlly
//Maps are already sortded

	
//Step3: Write to file
/********************
*	O-values File	*
*********************/
string o_file_new="o_values_noCC.txt";
cout<<"\nWrite Observed Read Frequencies to: \n\t"<<o_file_new<<endl;
ofstream o_stream_new;
o_stream_new.open(o_file_new.c_str());
	if(!o_stream_new){
		cout<<"Unable to open" <<o_file_new<<endl;
		exit(1);
	}
/**********************/

	
/********************
*	D-values File	*
*********************/
string d_file_new2="d_values_noCC.txt";
cout<<"\nWrite d_t,r values to: \n\t"<<d_file_new2<<endl;

ofstream d_stream_new2;
d_stream_new2.open(d_file_new2.c_str());

	if(!d_stream_new2){
		cout<<"Unable to open" <<d_file_new2<<endl;
		exit(1);
	}
	
	//for each read class
	for (const auto &read_class : obsReadsClasses) {
	//read classes in ReadsClasses are also sorted
	
		//for each transcript in T (T contains all transcripts)
		for (const auto &tr : d_values_noCC){
			d_stream_new2<<d_values_noCC[tr.first][read_class.first]<<"\t";
			
			#if DEBUG
				//This is only for noCC option :)
				if (d_values_noCC[tr.first][read_class.first]<0)
				{
					cout<<"d_values_noCC[tr.first][read_class.first] = "<<d_values_noCC[tr.first][read_class.first]<<endl;
					exit(7);
				}
			#endif
			
		}
		d_stream_new2<<"\n";
		
		//Create also the o_values file
		if(least_square_type == 0)
		{
			//0 - Relative -- i.e., we have already divided d/o -- just add 1
			o_stream_new<<"1"<<endl;
		}
		else
			o_stream_new<<new_o[read_class.first]<<endl;
			
	}
	
o_stream_new.close();	
d_stream_new2.close();			
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
}//END: else no CC
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
		
cout<<"\nmain: Done"<<endl;
}//end main	