//Compile command:
//g++ simreg.cpp -fopenmp -o simreg 

#include <stdlib.h>     /* system, NULL, EXIT_FAILURE */
#include <cmath>

#include <string>
#include <iostream>
#include <algorithm>
#include "tclap/CmdLine.h"

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>

#define DEBUG 0

using namespace TCLAP;
using namespace std;

int getdir (string dir, vector<string> &directory)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
		//for some reason it also returns the . and .. (which we are not interested ) 
        if ((string(dirp->d_name).compare(".") != 0) && (string(dirp->d_name).compare("..") != 0)){
			directory.push_back(string(dirp->d_name));
		}
		//cout<<string(dirp->d_name)<<endl;
		//exit(7);
    }
    closedir(dp);
    return 0;
}

int main(int argc, char** argv)
{

cout<<"\nRunning "<<argv[0]<<endl;

    	// Wrap everything in a try block.  Do this every time, 
	// because exceptions will be thrown for problems. 
	try {  

	// Define the command line object.
	CmdLine cmd("Command description message", ' ', "0.9");

	// Define a value argument and add it to the command line.
	ValueArg<string> gtfarg("G","GTF","Path to GTF File",true,"homer","string");
	cmd.add( gtfarg );
	
	ValueArg<string> faarg("F","FA","Path to FA File",true,"homer","string");
	cmd.add( faarg );
	
	ValueArg<string> ccpatharg("C","CC_PATH","Path to Connected Components",true,"homer","string");
	cmd.add( ccpatharg );
	
	ValueArg<string> dirpatharg("s","SCRIPT_DIR","Path to Directory with all scripts",true,"homer","string");
	cmd.add( dirpatharg );
	
	ValueArg<string> meanarg("m","mean","-m <int> Fragment length mean",true,"homer","string");
	cmd.add( meanarg );
	
	ValueArg<string> devarg("d","deviation","-d <int> Fragment length standard deviation",true,"homer","string");
	cmd.add( devarg );
	
	ValueArg<string> precisarg("t","tuning","-t <int> Tuning(Precision)",true,"homer","string");
	cmd.add( precisarg );

	// Define a switch and add it to the command line.
	SwitchArg simSwitch("g","grinder","Use Grinder for MonteCarlo reads", false);
	cmd.add( simSwitch );

	// Parse the args.
	cmd.parse( argc, argv );

	// Get the value parsed by each arg. 
	string gtf=gtfarg.getValue();
	string fa=faarg.getValue();
	string CC_PATH=ccpatharg.getValue();
	string SCRIPT_DIR=dirpatharg.getValue();
	string tuning=precisarg.getValue();
	
	bool grinder = simSwitch.getValue();

	string mean=meanarg.getValue();
	string dev=devarg.getValue();
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Do what you intend too...


			#if DEBUG
				cout<<"GTF File = "<<gtf<<endl;
				cout<<"FA File = "<<fa<<endl;
				cout<<"CC Path = "<<CC_PATH<<endl;
				cout<<"SCRIPT_DIR = "<<SCRIPT_DIR<<endl;
				cout<<"precision = "<<tuning<<endl;
				//exit(7);
			#endif
	
	//Run this for loop for each component
	//So get the components names from the path
	
    vector<string> directory = vector<string>();

    getdir(CC_PATH,directory);

		#if DEBUG
			/*for (unsigned int i = 0;i < directory.size(); i++) {
				cout << directory[i] << endl;
			}
			*/
			//exit(7);
		#endif
    

	#pragma omp parallel for num_threads(30) schedule(dynamic)
	for(int i=0; i<directory.size(); i++)
	{
	
		//int i=0;
		//directory.push_back("1");
		
		//stringstream cpath_command;
		//cpath_command<<CC_PATH<<"/"<<i<<"*";
			
			#if DEBUG
				cout<<"directory["<<i<<"] = "<<directory[i]<<endl;
				//cout<<"cpath_command = "<<cpath_command.str().c_str()<<endl;
				//exit(7);
			#endif
			
		
		//string cpath;
		//cpath=system(cpath_command.str().c_str());
		//cout<<cpath<<endl;
		
		//exit(7);
		
		//for each component run simreg.sh
		//${SCRIPT_DIR}/scripts/simReg.sh -m $mean -d $deviation -l $read_length -r $rpf -G $GTF_File -F $FA_File -s $SCRIPT_DIR -C $CC_Path
		
		stringstream bash_command;
		
		if(grinder)
			bash_command<<"nice "<<SCRIPT_DIR<<"/scripts/simReg_cc.sh -g -m "<<mean<<" -d "<<dev<<" -l 100 -r 2 -G "<<gtf<<" -F "<<fa<<" -s "<<SCRIPT_DIR<<" -C "<<CC_PATH<<" -c "<<directory[i];
		else
			bash_command<<"nice "<<SCRIPT_DIR<<"/scripts/simReg_cc.sh -m "<<mean<<" -d "<<dev<<" -t "<<tuning<<" -l 100 -r 2 -G "<<gtf<<" -F "<<fa<<" -s "<<SCRIPT_DIR<<" -C "<<CC_PATH<<" -c "<<0;//directory[i];
		
			#if DEBUG
				cout<<"bash_command = "<<bash_command.str().c_str()<<endl;
				//exit(7);
			#endif
		
		system(bash_command.str().c_str());
		
		//exit(7);
	}
		
		
	} catch (ArgException &e)  // catch any exceptions
	{ cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }

	

}



