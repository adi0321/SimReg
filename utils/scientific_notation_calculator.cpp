//Reads in a a number (or string), and converts it to an float

/*Compile command:
g++ scientific_notation_calculator.cpp -o scientific_notation_calculator
*/

#include <cstdlib>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace std;

#define HELPMESSAGE "too few parameters: \n\
Input:\n\
argv[1] - gene frequency \n\
argv[2] - number of transcripts (for this gene) \n\
argv[3] - 0 for Random - 1 for uniform \n\
"

int main(int argc,char *argv[])
{
	
	if(argc<4){
		cout<<HELPMESSAGE<<endl;
		exit(1);
	}

float gene_freq=atof(argv[1]);
int nt = atoi(argv[2]); // number of transcripts
int option = atoi(argv[3]); // 0-Random, 1-Uniform

float t[nt];

	if(option==0)
	{
	/* initialize random seed: */
	srand (time(NULL));

		//assign random number and compute total sum
		float sum=0.0;
		for(int i=0;i<nt;i++)
		{
		t[i]=rand();
		sum+=t[i];
		}

		//normalize and print on the screen
		for(int i=0;i<nt;i++)
		{
			t[i]=(t[i]/sum)*gene_freq;
			cout<<t[i]<<endl;
		}
	}
	else
	{
		float tr_freq = gene_freq/nt;
		for(int i=0;i<nt;i++)
			cout<<tr_freq;
	}
}
