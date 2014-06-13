#include "norm.h"
#include <string>
#include <map>

//normalize vector
void norm(vector<float> &a){
	float sum=0;
	for(int i=0;i<a.size();i++)
		sum+=a[i];
	for(int i=0;i<a.size();i++)
		a[i]/=sum;
}
void norm(vector<double> &a){
	double sum=0;
	for(int i=0;i<a.size();i++)
		sum+=a[i];
	for(int i=0;i<a.size();i++)
		a[i]/=sum;
}

void norm(map<string, double> &a){
	double sum=0.0;
	
	for (const auto &it : a)
		sum+=it.second;
	
	for (const auto &it : a)
		a[it.first] = it.second / sum;
}

void norm(map<string, map<string, double> > &a){
        double sum=0.0;

        for (const auto &g : a)
		for (const auto &t : g.second)
                	sum+=t.second;

        for (const auto &g : a)
		for (const auto &t : g.second)
                	a[g.first][t.first] = t.second / sum;
}

