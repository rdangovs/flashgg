#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector> 
#include <fstream> 

using namespace std; 

struct Plot { 
	int numBins; 
	int category; 

	bool isXLog; 
	bool isYLog; 
	
	string title; 
	string cuts;
	string param; 
	string xLabel; 
	string yLabel; 
	
	float minX; 
	float maxX;
};

string writeOn = "/afs/cern.ch/user/r/rumenrd/public/Plots/"; 


