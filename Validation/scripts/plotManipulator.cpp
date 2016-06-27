/*
Author: Rumen Dangovski;
Date: 27.06.2016;
Build on code by Jack Wright.
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector> 
#include <fstream> 
#include <iomanip> 
#include <TCut.h>

using namespace std; 

const int topLIMIT = (1 << 20); 
const int bottomLIMIT = - (1 << 20); 

struct Plot { 
	int numBins; 
	int category; 

	bool isXLog; 
	bool isYLog; 
	
	string title; 
	string cuts;
	string var; 
	string xLabel; 
	string yLabel;
	//A formal expression of the histogram. Might be a useful feature.
	string expression; 
	
	float minX; 
	float maxX;

	void init () { 
		numBins = -1; 
		category = -1; 
		isXLog = false; 
		isYLog = false; 
		title = ""; 
		cuts = "";
	 	var = "";
		xLabel = ""; 
		yLabel = ""; 
		expression = ""; 
		minX = topLIMIT; 
		maxX = bottomLIMIT; 	
	}
	
};

string writeOn = "/afs/cern.ch/user/r/rumenrd/public/Plots/"; 

//The following function draws the fancy plots. 
void makePlot (TTree *tree, Plot &info) { 
	TCanvas c2 ("c2", "c2", 500, 500); 
	TCut cuts; 

	if (info.category >= 0 && !info.cuts.empty ()) 
		cuts = Form ("(%s) && (category==%d)", info.cuts.c_str(), info.category); 
	else if (info.category >= 0 && info.cuts.empty ())  
		cuts = Form ("(category==%d)", info.category); 
	else 
		cuts = ""; 

	tree->Draw (Form ("%s>>h(%d,%f,%f)", info.var.c_str(), info.numBins, info.minX, info.maxX), cuts); 
	TH1F *h = (TH1F*)gPad->GetPrimitive("h");

	TH1F *hist = new TH1F(*h);
	hist->SetName(Form("%s_hist",info.var.c_str()));
	hist->SetTitle(Form("%s Cat%d",info.title.c_str(),info.category));
	hist->GetXaxis()->SetTitle(info.xLabel.c_str());
	hist->GetYaxis()->SetTitle(info.yLabel.c_str());

	hist->SetMarkerColor(kBlack);
	hist->SetStats(kFALSE);
	hist->SetMarkerStyle(20);
	hist->SetLineWidth(2);
										    
	gStyle->SetEndErrorSize(3);
	hist->SetMarkerStyle(20);
												    
	hist->Draw("E1");
	if (info.isXLog) c2.SetLogx(); 
	if (info.isYLog) c2.SetLogy();

	c2.Print(Form("%s%s_cat%d_plot.pdf",writeOn.c_str(),info.var.c_str(),info.category));
	c2.Print(Form("%s%s_cat%d_plot.png",writeOn.c_str(),info.var.c_str(),info.category));
}	

void plotManipulator () { 
	cout << "Hello!" << endl;
}


