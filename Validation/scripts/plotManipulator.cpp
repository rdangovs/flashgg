/*
Author: Rumen Dangovski;
Date: 28.06.2016;
Inspired by code by Jack Wright.
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector> 
#include <fstream> 
#include <iomanip> 
#include <TCut.h>

using namespace std; 

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
		minX = 100000000; 
		maxX = -100000000; 	
	}
	
};

//predefines functions and a writing directory
string writeOn = "/afs/cern.ch/user/r/rdangovs/public/Plots2/"; 
vector <Plot> getPlotDetails (); 
void makePlot(TTree* tree, Plot &info);

vector <TH1F*> histograms;

void plotManipulator () {
	cout << "-----------------------------------------" << endl << "Wellcome to plotManipulator!" << endl << "-----------------------------------------" << endl; 

	gROOT->SetBatch(1); 
	TCanvas c1 ("c1", "c1", 500, 500); 

	//chooses the tree to work with
	TFile *_file0 = TFile::Open ("/afs/cern.ch/user/r/rdangovs/EXO_8_0_3_patch2/src/DataCollection/backgroundCut2.root"); 
	//TFile *_file0 = TFile::Open ("/afs/cern.ch/work/j/jwright/public/Louie/output_extras.root");
	TTree *tree = (TTree*)_file0->Get("flashggEXOValidationTreeMaker/diphotonTree_");
	
	gStyle->SetOptStat(11111);
	gStyle->SetOptFit(11111);

	vector <Plot> variables = getPlotDetails();
	
	for (int i = 0; i < variables.size (); ++ i) 
		makePlot (tree, variables[i]); 
	
	//std::cout << histograms.size () << std::endl;
	TH1F hans = (*histograms[0]) + (*histograms[1]); 
	/*
	string input; 
	while (true) { 
		cin >> input; 
		if (input == "end" || input == "END" || input == "End" || input == "E") break; 
	}*/
}


//The following function draws the fancy plots. 
void makePlot (TTree *tree, Plot &info) { 
	//need to fix for segmentation faults 
	
	//sets a canvas and an environment for the cuts 
	TCanvas c2 ("c2", "c2", 500, 500); 
	TCut cuts; 
	
	//clears things for the cuts
	if (info.category >= 0 && !info.cuts.empty ()) 
		cuts = Form ("(%s) && (category==%d)", info.cuts.c_str(), info.category); 
	else if (info.category >= 0 && info.cuts.empty ())  
		cuts = Form ("(category==%d)", info.category); 
	else 
		cuts = ""; 
	
	tree->Draw (Form ("%s>>h(%d,%f,%f)", info.var.c_str(), info.numBins, info.minX, info.maxX), cuts); 
	TH1F *h = (TH1F*)gPad->GetPrimitive("h");
	
	//starts a new histogram
	TH1F *hist = new TH1F(*h);
	hist->SetName(Form("%s_hist",info.var.c_str()));
	hist->SetTitle(Form("%s Cat%d",info.title.c_str(),info.category));
	hist->GetXaxis()->SetTitle(info.xLabel.c_str());
	hist->GetYaxis()->SetTitle(info.yLabel.c_str());

	hist->SetMarkerColor(kBlack);
	hist->SetStats(kFALSE);
	hist->SetMarkerStyle(27);
	hist->SetLineWidth(2);
										    
	gStyle->SetEndErrorSize(3);
	hist->SetMarkerStyle(20);
	
	//void SetHistFillColor(Color_t color = 30);
	hist->SetFillColor(30);
	hist->Draw("CONTE1bar");

	if (info.isXLog) c2.SetLogx(); 
	if (info.isYLog) c2.SetLogy();

	
	c2.Print(Form("%s%s_cat%d_plot.pdf",writeOn.c_str(), info.var.c_str(),info.category));
	c2.Print(Form("%s%s_cat%d_plot.png",writeOn.c_str(), info.var.c_str(),info.category));
	
	if ((info.var == "leadPt" || info.var == "subLeadPt") && (info.category == 0)) {
	//		std::cout << "YES" << " " << info.var.c_str() << std::endl;
			histograms.push_back (hist); 
	}

}	

vector<Plot> getPlotDetails(){

    vector<Plot> plotInfo;
    Plot info;
    
    //mgg
    info.var = "mgg";
    info.numBins = 69;
    info.minX = 230;
    info.maxX = 1610;
    info.category = 0;
    info.title = "m_{#gamma#gamma}";
    info.xLabel = "m_{#gamma#gamma} (GeV)";
    info.yLabel = "dN/dm_{#gamma#gamma}";
    info.isXLog = false;
    info.isYLog = true;
    plotInfo.push_back(info);

    info.category = 1;
    info.numBins = 63;
    info.minX = 320;
    info.maxX = 1580;
    plotInfo.push_back(info);

    //pT(lead)
    info.var = "leadPt";
    info.numBins = 34;
    info.minX = 75;
    info.maxX = 1200;
    info.category = 0;
    info.title = "Leading Photon p_{T}";
    info.xLabel = "p_T (GeV)";
    info.yLabel = "dN/dp_{T}";
    info.isXLog = false;
    info.isYLog = false;
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //pT(subLead);
    info.var = "subLeadPt";
    info.numBins = 34;
    info.minX = 75;
    info.maxX = 1200;
    info.category = 0;
    info.title = "Subleading Photon p_{T}";
    info.xLabel = "p_T (GeV)";
    info.yLabel = "dN/dp_{T}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead R9
    info.var = "leadR9";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 1;
    info.category = 0;
    info.title = "Leading Photon R_{9}";
    info.xLabel = "R_{9}";
    info.yLabel = "dN/dR_{9}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead R9
    info.var = "subLeadR9";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 1;
    info.category = 0;
    info.title = "Subleading Photon R_{9}";
    info.xLabel = "R_{9}";
    info.yLabel = "dN/dR_{9}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Supercluster eta
    info.var = "leadEtaSC";
    info.numBins = 34;
    info.minX = -2.5;
    info.maxX = 2.5;
    info.category = 0;
    info.title = "Leading Photon Supercluster #eta";
    info.xLabel = "#eta_{SC}";
    info.yLabel = "dN/#eta_{SC}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Supercluster eta
    info.var = "subLeadEtaSC";
    info.numBins = 34;
    info.minX = -2.5;
    info.maxX = 2.5;
    info.category = 0;
    info.title = "Subleading Photon Supercluster #eta";
    info.xLabel = "#eta_{SC}";
    info.yLabel = "dN/#eta_{SC}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Supercluster phi
    info.var = "leadPhiSC";
    info.numBins = 34;
    info.minX = -3.1415;
    info.maxX = 3.1415;
    info.category = 0;
    info.title = "Leading Photon Supercluster #phi";
    info.xLabel = "#phi_{SC}";
    info.yLabel = "dN/#phi_{SC}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Supercluster phi
    info.var = "subLeadPhiSC";
    info.numBins = 34;
    info.minX = -3.1415;
    info.maxX = 3.1415;
    info.category = 0;
    info.title = "Subleading Photon Supercluster #phi";
    info.xLabel = "#phi_{SC}";
    info.yLabel = "dN/#phi_{SC}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadChargedHadronIso";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 5;
    info.category = 0;
    info.title = "Leading Photon Charged Hadron Isolation";
    info.xLabel = "C.H.I.";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadChargedHadronIso";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 5;
    info.category = 0;
    info.title = "Subleading Photon Charged Hadron Isolation";
    info.xLabel = "CHI";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadPfPhoIso03";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 6;
    info.category = 0;
    info.title = "Leading Photon PF Photon Isolation 03";
    info.xLabel = "PFPI";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadPfPhoIso03";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 4;
    info.category = 0;
    info.title = "Subleading Photon PF Photon Isolation 03";
    info.xLabel = "PFPI";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadFull5x5_sigmaIetaIeta";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 0.04;
    info.category = 0;
    info.title = "Leading Photon #sigma_{i#eta i#eta}";
    info.xLabel = "#sigma_{i#eta i#eta}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadFull5x5_sigmaIetaIeta";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 0.04;
    info.category = 0;
    info.title = "Subleading Photon #sigma_{i#eta i#eta}";
    info.xLabel = "#sigma_{i#eta i#eta}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadFull5x5_r9";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 1.1;
    info.category = 0;
    info.title = "Leading Photon Full 5x5 R_{9}";
    info.xLabel = "R_{9}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadFull5x5_r9";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 1.1;
    info.category = 0;
    info.title = "Subleading Photon Full 5x5 R_{9}";
    info.xLabel = "R_{9}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadHadronicOverEm";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 0.06;
    info.category = 0;
    info.title = "Leading Photon Hadronic over EM";
    info.xLabel = "HoEM";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadHadronicOverEm";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 0.06;
    info.category = 0;
    info.title = "Subleading Photon Hadronic over EM";
    info.xLabel = "HoEM";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);
	
	/*
    //Jet Multiplicity pt > 20
    info.var = "jetMultiplicity_EGT20";
    info.numBins = 10;
    info.minX = 0;
    info.maxX = 10;
    info.category = 0;
    info.title = "Jet Multiplicity p_{T} > 20 GeV";
    info.xLabel = "N";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Jet Multiplicity pt > 30
    info.var = "jetMultiplicity_EGT30";
    info.numBins = 10;
    info.minX = 0;
    info.maxX = 10;
    info.category = 0;
    info.title = "Jet Multiplicity p_{T} > 30 GeV";
    info.xLabel = "N";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Jet Multiplicity pt > 40
    info.var = "jetMultiplicity_EGT40";
    info.numBins = 10;
    info.minX = 0;
    info.maxX = 10;
    info.category = 0;
    info.title = "Jet Multiplicity p_{T} > 40 GeV";
    info.xLabel = "N";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Leading Jet Pt
    info.cuts = "dijetLeadPt > 0";

    info.var = "dijetLeadPt";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 1000;
    info.category = 0;
    info.title = "Leading Jet p_{T}";
    info.xLabel = "p_{T} (GeV)";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Subleading Jet Pt
    info.var = "dijetSubleadPt";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 350;
    info.category = 0;
    info.title = "Subleading Jet p_{T}";
    info.xLabel = "p_{T} (GeV)";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Leading Jet Eta
    info.var = "dijetLeadEta";
    info.numBins = 34;
    info.minX = -5;
    info.maxX = 5;
    info.category = 0;
    info.title = "Leading Jet #eta";
    info.xLabel = "#eta";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Subleading Jet Eta
    info.var = "dijetSubleadEta";
    info.numBins = 34;
    info.minX = -5;
    info.maxX = 5;
    info.category = 0;
    info.title = "Leading Jet #eta";
    info.xLabel = "#eta";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Dijet Invariant Mass
    info.var = "dijetMass";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 2000;
    info.category = 0;
    info.title = "Dijet Invariant Mass";
    info.xLabel = "Inv. Mass (GeV)";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Dijet Delta Eta
    info.var = "dijetDeltaEta";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 8;
    info.category = 0;
    info.title = "Dijet Pseudorapidity Gap";
    info.xLabel = "#Delta#eta";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Dijet Zeppenfeld
    info.var = "dijetZeppenfeld";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 9;
    info.category = 0;
    info.title = "Dijet Zeppenfeld";
    info.xLabel = "#eta^{*}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Dijet Delta Phi jj
    info.var = "dijetDeltaPhi_jj";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 3.14159;
    info.category = 0;
    info.title = "Azimuthal Gap between Lead/Sublead Jets";
    info.xLabel = "#Delta#phi_{jj}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Dijet Delta Phi jj
    info.var = "dijetDeltaPhi_ggjj";
    info.numBins = 34;
    info.minX = 0;
    info.maxX = 3.14159;
    info.category = 0;
    info.title = "Azimuthal Gap between Diphoton and Dijet";
    info.xLabel = "#Delta#phi_{#gamma#gamma{jj}}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);
	*/

    return plotInfo;
}



