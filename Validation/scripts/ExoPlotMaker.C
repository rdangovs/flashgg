
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

struct PlotInfo {

    string var;
    unsigned nBins;
    float xMin;
    float xMax;
    int category;
    string title;
    string xLabel;
    string yLabel;

};

vector<PlotInfo> getPlotDetails();
void makeFancyPlot(TTree* tree, PlotInfo &info);
void makeSimplePlot(TTree* tree, string var, int cat);

void ExoPlotMaker(){

    TCanvas c1("c1");

    TFile *_file0 = TFile::Open("output.root");
    TTree *tree = (TTree*)_file0->Get("flashggEXOValidationTreeMaker/diphotonTree_");
    gStyle->SetOptStat(11111);
    gStyle->SetOptFit(11111);

    vector<PlotInfo> plotInfo = getPlotDetails();
    for (unsigned i=0;i<plotInfo.size();i++){
        makeFancyPlot(tree,plotInfo[i]);
    }

}

void makeFancyPlot(TTree* tree, PlotInfo &info){

    TCanvas c2("c2");
    TCut cut = Form("category==%d",info.category);
    tree->Draw(Form("%s>>h(%d,%f,%f)",info.var.c_str(),info.nBins,info.xMin,info.xMax),cut);
    TH1F *h = (TH1F*)gPad->GetPrimitive("h");

    TH1F *hist = new TH1F(*h);
    hist->SetName(Form("%s_hist",info.var.c_str()));
    hist->SetTitle(Form("%s Cat%d",info.title.c_str(),info.category));
    hist->GetXaxis()->SetTitle(info.xLabel.c_str());
    hist->GetYaxis()->SetTitle(info.yLabel.c_str());

    gStyle->SetOptStat(11111);
    gStyle->SetOptFit(11111);
    hist->SetMarkerColor(kBlack);
    hist->SetStats(kFALSE);

    hist->Draw("EP");
    c2.Print(Form("plots/%s_cat%d_plot.pdf",info.var.c_str(),info.category));

}

void makeSimplePlot(TTree* tree, string var, int cat){

    TCanvas c2("c2");
    tree->Draw(var.c_str(),Form("category==%d",cat));
    c2.Print(Form("plots/%s_cat%d_SimplePlot.pdf",var.c_str(),cat));

}

vector<PlotInfo> getPlotDetails(){

    vector<PlotInfo> plotInfo;
    PlotInfo info;
    //mgg
    info.var = "mgg";
    info.nBins = 34;
    info.xMin = 230;
    info.xMax = 910;
    info.category = 0;
    info.title = "m_{#gamma#gamma}";
    info.xLabel = "m_{#gamma#gamma} (GeV)";
    info.yLabel = "dN/dm_{#gamma#gamma}";
    plotInfo.push_back(info);

    info.category = 1;
    info.nBins = 25;
    info.xMin = 320;
    info.xMax = 910;
    plotInfo.push_back(info);

    //pT(lead)
    info.var = "leadPt";
    info.nBins = 34;
    info.xMin = 75;
    info.xMax = 1200;
    info.category = 0;
    info.title = "Leading Photon p_{T}";
    info.xLabel = "p_T (GeV)";
    info.yLabel = "dN/dp_{T}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //pT(subLead);
    info.var = "subLeadPt";
    info.nBins = 34;
    info.xMin = 75;
    info.xMax = 700;
    info.category = 0;
    info.title = "Subleading Photon p_{T}";
    info.xLabel = "p_T (GeV)";
    info.yLabel = "dN/dp_{T}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead R9
    info.var = "leadR9";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 1;
    info.category = 0;
    info.title = "Leading Photon R_{9}";
    info.xLabel = "R_{9}";
    info.yLabel = "dN/dR_{9}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead R9
    info.var = "subLeadR9";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 1;
    info.category = 0;
    info.title = "Subleading Photon R_{9}";
    info.xLabel = "R_{9}";
    info.yLabel = "dN/dR_{9}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Supercluster eta
    info.var = "leadEtaSC";
    info.nBins = 34;
    info.xMin = -2.5;
    info.xMax = 2.5;
    info.category = 0;
    info.title = "Leading Photon Supercluster #eta";
    info.xLabel = "#eta_{SC}";
    info.yLabel = "dN/#eta_{SC}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Supercluster eta
    info.var = "subLeadEtaSC";
    info.nBins = 34;
    info.xMin = -2.5;
    info.xMax = 2.5;
    info.category = 0;
    info.title = "Subleading Photon Supercluster #eta";
    info.xLabel = "#eta_{SC}";
    info.yLabel = "dN/#eta_{SC}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);






    //Lead Supercluster phi
    info.var = "leadPhiSC";
    info.nBins = 34;
    info.xMin = -3.1415;
    info.xMax = 3.1415;
    info.category = 0;
    info.title = "Leading Photon Supercluster #phi";
    info.xLabel = "#phi_{SC}";
    info.yLabel = "dN/#phi_{SC}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Supercluster phi
    info.var = "subLeadPhiSC";
    info.nBins = 34;
    info.xMin = -3.1415;
    info.xMax = 3.1415;
    info.category = 0;
    info.title = "Subleading Photon Supercluster #phi";
    info.xLabel = "#phi_{SC}";
    info.yLabel = "dN/#phi_{SC}";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadChargedHadronIso";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 400;
    info.category = 0;
    info.title = "Leading Photon Charged Hadron Isolation";
    info.xLabel = "C.H.I.";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadChargedHadronIso";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 300;
    info.category = 0;
    info.title = "Subleading Photon Charged Hadron Isolation";
    info.xLabel = "CHI";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadPfPhoIso03";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 150;
    info.category = 0;
    info.title = "Leading Photon PF Photon Isolation 03";
    info.xLabel = "PFPI";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadPfPhoIso03";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 150;
    info.category = 0;
    info.title = "Subleading Photon PF Photon Isolation 03";
    info.xLabel = "PFPI";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadFull5x5_sigmaIetaIeta";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 0.07;
    info.category = 0;
    info.title = "Leading Photon #sigma_{{i#eta}{i#eta}}";
    info.xLabel = "#sigma_{{i#eta}{i#eta}}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadFull5x5_sigmaIetaIeta";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 0.07;
    info.category = 0;
    info.title = "Subleading Photon #sigma_{{i#eta}{i#eta}}";
    info.xLabel = "#sigma_{{i#eta}{i#eta}}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadFull5x5_r9";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 1.3;
    info.category = 0;
    info.title = "Leading Photon Full 5x5 R_{9}";
    info.xLabel = "R_{9}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadFull5x5_r9";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 1.3;
    info.category = 0;
    info.title = "Subleading Photon Full 5x5 R_{9}";
    info.xLabel = "R_{9}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadHadronicOverEm";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 1.5;
    info.category = 0;
    info.title = "Leading Photon Hadronic over EM";
    info.xLabel = "HoEM";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadHadronicOverEm";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 1.5;
    info.category = 0;
    info.title = "Subleading Photon Hadronic over EM";
    info.xLabel = "HoEM";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    return plotInfo;
}


