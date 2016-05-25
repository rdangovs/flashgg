
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

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
    bool logXPlot;
    bool logYPlot;

};

string outputDir = "/afs/cern.ch/user/j/jwright/www/EXO_Plots_17_05/";

vector<PlotInfo> getPlotDetails();
void makeFancyPlot(TTree* tree, PlotInfo &info);
void makeSimplePlot(TTree* tree, string var, int cat);
void compareWithEXO(TTree *tree, PlotInfo &info);

void ExoPlotMaker(){

    gROOT->SetBatch(1);
    TCanvas c1("c1","c1",500,500);

    TFile *_file0 = TFile::Open("/afs/cern.ch/work/j/jwright/public/Louie/output_17_05.root");
    TTree *tree = (TTree*)_file0->Get("flashggEXOValidationTreeMaker/diphotonTree_");
    gStyle->SetOptStat(11111);
    gStyle->SetOptFit(11111);

    vector<PlotInfo> plotInfo = getPlotDetails();
    for (unsigned i=0;i<plotInfo.size();i++){
        makeFancyPlot(tree,plotInfo[i]);
    }
    compareWithEXO(tree,plotInfo[0]);
    compareWithEXO(tree,plotInfo[1]);

}

void compareWithEXO(TTree *tree, PlotInfo &info){

    string cat;
    if (info.category == 0){
        cat = "EBEB";
    }else{
        cat = "EBEE";
    }

    ofstream file;
    string filePath = Form("%sComparison%s.txt",outputDir.c_str(),cat.c_str());
    file.open(filePath.c_str());

    file << "Comparing " << info.var << " plots in category " << cat << endl;
    if (file.is_open()) {
        cout << "Successfully opened " << filePath << endl;
    }else{
        cout << "Failed to open " << filePath << endl;
    }

    TFile *_file0 = TFile::Open("/afs/cern.ch/user/m/musella/public/workspace/exo/full_analysis_spring15_7415v2_sync_v5_data_ecorr_cic2_final_ws.root");
    TTree *t_EBEB=(TTree*)_file0->Get(Form("tree_data_cic2_%s",cat.c_str()));
    t_EBEB->Draw(Form("%s>>h(%d,%f,%f)",info.var.c_str(),info.nBins,info.xMin,info.xMax));
    TH1F *h = (TH1F*)gPad->GetPrimitive("h");

    TCut cut = Form("category==%d",info.category);
    tree->Draw(Form("%s>>histNew(%d,%f,%f)",info.var.c_str(),info.nBins,info.xMin,info.xMax),cut);
    TH1F *histNew = (TH1F*)gPad->GetPrimitive("histNew");

    file << "Pasquale's has " << h->Integral() << " events, new has " << histNew->Integral() << endl;

    file << setw(24) << "Pasquale" << setw(24) << "Us" << setw(24) << "Compare" << endl;
    file << setw(12) << "Centre" << setw(12) << "Height";
    file << setw(12) << "Centre" << setw(12) << "Height";
    file << setw(12) << "Difference" << setw(12) << "\% diff";
    file << endl;
    for(unsigned i=1;i<=h->GetNbinsX();i++){
        file << setw(12) << h->GetBinCenter(i);
        file << setw(12) << h->GetBinContent(i);
        file << setw(12) << histNew->GetBinCenter(i);
        file << setw(12) << histNew->GetBinContent(i);
        file << setw(12) << histNew->GetBinContent(i) - h->GetBinContent(i);
        file << setw(12) << 100*(histNew->GetBinContent(i) - h->GetBinContent(i))/h->GetBinContent(i);
        file << endl;
    }
    TCanvas c2("c2");
    if (info.logXPlot){
        c2.SetLogx();
    }
    if (info.logYPlot){
        c2.SetLogy();
    }
    h->SetLineColor(kBlue);
    histNew->SetLineColor(kRed);
    h->Draw();
    histNew->Draw("same");
    c2.Print(Form("%s%s_cat%d_comparisonplot.pdf",outputDir.c_str(),info.var.c_str(),info.category));
    c2.Print(Form("%s%s_cat%d_comparisonplot.png",outputDir.c_str(),info.var.c_str(),info.category));

    c2.Clear();
    h->SetMarkerColor(kBlack);
    h->SetStats(kFALSE);
    h->SetMarkerStyle(20);
    h->SetLineWidth(2);
    h->Draw("E1");
    c2.Print(Form("%s%s_cat%d_Pasquale.pdf",outputDir.c_str(),info.var.c_str(),info.category));
    c2.Print(Form("%s%s_cat%d_Pasquale.png",outputDir.c_str(),info.var.c_str(),info.category));

    file.close();
}












void makeFancyPlot(TTree* tree, PlotInfo &info){

    TCanvas c2("c2","c2",500,500);
    TCut cut = Form("category==%d",info.category);
    tree->Draw(Form("%s>>h(%d,%f,%f)",info.var.c_str(),info.nBins,info.xMin,info.xMax),cut);
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
    if (info.logXPlot){
        c2.SetLogx();
    }
    if (info.logYPlot){
        c2.SetLogy();
    }

    c2.Print(Form("%s%s_cat%d_plot.pdf",outputDir.c_str(),info.var.c_str(),info.category));
    c2.Print(Form("%s%s_cat%d_plot.png",outputDir.c_str(),info.var.c_str(),info.category));
}

vector<PlotInfo> getPlotDetails(){

    vector<PlotInfo> plotInfo;
    PlotInfo info;
    //mgg
    info.var = "mgg";
    info.nBins = 69;
    info.xMin = 230;
    info.xMax = 1610;
    info.category = 0;
    info.title = "m_{#gamma#gamma}";
    info.xLabel = "m_{#gamma#gamma} (GeV)";
    info.yLabel = "dN/dm_{#gamma#gamma}";
    info.logXPlot = false;
    info.logYPlot = true;
    plotInfo.push_back(info);

    info.category = 1;
    info.nBins = 63;
    info.xMin = 320;
    info.xMax = 1580;
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
    info.logXPlot = false;
    info.logYPlot = false;
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
    info.xMax = 5;
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
    info.xMax = 5;
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
    info.xMax = 6;
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
    info.xMax = 4;
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
    info.xMax = 0.04;
    info.category = 0;
    info.title = "Leading Photon #sigma_{i#eta i#eta}";
    info.xLabel = "#sigma_{i#eta i#eta}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Sublead Charged Hadron Isolation
    info.var = "subLeadFull5x5_sigmaIetaIeta";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 0.04;
    info.category = 0;
    info.title = "Subleading Photon #sigma_{i#eta i#eta}";
    info.xLabel = "#sigma_{i#eta i#eta}";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    //Lead Charged Hadron Isolation
    info.var = "leadFull5x5_r9";
    info.nBins = 34;
    info.xMin = 0;
    info.xMax = 1.1;
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
    info.xMax = 1.1;
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
    info.xMax = 0.06;
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
    info.xMax = 0.06;
    info.category = 0;
    info.title = "Subleading Photon Hadronic over EM";
    info.xLabel = "HoEM";
    info.yLabel = "";
    plotInfo.push_back(info);

    info.category = 1;
    plotInfo.push_back(info);

    return plotInfo;
}


