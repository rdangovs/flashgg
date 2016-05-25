#!/usr/bin/env python
# Simple script to make Effiency X Acceptance plot from Binned Baseline/Massfac analysis
# run with python makeEffAcc.py CMS-HGG.root
import ROOT as r
import sys
import re
import string
import random

r.gSystem.Load("Signal/lib/libSimultaneousSignalFit.so")
r.gSystem.Load("libHiggsAnalysisCombinedLimit")
r.gSystem.Load("libdiphotonsUtils")
GetBR = lambda x : norm.GetBR(float(x))
GetXsection = lambda x : norm.GetXsection(float(x))
GetProcXsection = lambda x,y : norm.GetXsection(x,y)

r.gROOT.SetBatch(1)
infileName1=sys.argv[1]
infileName2=sys.argv[2]
f1 = r.TFile(infileName1)
f2 = r.TFile(infileName2)
ws = f1.Get("wtemplates")
mggEBEB = wtemplates.var("mggEBEB")
a = wtemplates.data("data_EBEB")

TList *list = new TList; 
list->Add(t1); 
list->Add(t2); 
t1 = f2Get("cic/trees/"); 
TTree *newtree = TTree::MergeTrees(list); 
newtree->SetName("T"); 
newtree->Write(); 

#outws = r.RooWorkspace("multipdf_reduced");
#data = ws.data("roohist_data_mass_UntaggedTag_2")
#dataReduced = data.emptyClone(data.GetName()+"_reduced",data.GetName()+"_reduced")
#random = r.TRandom3()
#random.Rndm()
#rArray=[]
#random.RndmArray(data.numEntries(),rArray)
#for i in range(0,data.numEntries()):
#  print "bin ", i, " data.get(i) " , data.get(i) , " data.weight() ", data.weight() 
#  dataReduced.add(data.get(i),int(data.weight()*0.1))
#  #if rArray[i] > 0.1:
#    #dataReduced.add(data.get(i),data.weight())
#print "Original Dataset"
#data.Print()
#print "Reduced  Dataset"
#dataReduced.Print()
##outws.r.RooWorkspace.import(dataReduced)
#getattr(outws,'import')(dataReduced)
#outfile.cd()
#outws.Write()
#outfile.Close()
#exit(1)
#print " NOW MULTIDIMFIT WS"
#ws.getSnapshot("MultiDimFit").Print("V")
#print " NOW CLEAN WS"

ws.Print("V")
#ws.obj("shapeBkg_bkg_ch1_EBEB").Print("V")
#print " current pdf index "

#print ws.obj("shapeBkg_bkg_ch1_EBEB").getCurrentIndex()
#print " NOW NORMAL WS"

#ws.Print("V")

exit(1)

