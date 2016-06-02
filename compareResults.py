#!/usr/bin/env python
# Simple script to make Effiency X Acceptance plot from Binned Baseline/Massfac analysis
# run with python makeEffAcc.py CMS-HGG.root
import ROOT as r
import sys
import re
import string
import random
import math 

r.gSystem.Load("Signal/lib/libSimultaneousSignalFit.so")
r.gSystem.Load("libHiggsAnalysisCombinedLimit")
r.gSystem.Load("libdiphotonsUtils")
GetBR = lambda x : norm.GetBR(float(x))
GetXsection = lambda x : norm.GetXsection(float(x))
GetProcXsection = lambda x,y : norm.GetXsection(x,y)

r.gROOT.SetBatch(1)
rdsList=[]
th1List=[]
mggEBEB =None
print "arguments " , sys.argv[1:] 
for infileName in sys.argv[1:]: 
  print " processing file ", infileName.split("/")[-1]
  #infileName1=sys.argv[1]
  #infileName2=sys.argv[2]
  #f1 = r.TFile(infileName1)
  #f2 = r.TFile(infileName2)
  f = r.TFile(infileName)
  ws = f.Get("wtemplates")
  #print "WS ", ws, " hasattr(ws,'var')", hasattr(ws,'var')
  hasWS=False
  if (hasattr(ws,'var')): 
    #print "setting hasWS to true"
    hasWS=True
  if (hasWS):
    #print "using the hasWS worksflow"
    mggEBEB = ws.var("mggEBEB")
    a = ws.data("data_EBEB")
    rdsList.append(a)
    for i in range(0,a.numEntries()):
      m = a.get(i).getRealValue("mggEBEB")
      if (m > 1300):
      #if ((abs(m - 236)<10)):
        print "For file:", infileName.split("/")[-1], "  (hasWS)  m = ", m
  else:
    t = f.Get("cic/trees/Data_13TeV_EBEB");
    if (t==None) : t = f.Get("cic/trees/test_13TeV_EBEB");
    if (t==None) : t = f.Get("flashggEXOValidationTreeMaker/diphotonTree_") 
    print "test1 t.GetEntries() ",t.GetEntries()
    for i in range(0,t.GetEntries()):
      t.GetEntry(i)
      m=-999
      massName=""
      if hasattr(t,'mass'):
        m = t.mass
        massName='mass'
      if hasattr(t,'mgg'):
        m = t.mgg
        massName='mgg'
      if (m >1300): 
      #if ((abs(m - 236)<10)):
        print "For file:", infileName.split("/")[-1], " (!hasWS)  m = ", m
    print "%s>>htemp(69,230,1610)"%massName
    t.Draw("%s>>htemp(69,230,1610)"%massName,"category==0")
    htemp = r.gROOT.FindObject("htemp")
    print " htemp", htemp, htemp.GetEntries()
    th1List.append(htemp.Clone())

c1 = r.TCanvas("c1","c1",500,500)
mggEBEB.setMin(230)
mggEBEB.setMax(1610)
mggEBEB.setBins(69)

frame = mggEBEB.frame()
print "f ",frame, " mggEBEB ", mggEBEB
for a in rdsList:
  a.plotOn(frame)
for b in th1List:
  print " htemp", b, b.GetEntries()
  c = r.RooDataHist ("c","c",r.RooArgList(mggEBEB),b)
  #b.Draw("same")
  c.plotOn(frame,r.RooFit.MarkerColor(r.kRed),r.RooFit.LineColor(r.kRed))

c1.SetLogy()
frame.Draw()
c1.SaveAs("test.pdf")



exit(1)

