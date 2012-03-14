#define DOPIDTABLE

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "TEventList.h"
#include "setTDRStyle_modified.C"
#include "cuts.C"
#include "efficiencyPlots.h"
#include "canvaspager.h"

#ifdef DOPIDTABLE
#include "/shome/meier_f1/CMSSW/CMSSW_3_8_6/src/AnalysisDataFormats/HeavyFlavorObjects/rootio/PidTable.hh"
#endif

// Usage
// root [0] gSystem->Load("lib/libUtil.so");
// root [1] gSystem->Load("lib/libAnaClasses.so");
// root [2] .L efficiencyCuts_01.C+
// root [3] efficiencyL0_01("pathToFile")


TCanvas *c, *c2;

void efficiencyCuts_01(std::string fullPath, bool doOverview = true)
{
    const std::string outpdf("efficiencyCuts_01plot");
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    if (doOverview)
	c = new TCanvas("c1","c1",1400,1200);
    else
	c = new TCanvas("c1","c1",600,600);
    const unsigned int nPadX = doOverview ? 3 : 1;
    const unsigned int nPadY = doOverview ? 2 : 1;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTrees
    TTree* tree = (TTree*) f->Get("events");
    if (tree==0)
    {
	cout << "Tree events not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0) cout << "Got TTree evemts with " << tree->GetEntries() << " entries" << endl;

    // General cuts
    Cuts cutsAcc, cutsAll;
    //cuts.selectCut("an03","acc03","HLT_matched_01");
    cutsAcc.selectCut("acc03","HLT_matched_01");
    cutsAll.selectCut("an05","acc03","HLT_matched_01");
    std::string addcut = "&&isSig==1&&isMCmatch==1";

    // As the cutsAcc are included in cutsAll, we can speed up by making a TEventList
    tree->Draw(">>lst",(cutsAcc.getCut()+addcut).c_str());
    TEventList* lst;
    lst = (TEventList*)gDirectory->Get("lst");
    tree->SetEventList(lst);

    // Do plots
    int canvasCdCounter(0), canvasPageCounter(0);
    const int nbinseta(5);
    const int nbinspt(5);
    const double etamax(2.6);
    const double ptmax(40.);
    //const double ptbins[] = {0,1,3,5,10,15,40};
    const double ptbins[] = {0,4.5,6.0,8,11,18,50};
    const int ptbins_size = (sizeof(ptbins)/sizeof(double));
    std::vector<double> ptbinvec(ptbins,ptbins+ptbins_size);
    //const double etabins[] = {0.0,0.5,1.0,1.4,1.8,2.2};
    const double etabins[] = {0.0,1.8,2.5};
    const int etabins_size = (sizeof(etabins)/sizeof(double));
    std::vector<double> etabinvec(etabins,etabins+etabins_size);
    //const int nbinseta(2);
    //const double etamax(2.4);
    //const int nbinspt(1);
    const bool etabetrag(true);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doPlot2d((TPad*)c->cd(canvasCdCounter),tree,"h1", "ptlb:TMath::Abs(etalb)", (cutsAcc.getCut()+addcut).c_str(), etabinvec, ptbinvec,"#Lambda_{b} before cuts","|#eta(#Lambda_{b})|","p_{T}(#Lambda){b}","","GeV/c");
    setTH2params((TPad*)c->cd(canvasCdCounter),(TH2F*)gDirectory->GetList()->FindObject("h1"),true);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doPlot2d((TPad*)c->cd(canvasCdCounter),tree,"h2", "ptlb:TMath::Abs(etalb)", (cutsAll.getCut()+addcut).c_str(), etabinvec, ptbinvec,"#Lambda_{b} after cuts","|#eta(#Lambda_{b})|","p_{T}(#Lambda_{b})","","GeV/c");
    setTH2params((TPad*)c->cd(canvasCdCounter),(TH2F*)gDirectory->GetList()->FindObject("h2"),true);

#ifdef DOPIDTABLE
    // now we make a pidTable out of this
    PidTable *pid = new PidTable(1);
    pid->readFromHist(gDirectory, "h2", "h1");
    pid->dumpToFile("pid_cuts.dat");
#endif

    // make a clone for the ratio
    {
	TH2F *hratio = (TH2F*)gDirectory->GetList()->FindObject("h2")->Clone("h2Dratio");
    }
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doPlotRatio2d((TPad*)c->cd(canvasCdCounter),"#Lambda_{b} cut efficiency","h1","h2Dratio");

    // 1D ratio plots
    const double ptbinsratio[] = {0,4.5,6.5,8,10,12,14,16,18,22,40};
    const int ptbinsratio_size = (sizeof(ptbinsratio)/sizeof(double));
    std::vector<double> ptbinratiovec(ptbinsratio,ptbinsratio+ptbinsratio_size);
    const double etabinsratio[] = {0,0.4,0.8,1.2,1.5,1.7,1.9,2.1,2.25,2.4,2.5,2.6};
    const int etabinsratio_size = (sizeof(etabinsratio)/sizeof(double));
    std::vector<double> etabinratiovec(etabinsratio,etabinsratio+etabinsratio_size);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tree,tree,"hpt","ptlb","ptlb",(cutsAcc.getCut()+addcut).c_str(), (cutsAll.getCut()+addcut).c_str(),ptbinratiovec,"#Lambda_{b} cut efficiency for |#eta|<2.4","p_{T}","GeV/c");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tree,tree,"heta","TMath::Abs(etalb)","TMath::Abs(etalb)",(cutsAcc.getCut()+addcut).c_str(), (cutsAll.getCut()+addcut).c_str(),etabinratiovec,"#Lambda_{b} cut efficiency for 3<p_{T}<40","|#eta|","");

    // finalize current page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
}

