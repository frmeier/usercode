#define noDOPIDTABLE

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
// root [2] .L efficiencyL0_01.C+
// root [3] efficiencyL0_01("pathToFile")


TCanvas *c, *c2;

void efficiencyL0_02(std::string fullPath, bool doOverview = true)
{
    const std::string outpdf("efficiencyL0_02plot");
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    if (doOverview)
	c = new TCanvas("c1","c1",1400,1200);
    else
	c = new TCanvas("c1","c1",600,600);
    const unsigned int nPadX = doOverview ? 2 : 1;
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
    //std::string cutacc = "genL0vtxR>1&&genL0vtxR<35&&TMath::Abs(genL0vtxZ)<100";
    std::string cutacc = "genL0pt>2&&genL0pt<4";
    //std::string cutacc = "";

    // Do plots
    int canvasCdCounter(0), canvasPageCounter(0);
    const int nbinseta(5);
    const int nbinspt(5);
    const double etamax(2.6);
    const double ptmax(40.);
    //const double rbins[] = {0,1,3,5,10,15,40};
    //const double rbins[] = {0,0.5,1,2,4,8,16,32,64};
    const double rbins[] = {0,.1,.2,.3,.4,0.5,1,2,4,8,16,64};
    const int rbins_size = (sizeof(rbins)/sizeof(double));
    std::vector<double> rbinvec(rbins,rbins+rbins_size);
    //const double zbins[] = {0.0,0.5,1.0,1.4,1.8,2.2};
    const double zbins[] = {0.0,10,20,40,60,80,100};
    const int zbins_size = (sizeof(zbins)/sizeof(double));
    std::vector<double> zbinvec(zbins,zbins+zbins_size);
    //const int nbinseta(2);
    //const double etamax(2.4);
    //const int nbinspt(1);
    const bool etabetrag(true);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doPlot2d((TPad*)c->cd(canvasCdCounter),tree,"h1", "genL0vtxR:TMath::Abs(genL0vtxZ)", cutacc.c_str(), zbinvec, rbinvec,"Generated #Lambda","|z(#Lambda)|","r(#Lambda)","cm","cm");
    setTH2params((TPad*)c->cd(canvasCdCounter),(TH2F*)gDirectory->GetList()->FindObject("h1"),true);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doPlot2d((TPad*)c->cd(canvasCdCounter),tree,"h2", "genL0vtxR:TMath::Abs(genL0vtxZ)", (cutacc+(cutacc.size()>0?"&&":"")+"L0matched==1").c_str(), zbinvec, rbinvec,"Generated #Lambda, found in reco","|z(#Lambda)|","r(#Lambda)","cm","cm");
    setTH2params((TPad*)c->cd(canvasCdCounter),(TH2F*)gDirectory->GetList()->FindObject("h2"),true);

#ifdef DOPIDTABLE
    // now we make a pidTable out of this
    PidTable *pid = new PidTable(1);
    pid->readFromHist(gDirectory, "h2", "h1");
    pid->dumpToFile("pid_lambda0.dat");
#endif

    // make a clone for the ratio
    {
	TH2F *hratio = (TH2F*)gDirectory->GetList()->FindObject("h2")->Clone("h2Dratio");
    }
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doPlotRatio2d((TPad*)c->cd(canvasCdCounter),"#Lambda reco efficiency","h1","h2Dratio");

    // Now we do a 1D histo for each x-bin
    //plot1Dfrom2DforeachXbin((TH2F*)gDirectory->GetList()->FindObject("h1"));
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    plot1Dfrom2DforeachXbin((TPad*)c->cd(canvasCdCounter),"h2Dratio");

    // 1D ratio plots
    const double rbinsratio[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,5,6,8,10,40,120};
    const int rbinsratio_size = (sizeof(rbinsratio)/sizeof(double));
    std::vector<double> rbinratiovec(rbinsratio,rbinsratio+rbinsratio_size);
    const double zbinsratio[] = {0.0,10,20,40,60,80,100};
    const int zbinsratio_size = (sizeof(zbinsratio)/sizeof(double));
    std::vector<double> zbinratiovec(zbinsratio,zbinsratio+zbinsratio_size);

    std::string cutaccZ = (cutacc.size()>0?cutacc+"&&":"")+ "TMath::Abs(genL0vtxZ)<100";
    std::string cutaccR = (cutacc.size()>0?cutacc+"&&":"")+ "genL0vtxR>1&&genL0vtxR<35";

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tree,tree,"hvtxR","genL0vtxR","genL0vtxR",cutaccZ.c_str(), (cutaccZ+(cutaccZ.size()>0?"&&":"")+"L0matched==1").c_str(),25,0,40,"#Lambda reco efficiency for |z|<100","r","cm");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tree,tree,"hvtxRzoom","genL0vtxR","genL0vtxR",cutaccZ.c_str(), (cutaccZ+(cutaccZ.size()>0?"&&":"")+"L0matched==1").c_str(),25,0,8,"#Lambda reco efficiency for |z|<100","r","cm");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tree,tree,"hvtxZ","TMath::Abs(genL0vtxZ)","TMath::Abs(genL0vtxZ)",cutaccR.c_str(), (cutaccR+(cutaccR.size()>0?"&&":"")+"L0matched==1").c_str(),25,0,100,"#Lambda reco efficiency for 2<r<40","|z|","cm");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tree,tree,"hvtxZzoom","TMath::Abs(genL0vtxZ)","TMath::Abs(genL0vtxZ)",cutaccR.c_str(), (cutaccR+(cutaccR.size()>0?"&&":"")+"L0matched==1").c_str(),25,0,20,"#Lambda reco efficiency for 2<r<40","|z|","cm");

    // finalize current page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
}

