#include <iostream>
#include <string>
#include <sstream>
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

TCanvas *c, *c2;

void doPlot1d(TTree *t, std::string hname, std::string todraw, std::string cut,
	int nBins, double min, double max, std::string title, std::string titleX, std::string unitX)
{
    const std::string plotstring = todraw + ">>" + hname + "(" +
	toString(nBins) + "," + toString(min) + "," + toString(max) + ")";
    t->Draw(plotstring.c_str(), cut.c_str());
    TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(hname.c_str());
    h->SetTitle(title.c_str());
    h->SetStats(0);
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" ["+unitX+"]").c_str() : titleX.c_str());
    const double binsize = (max-min)/nBins;
    h->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    return; 
}

void doPlot1d(TTree *t, std::string hname, std::string todraw, 
	int nBins, double min, double max, std::string title, std::string titleX, std::string unitX)
{
    doPlot1d(t,hname,todraw,"",nBins,min,max,title,titleX,unitX);
}

void repositionStatbox(std::string name)
{
    // Reposition and resize statbox
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(name.c_str());
    TPaveStats *st = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.63);
    st->SetY1NDC(0.72);
    st->SetX2NDC(0.99);
    st->SetY2NDC(0.99);
    st->Draw();
}

void canvaspager(TCanvas *c, std::string filename, const int &ncd, int &cdcounter, int &pagecounter)
{
    cdcounter++;
    if (cdcounter > ncd)
    {
	c->SaveAs((filename + toString(pagecounter) + ".pdf").c_str());
	c->SaveAs((filename + toString(pagecounter) + ".png").c_str());
	pagecounter++;
	cdcounter=1;
	c->Clear("D");
    }
}

void dataPlotsJpL0(std::string fullPath)
{
    const std::string outpdf("dataPlotsJpL0plots");
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c1","c1",1000,600);
    int canvasCdCounter(0), canvasPageCounter(0);
    const unsigned int nPadX = 1;
    const unsigned int nPadY = 1;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    for(unsigned int i=1; i<=nPads; i++)
    {
	TPad* pad= (TPad*)c->cd(i);
	pad->SetTopMargin(0.10);
    }
    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTree
    TTree* t = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;

    // General cuts
    Cuts cuts;
    //cuts.selectCut("an04","acc03","HLT_matched_01");
    cuts.selectCut("an05exp","acc03","HLT_matched_2011_01");
    std::string cutsAnalysis = cuts.getCut();

    // Do plots
    const double mjp(3.097), mjpwindow(0.3);
    const double mjplo(mjp-mjpwindow), mjphi(mjp+mjpwindow);
    const double ml0(1.116), ml0window(0.030);
    const double ml0lo(ml0-ml0window), ml0hi(ml0+ml0window);

    const int nBins(20);

    // create an additional cut for the mlb window
    Cuts mlbcut;
    mlbcut.selectCut("mlbWindow01");

    // Apply cuts
    //t->Draw(">>lst",cutsgen.c_str());
    //TEventList* lst;
    //lst = (TEventList*)gDirectory->Get("lst");
    //t->SetEventList(lst);

    TFile *fout = new TFile("histosdataPlotsJpL0.root","recreate");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hmjp", "mjp", (cuts.cs.getCutExceptOne(cuts.parvec,cuts.cs.getCutPos("mjp"))+"&&"+mlbcut.getCut()).c_str(),nBins,mjplo,mjphi,"","m(#mu#mu)","GeV/c^{2}");
    c->Update();

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hml0", "ml0", (cuts.cs.getCutExceptOne(cuts.parvec,cuts.cs.getCutPos("ml0"))+"&&"+mlbcut.getCut()).c_str(),nBins,ml0lo,ml0hi,"","m(p#pi)","GeV/c^{2}");
    c->Update();

    /*
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hmjp", "mjp", cuts.cs.getCutExceptOne(cuts.parvec,cuts.cs.getCutPos("mjp")),nBins,mjplo,mjphi,"","m(#mu#mu)","GeV/c^{2}");
    c->Update();

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hml0", "ml0", cuts.cs.getCutExceptOne(cuts.parvec,cuts.cs.getCutPos("ml0")),nBins,ml0lo,ml0hi,"","m(p#pi)","GeV/c^{2}");
    c->Update();
    */

    // finalize current page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".png").c_str());

    fout->Write();
}

