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

////////////////////////////////////////////////////////////
// makes efficiency plots vy r and z of L0 decay vertex
// in slices of pt

// Usage
// root [0] gSystem->Load("lib/libUtil.so");
// root [1] gSystem->Load("lib/libAnaClasses.so");
// root [2] .L efficiencyL0_01.C+
// root [3] efficiencyL0_01("pathToFile")


TCanvas *c, *c2;

std::string makePtCut(std::string cutname, double lo, double hi)
{
    return cutname + ">" + toString(lo) + "&&" + cutname + "<=" + toString(hi);
}

std::string makePtCutTitle(double lo, double hi)
{
    return toString(lo) + " < p_{T} #leq " + toString(hi) + " GeV/c";
}

void efficiencyL0_03(std::string fullPath, bool doOverview = true)
{
    const std::string outpdf("efficiencyL0_03plot");
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    if (doOverview)
	c = new TCanvas("c1","c1",1400,1200);
    else
	c = new TCanvas("c1","c1",600,600);
    const unsigned int nPadX = doOverview ? 4 : 1;
    const unsigned int nPadY = doOverview ? 4 : 1;
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
    //std::string cutacc = "genL0pt>2&&genL0pt<4";
    std::string cutacc = "";

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

    // Now we do a 1D histo for each x-bin
    //plot1Dfrom2DforeachXbin((TH2F*)gDirectory->GetList()->FindObject("h1"));
    //canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    //plot1Dfrom2DforeachXbin((TPad*)c->cd(canvasCdCounter),"h2Dratio");

    // 1D ratio plots
    /*
    const double rbinsratio[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,5,6,8,10,40,120};
    const int rbinsratio_size = (sizeof(rbinsratio)/sizeof(double));
    std::vector<double> rbinratiovec(rbinsratio,rbinsratio+rbinsratio_size);
    const double zbinsratio[] = {0.0,10,20,40,60,80,100};
    const int zbinsratio_size = (sizeof(zbinsratio)/sizeof(double));
    std::vector<double> zbinratiovec(zbinsratio,zbinsratio+zbinsratio_size);
    */
    //const double pTsclices[] = {0,1,2,3,4,6,50};
    const double pTsclices[] = {2,50};
    const int pTsclices_size = (sizeof(pTsclices)/sizeof(double));
    std::vector<double> pTsclicesVec(pTsclices,pTsclices+pTsclices_size);

    std::string cutaccZ = (cutacc.size()>0?cutacc+"&&":"")+ "TMath::Abs(genL0vtxZ)<100";
    std::string cutaccR = (cutacc.size()>0?cutacc+"&&":"")+ "genL0vtxR>1&&genL0vtxR<35";

    double maxvalR(0), maxvalZ(0);;

    for(int i = 0; i < pTsclices_size-1; i++)
    {
	const std::string cutPt = makePtCut("genL0pt", pTsclices[i], pTsclices[i+1]);
	const std::string cutPtTitle = makePtCutTitle(pTsclices[i], pTsclices[i+1]);
	cout << i << ": Doing plot for " << cutPtTitle <<  endl;
	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tree,tree,"hvtxR"+toString(i),"genL0vtxR","genL0vtxR",(cutaccZ+"&&"+cutPt).c_str(), (cutaccZ+(cutaccZ.size()>0?"&&":"")+"L0matched==1&&"+cutPt).c_str(),25,0,40,"#Lambda eff, |z|<100, "+cutPtTitle,"r","cm");
	const double curMaxValR = ((TH1F*)gDirectory->GetList()->FindObject(("hvtxR"+toString(i)+"2").c_str()))->GetMaximum();
	if (curMaxValR > maxvalR) maxvalR = curMaxValR;

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tree,tree,"hvtxZ"+toString(i),"TMath::Abs(genL0vtxZ)","TMath::Abs(genL0vtxZ)",cutaccR+"&&"+cutPt, (cutaccR+(cutaccR.size()>0?"&&":"")+"L0matched==1&&"+cutPt).c_str(),25,0,100,"#Lambda eff, 2<r<40, "+cutPtTitle,"|z|","cm");
	const double curMaxValZ = ((TH1F*)gDirectory->GetList()->FindObject(("hvtxZ"+toString(i)+"2").c_str()))->GetMaximum();
	if (curMaxValZ > maxvalZ) maxvalZ = curMaxValZ;
    }

    // overlay the efficiencies
    const int colarr[] = {1,2,3,4,6,7,8,9};
    const int colarr_size = (sizeof(colarr)/sizeof(int));
    { // R
	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	(TPad*)c->cd(canvasCdCounter);
	TLegend *leg = new TLegend(0.4,0.6,0.89,0.89);
	for(int i = 0; i < pTsclices_size-1; i++)
	{
	    TH1F* h = (TH1F*)((TH1F*)gDirectory->GetList()->FindObject(("hvtxR"+toString(i)+"2").c_str()))->Clone(("hvtxR"+toString(i)).c_str());
	    h->SetMaximum(maxvalR);
	    h->SetMarkerColor(colarr[i%colarr_size]);
	    h->GetYaxis()->SetTitleOffset(1.5);
	    h->GetYaxis()->SetTitle("Efficiency");
	    h->GetXaxis()->SetTitleOffset(1.5);
	    h->SetLineColor(colarr[i%colarr_size]);
	    h->DrawCopy(i==0?"":"same");
	    leg->AddEntry(h,(makePtCutTitle(pTsclices[i], pTsclices[i+1])).c_str());
	}
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->Draw();
    }
    { // Z
	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	TPad *pad = (TPad*)c->cd(canvasCdCounter);
	TLegend *leg = new TLegend(0.4,0.6,0.89,0.89);
	for(int i = 0; i < pTsclices_size-1; i++)
	{
	    TH1F* h = (TH1F*)((TH1F*)gDirectory->GetList()->FindObject(("hvtxZ"+toString(i)+"2").c_str()))->Clone(("hvtxZ"+toString(i)).c_str());
	    h->SetMaximum(maxvalR);
	    h->SetMarkerColor(colarr[i%colarr_size]);
	    h->SetLineColor(colarr[i%colarr_size]);
	    h->GetYaxis()->SetTitleOffset(1.5);
	    h->GetYaxis()->SetTitle("Efficiency");
	    h->GetXaxis()->SetTitleOffset(1.5);
	    h->Draw(i==0?"":"same");
	    leg->AddEntry(h,(makePtCutTitle(pTsclices[i], pTsclices[i+1])).c_str());
	}
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->Draw();
    }

    // finalize current page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
}

