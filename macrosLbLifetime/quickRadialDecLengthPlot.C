#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"

#include "setTDRStyle_modified.C"
#include "Cuts.C"
#include "utils.h"

void doQuickRadialDecLengthPlot(string filename, double declength, int nBins = 100)
{
    setTDRStyle();

    TFile *f = TFile::Open(filename.c_str());
    TTree *events = (TTree*)f->Get("events");
    TTree *genevents = (TTree*)f->Get("genevents");

    const double lo = 0;
    const double hi = declength;
    genevents->Draw(makePlotsString("vrbcPV","h1",nBins,lo,hi).c_str(),"","");
    genevents->Draw(makePlotsString("vrrsPV","h2",nBins,lo,hi).c_str(),"","same");

    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("h1");
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject("h2");

    h1->SetLineColor(2);
    h2->SetLineColor(4);
    h1->GetXaxis()->SetTitle(valueWithUnit("radial decay length","cm").c_str());
    h1->GetYaxis()->SetTitle(entriesPerBin(nBins,lo,hi,"cm").c_str());
    h2->GetXaxis()->SetTitle(valueWithUnit("radial decay length","cm").c_str());
    h2->GetYaxis()->SetTitle(entriesPerBin(nBins,lo,hi,"cm").c_str());
    h1->SetTitle("");
    h2->SetTitle("");
    gPad->SetLogy(true);

    TLegend *leg = new TLegend(.70, .80, .87, .90);
    leg->AddEntry(h1, "r(B^{0})", "lm");
    leg->AddEntry(h2, "r(V^{0})", "lm");
    leg->Draw();
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
}

void drawFinalPlot()
{
    doQuickRadialDecLengthPlot("../data/run546.root", 4,200);
    TPad *npad = new TPad("npad","",0.6, 0.25, 0.9, 0.55);
    npad->cd();
    doQuickRadialDecLengthPlot("../data/run546.root", 40,2000);
    npad->Draw();
}

