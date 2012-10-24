#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "setTDRStyle_modified.C"
#include "TLatex.h"

#include "Cuts.C"
#include "utils.h"

void quickPvLipPlot(string filename, string addCut, bool isB0 = false, string filenamePDF = "")
{
    setTDRStyle();
    TFile *f = TFile::Open(filename.c_str());
    TTree *events = (TTree*)f->Get("events");

    const int nBins = 100;
    const double lo = -.1;
    const double hi = .1;
    string plotstring = "(100,-.1,.1)";
    Cuts c;
    if (!isB0) c.selectCut("acc06Lb", "muSoft", "lb14", "HLT_jpsiBarrel", "HLT_matched");
    else c.selectCut("acc06B0", "muSoft", "B008", "HLT_jpsiBarrel", "HLT_matched");
    string fullCut = addCut.size()==0 ? c.getCut() : c.getCut() + "&&" + addCut;
    events->Draw(makePlotsString("PvLip","h1",nBins,lo,hi).c_str(), fullCut.c_str());
    events->Draw(makePlotsString("PvLip2","h2",nBins,lo,hi).c_str(), fullCut.c_str(),"same");

    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("h1");
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject("h2");

    h1->SetTitle("");
    h1->SetLineColor(2);
    h1->SetFillColor(2);
    h1->SetFillStyle(3005);
    h1->GetXaxis()->SetNdivisions(505);
    h1->GetXaxis()->SetTitle(valueWithUnit("longitudinal i.p.","cm").c_str());
    h1->GetYaxis()->SetTitle(entriesPerBin(nBins, lo, hi,"cm").c_str());

    h2->SetTitle("");
    h2->SetLineColor(4);
    h2->SetFillColor(4);
    h2->SetFillStyle(3004);
    h2->GetXaxis()->SetNdivisions(505);

    gPad->SetLogy();

    TLegend *leg = new TLegend(.15, .80, .50, .90);
    leg->AddEntry(h1, "lip to closest PV", "f");
    leg->AddEntry(h2, "lip to second closest PV", "f");
    leg->Draw();
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    TLatex* particlelabel = writeTLatex(isB0 ? "B^{0}" : "#Lambda_{b}", .82, .82);
    particlelabel->Draw();

    if (filenamePDF.size() != 0) gPad->SaveAs(filenamePDF.c_str());
}

void mkSomePvLipPlots()
{
    setTDRStyle();
    quickPvLipPlot("../data/run479.root", "", true, "r479_lip_cuts.pdf");
    quickPvLipPlot("../data/run479.root", "(mbc<5.25||mbc>5.308)", true, "r479_lip_cuts_sigwindow.pdf");
    quickPvLipPlot("../data/run548.root", "", false, "r548_lip_cuts.pdf");
    quickPvLipPlot("../data/run548.root", "(mbc<5.582||mbc>5.654)", false, "r548_lip_cuts_sigwindow.pdf");
}

