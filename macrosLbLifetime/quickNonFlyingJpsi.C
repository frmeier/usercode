
#include <string>
#include <iostream>

#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TChain.h"

#include "utils.h"
#include "setTDRStyle_modified.C"

using std::string;
using std::cout;
using std::endl;

void quickNonFlyingJpsi()
{
    setTDRStyle();

    //TFile *f1 = TFile::Open("run323.root");
    //TFile *f2 = TFile::Open("run324.root");
    //TFile *f3 = TFile::Open("run325.root");

    TChain *t1 = new TChain("events"); t1->Add("../data/run323.root");
    TChain *t2 = new TChain("events"); t2->Add("../data/run324.root");
    TChain *t3 = new TChain("events"); t3->Add("../data/run325.root");

    string cut = "mlb>5.52&&mlb<5.72";
    string toDraw = "ct3dlb";
    string binning = "(100,-10e-12,30e-12)";
    t1->Draw((toDraw+">>h1"+binning).c_str(), cut.c_str());
    t2->Draw((toDraw+">>h2"+binning).c_str(), cut.c_str(), "same");
    t3->Draw((toDraw+">>h3"+binning).c_str(), cut.c_str(), "same");

    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("h1");
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject("h2");
    TH1F *h3 = (TH1F*)gDirectory->GetList()->FindObject("h3");

    h1->SetFillColor(0);
    h2->SetFillColor(0);
    h3->SetFillColor(0);

    h1->SetLineColor(1);
    h2->SetLineColor(2);
    h3->SetLineColor(3);
}

