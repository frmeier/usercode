#include <string>
#include <iostream>

#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"

#include "setTDRStyle_modified.C"

using std::string;
using std::cout;
using std::endl;

void quickLifetimeFit(string filename, bool isB0)
{
    setTDRStyle();

    TFile *_file0 = TFile::Open(filename.c_str());

    TTree *tree = (TTree*)gDirectory->Get("fittree");

    TCanvas *c = new TCanvas("c","c",600,800);
    c->Divide(2,3);
    int canvasctr(0);

    string strSigWindow, strBgrWindow;
    double SoverSB;

    if (isB0)
    {
	strSigWindow = "mass>5.25&&mass<5.31";
	strBgrWindow = "mass>5.35&&mass<5.55";
	SoverSB = 0.90; // from another fit
    }
    else
    {
	strSigWindow = "mass>5.60&&mass<5.64";
	strBgrWindow = "mass>5.75&&mass<6.30";
	SoverSB = 0.60; // from another fit
    }

    c->cd(++canvasctr);
    tree->Draw("mass>>hmsigwindow", strSigWindow.c_str());
    TH1F *hmsigwindow = (TH1F*)gDirectory->GetList()->FindObject("hmsigwindow");
    const int nSigWindow = hmsigwindow->GetEntries();

    c->cd(++canvasctr);
    tree->Draw("mass>>hbgrwindow", strBgrWindow.c_str());
    TH1F *hbgrwindow = (TH1F*)gDirectory->GetList()->FindObject("hbgrwindow");
    const int nBgrWindow = hbgrwindow->GetEntries();

    const double nSig = nSigWindow * SoverSB;
    const double nBgr = nSigWindow * (1.-SoverSB);
    const double scale = (double)nBgr / (double)nBgrWindow;

    c->cd(++canvasctr);
    tree->Draw("t>>hsigt(34,-2e-12,15e-12)", strSigWindow.c_str());
    TH1F *hsigt = (TH1F*)gDirectory->GetList()->FindObject("hsigt");
    gPad->SetLogy();

    c->cd(++canvasctr);
    tree->Draw("t>>hbgrt(34,-2e-12,15e-12)", strBgrWindow.c_str());
    TH1F *hbgrt = (TH1F*)gDirectory->GetList()->FindObject("hbgrt");
    hbgrt->Scale(scale);
    gPad->SetLogy();

    c->cd(++canvasctr);
    hsigt->Draw();
    hsigt->SetLineColor(4);
    hbgrt->Draw("same");
    hbgrt->SetLineColor(2);
    gPad->SetLogy();

    c->cd(++canvasctr);
    TH1F *hdiff = new TH1F("hdiff","hdiff",34,-2e-12,15e-12);
    hdiff->Add(hsigt,hbgrt,1,-1);
    hdiff->Draw();
    gPad->SetLogy();

    c->SaveAs("quickFit.pdf");

    TCanvas *c2 = new TCanvas("c2","c2",1000,800);
    hdiff->Draw();
    gPad->SetLogy();

}

/*
root [2] fittree->Draw("mass","mass>5.60&&mass<5.64")
(Long64_t)944
root [3] fittree->Draw("mass","mass>5.75&&mass<6.3")
(Long64_t)3681
root [4] fittree->Draw("t","mass>5.60&&mass<5.64")
    (Long64_t)944
    root [5] fittree->Draw("t>>h1(20,-2e12,8e12)","mass>5.60&&mass<5.64")
    (Long64_t)944
    root [6] fittree->Draw("t>>h1(20,-2e-12,8e-12)","mass>5.60&&mass<5.64")
    (Long64_t)944
    root [7] fittree->Draw("t>>h1(30,-2e-12,13e-12)","mass>5.60&&mass<5.64")
    (Long64_t)944
    root [8] fittree->Draw("t>>h2(30,-2e-12,13e-12)","mass>5.75&&mass<6.30","same")
*/
