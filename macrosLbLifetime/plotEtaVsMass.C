
#include <string>
#include <iostream>

#include "TROOT.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TLegend.h"

#include "setTDRStyle_modified.C"

using std::string;
using std::cout;
using std::endl;

void plotEtaVsMeass(string filename, bool isB0)
{
    setTDRStyle();

    TFile *_file0 = TFile::Open(filename.c_str());

    TTree *tree = (TTree*)gDirectory->Get("events");

    TCanvas *c = new TCanvas("c","c",1000,800);
    c->Divide(2,2);

    c->cd(1);
    tree->Draw("mB0:TMath::Abs(etaB0)>>hprof(24,0,2.4,20,5.18,5.38)","isMCmatch==1","prof");

    TProfile *hprof = (TProfile*)gDirectory->GetList()->FindObject("hprof");
    if (isB0) hprof->GetXaxis()->SetTitle("#eta(B^{0})");
    else hprof->GetXaxis()->SetTitle("#eta(#Lambda_{b})");
    if (isB0) hprof->GetYaxis()->SetTitle("m(B^{0})");
    else hprof->GetXaxis()->SetTitle("m(#Lambda_{b})");

    TH1D *hrms = hprof->ProjectionX("hrms","C=E");
    hrms->GetYaxis()->SetTitle("RMS");

    c->cd(2);
    hrms->Draw();

    c->cd(3);
    tree->Draw("TMath::Abs(etaB0)>>heta(24,0,2.4)","isMCmatch==1");

    c->cd(4);
    tree->Draw("mB0>>hmB01(20,5.18,5.38)","isMCmatch==1&&TMath::Abs(etaB0)>=0.0&&TMath::Abs(etaB0)<0.2");
    tree->Draw("mB0>>hmB02(20,5.18,5.38)","isMCmatch==1&&TMath::Abs(etaB0)>=1.0&&TMath::Abs(etaB0)<1.2","same");
    tree->Draw("mB0>>hmB03(20,5.18,5.38)","isMCmatch==1&&TMath::Abs(etaB0)>=2.0&&TMath::Abs(etaB0)<2.2","same");
    ((TH1D*)gDirectory->GetList()->FindObject("hmB01"))->SetLineColor(1);
    ((TH1D*)gDirectory->GetList()->FindObject("hmB02"))->SetLineColor(2);
    ((TH1D*)gDirectory->GetList()->FindObject("hmB03"))->SetLineColor(3);
    gPad->Update();

    TLegend *leg = new TLegend(0.25,0.75,0.48,0.85);
    leg->AddEntry((TH1D*)gDirectory->GetList()->FindObject("hmB01"), "0.0 #leq |#eta| < 0.2", "lp");
    leg->AddEntry((TH1D*)gDirectory->GetList()->FindObject("hmB02"), "1.0 #leq |#eta| < 1.2", "lp");
    leg->AddEntry((TH1D*)gDirectory->GetList()->FindObject("hmB03"), "2.0 #leq |#eta| < 2.2", "lp");
    leg->Draw();
}

