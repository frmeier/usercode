#include "TROOT.h"
#include "TPad.h"
#include "TTree.h"
#include "TH1F.h"
#include "TPad.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "utils.h"
#include <string>

void quickEffPlots(string filename)
{
    TFile *f = TFile::Open(filename.c_str());
    //TTree *events = (TTree*)f->Get("events");
    TTree *genevents = (TTree*)f->Get("genevents");

    const int nBins(100);
    const double lo_p(0), hi_p(20), lo_pi(0), hi_pi(10);
    const string drawstrP = "("+toString(nBins)+","+toString(lo_p)+","+toString(hi_p)+")";
    const string drawstrPi = "("+toString(nBins)+","+toString(lo_pi)+","+toString(hi_pi)+")";
    genevents->Draw(("phap>>h1_Lb_p"+drawstrP).c_str(), "qha1>0");
    gPad->Update();
    genevents->Draw(("phap>>h2_Lb_p"+drawstrP).c_str(), "qha1>0&&hasCand==1");
    gPad->Update();
    TH1F *h1_Lb_p = (TH1F*)gDirectory->Get("h1_Lb_p");
    TH1F *h2_Lb_p = (TH1F*)gDirectory->Get("h2_Lb_p");
    //TEfficiency eff_Lb_p(*h2, *h1);

    genevents->Draw(("pham>>h1_Lb_pi"+drawstrPi).c_str(), "qha1>0");
    gPad->Update();
    genevents->Draw(("pham>>h2_Lb_pi"+drawstrPi).c_str(), "qha1>0&&hasCand==1");
    gPad->Update();
    TH1F *h1_Lb_pi = (TH1F*)gDirectory->Get("h1_Lb_pi");
    TH1F *h2_Lb_pi = (TH1F*)gDirectory->Get("h2_Lb_pi");
    //TEfficiency eff_Lb_pi(*h2, *h1);

    genevents->Draw(("pham>>h1_Lbbar_p"+drawstrP).c_str(), "qha1<0");
    gPad->Update();
    genevents->Draw(("pham>>h2_Lbbar_p"+drawstrP).c_str(), "qha1<0&&hasCand==1");
    gPad->Update();
    TH1F *h1_Lbbar_p = (TH1F*)gDirectory->Get("h1_Lbbar_p");
    TH1F *h2_Lbbar_p = (TH1F*)gDirectory->Get("h2_Lbbar_p");
    //TEfficiency eff_Lbbar_p(*h2, *h1);

    genevents->Draw(("phap>>h1_Lbbar_pi"+drawstrPi).c_str(), "qha1<0");
    gPad->Update();
    genevents->Draw(("phap>>h2_Lbbar_pi"+drawstrPi).c_str(), "qha1<0&&hasCand==1");
    gPad->Update();
    TH1F *h1_Lbbar_pi = (TH1F*)gDirectory->Get("h1_Lbbar_pi");
    TH1F *h2_Lbbar_pi = (TH1F*)gDirectory->Get("h2_Lbbar_pi");
    //TEfficiency eff_Lbbar_pi(*h2, *h1);

    //eff_Lb_p.Draw("");
    //eff_Lb_p.Draw("same");
    //eff_Lb_p.GetPaintedHistogram()->SetLineColor(2);
    //eff_Lbbar_p.GetPaintedHistogram()->SetLineColor(4);
}

void quickEffPlotsB0(string filename)
{
    TFile *f = TFile::Open(filename.c_str());
    //TTree *events = (TTree*)f->Get("events");
    TTree *genevents = (TTree*)f->Get("genevents");

    const int nBins(100);
    const double lo_p(0), hi_p(20), lo_pi(0), hi_pi(10);
    const string drawstrP = "("+toString(nBins)+","+toString(lo_p)+","+toString(hi_p)+")";
    const string drawstrPi = "("+toString(nBins)+","+toString(lo_pi)+","+toString(hi_pi)+")";
    genevents->Draw(("phap>>h1_B0_pip"+drawstrPi).c_str(), "");
    gPad->Update();
    genevents->Draw(("phap>>h2_B0_pip"+drawstrPi).c_str(), "hasCand==1");
    gPad->Update();
    TH1F *h1_B0_p = (TH1F*)gDirectory->Get("h1_B0_pip");
    TH1F *h2_B0_p = (TH1F*)gDirectory->Get("h2_B0_pip");
    //TEfficiency eff_B0_p(*h2, *h1);

    genevents->Draw(("pham>>h1_B0_pim"+drawstrPi).c_str(), "");
    gPad->Update();
    genevents->Draw(("pham>>h2_B0_pim"+drawstrPi).c_str(), "hasCand==1");
    gPad->Update();
    TH1F *h1_B0_pi = (TH1F*)gDirectory->Get("h1_B0_pim");
    TH1F *h2_B0_pi = (TH1F*)gDirectory->Get("h2_B0_pim");
    //TEfficiency eff_B0_pi(*h2, *h1);

    //eff_B0_p.Draw("");
    //eff_B0_p.Draw("same");
    //eff_B0_p.GetPaintedHistogram()->SetLineColor(2);
    //eff_B0bar_p.GetPaintedHistogram()->SetLineColor(4);
}

