#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TPad.h"

#include "Cuts.C"
#include "utils.h"

using std::cout;
using std::endl;
using std::string;

void doQuickEff2d(TTree *tree, string name, string drawX, int nBinsX, double minX, double maxX, string drawY, int nBinsY, double minY, double maxY, string cutAll, string cutPass)
{
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    c->Divide(3,2);

    const string drawstring = drawY + ":" + drawX;
    const string binstring = "("+toString(nBinsX)+","+toString(minX)+","+toString(maxX)+","+toString(nBinsY)+","+toString(minY)+","+toString(maxY)+")";

    c->cd(1);
    tree->Draw((drawstring+">>hall"+binstring).c_str(), cutAll.c_str(), "colz");

    c->cd(2);
    tree->Draw((drawstring+">>hpass"+binstring).c_str(), cutPass.c_str(), "colz");

    TH2F* hpass = (TH2F*)gDirectory->GetList()->FindObject("hpass");
    TH2F* hall = (TH2F*)gDirectory->GetList()->FindObject("hall");
    TEfficiency *eff = new TEfficiency(*hpass, *hall);

    c->cd(3);
    eff->Draw("surf3");
    gPad->Update();
    eff->GetPaintedHistogram()->GetXaxis()->SetTitle(drawX.c_str());
    eff->GetPaintedHistogram()->GetXaxis()->SetTitleOffset(1.8);
    eff->GetPaintedHistogram()->GetYaxis()->SetTitle(drawY.c_str());
    eff->GetPaintedHistogram()->GetYaxis()->SetTitleOffset(1.8);
    eff->GetPaintedHistogram()->Sumw2();

    c->cd(4);
    TH1D *effX = eff->GetPaintedHistogram()->ProjectionX(("heff_X"+name).c_str(), 0, -1, "e");
    effX->Scale(1./(double)nBinsX);
    effX->SetTitle(name.c_str());
    effX->Draw();
    c->cd(5);
    TH1D *effY = eff->GetPaintedHistogram()->ProjectionY(("heff_Y"+name).c_str(), 0, -1, "e");
    effY->Scale(1./(double)nBinsY);
    effY->SetTitle(name.c_str());
    effY->Draw();
    effY->Print("all");

    c->SaveAs((name+".pdf").c_str());
}

void doSomePlots(string filename)
{
    TFile *f = TFile::Open(filename.c_str());
    TTree *events = (TTree*)f->Get("events");
    TTree *genevents = (TTree*)f->Get("genevents");

    // Muon acceptance
    const string acc_mu1_1 = "TMath::Abs(etamu1)<1.2&&ptmu1>3.5";
    const string acc_mu1_2 = "TMath::Abs(etamu1)>=1.2&&TMath::Abs(etamu1)<1.6&&ptmu1>8-3.75*TMath::Abs(etamu1)";
    const string acc_mu1_3 = "TMath::Abs(etamu1)>=1.6&&TMath::Abs(etamu1)<2.4&&ptmu1>2.0";
    const string acc_mu2_1 = "TMath::Abs(etamu2)<1.2&&ptmu2>3.5";
    const string acc_mu2_2 = "TMath::Abs(etamu2)>=1.2&&TMath::Abs(etamu2)<1.6&&ptmu2>8-3.75*TMath::Abs(etamu2)";
    const string acc_mu2_3 = "TMath::Abs(etamu2)>=1.6&&TMath::Abs(etamu2)<2.4&&ptmu2>2.0";
    const string acc_mu = "(("+acc_mu1_1+")||("+acc_mu1_2+")||("+acc_mu1_3+"))&&(("+acc_mu2_1+")||("+acc_mu2_2+")||("+acc_mu2_3+"))";

    doQuickEff2d(genevents, "eff2d_muAcc", "d3dlb", 20, 0, .3, "plb", 20, 0, 30, "", acc_mu);
    doQuickEff2d(genevents, "eff2d_muAcc_wide", "d3dlb", 30, 0, 1, "plb", 30, 0, 100, "", acc_mu);

    doQuickEff2d(genevents, "eff2d_muAcc_cand", "d3dlb", 20, 0, .3, "plb", 20, 0, 30, "", acc_mu+"&&hasCand==1");
    doQuickEff2d(genevents, "eff2d_muAcc_cand_wide", "d3dlb", 30, 0, 1, "plb", 30, 0, 100, "", acc_mu+"&&hasCand==1");

    doQuickEff2d(genevents, "eff2d_cand", "d3dlb", 20, 0, .3, "plb", 20, 0, 30, "", "hasCand==1");
    doQuickEff2d(genevents, "eff2d_cand_wide", "d3dlb", 30, 0, 1, "plb", 30, 0, 100, "", "hasCand==1");


    Cuts cutLb;
    cutLb.selectCut("lb11", "acc05Lb", "muSoft");
    const string trigPassBarrel = "HLTmatch==1&&HLTokBarrelJpsi==1";
    const string trigPassDispl = "HLTmatch==1&&HLTokDisplJpsi==1";

    doQuickEff2d(events, "eff2d_ana", "d3dlbtruth", 20, 0, .3, "plb", 20, 0, 30, "", cutLb.getCut());
    doQuickEff2d(events, "eff2d_ana_wide", "d3dlbtruth", 30, 0, 1, "plb", 30, 0, 100, "", cutLb.getCut());

    doQuickEff2d(events, "eff2d_ana_barrel", "d3dlbtruth", 20, 0, .3, "plb", 20, 0, 30, "", concatCutString(cutLb.getCut(), trigPassBarrel));
    doQuickEff2d(events, "eff2d_ana_barrel_wide", "d3dlbtruth", 30, 0, 1, "plb", 30, 0, 100, "", concatCutString(cutLb.getCut(), trigPassBarrel));

    doQuickEff2d(events, "eff2d_ana_displ", "d3dlbtruth", 20, 0, .3, "plb", 20, 0, 30, "", concatCutString(cutLb.getCut(), trigPassDispl));
    doQuickEff2d(events, "eff2d_ana_displ_wide", "d3dlbtruth", 30, 0, 1, "plb", 30, 0, 100, "", concatCutString(cutLb.getCut(), trigPassDispl));
}

