#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TEventList.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include "utils.h"
#include "Cuts.C"
#include "setTDRStyle_modified.C"

using std::string;

void drawP(TTree *t, TH1F *h, int charge, string cut)
{
    const double MPION = 0.1396;
    int qha1;
    double ptha1, phiha1, etaha1;
    double ptha2, phiha2, etaha2;
    t->SetBranchAddress("rqha1", &qha1);
    t->SetBranchAddress("rptha1", &ptha1);
    t->SetBranchAddress("retaha1", &etaha1);
    t->SetBranchAddress("rphiha1", &phiha1);
    t->SetBranchAddress("rptha2", &ptha2);
    t->SetBranchAddress("retaha2", &etaha2);
    t->SetBranchAddress("rphiha2", &phiha2);
    TLorentzVector tlv;

    t->Draw(">>lst", cut.c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");

    const int N = lst->GetN();
    //const int N = 10000;
    for (int i=0; i!=N; i++)
    {
	t->GetEntry(lst->GetEntry(i));
	if(qha1 * charge > 0)
	    tlv.SetPtEtaPhiM(ptha1, etaha1, phiha1, MPION);
	else
	    tlv.SetPtEtaPhiM(ptha2, etaha2, phiha2, MPION);
	h->Fill(tlv.P());
    }
    delete lst;
}

void quickPipPimPlot(string filename)
{
    setTDRStyle();
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    TFile *f = TFile::Open(filename.c_str()); // Data
    TTree *events = (TTree*)f->Get("events");

    c->Divide(1,2);
    c->cd(1);

    Cuts cut;
    //cut.selectCut("lb12", "acc06Lb", "muSoft", "HLT_jpsiBarrel", "HLT_matched");
    cut.selectCut("B006", "acc06B0", "muSoft", "HLT_jpsiBarrel", "HLT_matched");
    //const string curCut = cut.getCut() + "&&mbc>5.24&&mbc<5.32";
    const string curCut = "isMCmatch==1&&mbc>5.24&&mbc<5.32";

    const int nbins = 180;
    const double lo = 0.5, hi = 5.0;
    //const string hstring = "("+toString(nbins)+","+toString(lo)+","+toString(hi)+")";

    TH1F *h1 = new TH1F("h1","h1", nbins, lo, hi);
    //events->Draw("(rqha1>0?rptha1:rptha2)>>h1");
    drawP(events, h1, +1, curCut);
    h1->Draw();
    h1->Sumw2();
    h1->SetMarkerColor(2);
    h1->SetLineColor(2);

    TH1F *h2 = new TH1F("h2","h2", nbins, lo, hi);
    //events->Draw("(rqha1<0?rptha1:rptha2)>>h2","","same");
    drawP(events, h2, -1, curCut);
    h2->Draw("same");
    h2->Sumw2();
    h2->SetMarkerColor(3);
    h2->SetLineColor(3);

    //TLegend *leg = new TLegend(0.70,0.70, 0.90,0.89);
    TLegend *leg = new TLegend(.15, .70, .25, .90);
    leg->AddEntry(h1, "#pi^{+}", "lm");
    leg->AddEntry(h2, "#pi^{-}", "lm");
    leg->Draw();
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    c->cd(2);
    TH1F *hdiv = new TH1F("hdiv","#pi^{-}/#pi^{+}", nbins, lo, hi);
    hdiv->Divide(h2, h1);
    hdiv->Draw("E");
}



