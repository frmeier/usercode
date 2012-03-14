#include <string>
#include <iostream>

#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TLegend.h"

#include "setTDRStyle_modified.C"

using std::string;
using std::cout;
using std::endl;

void studyEffCt(string filename)
{
    setTDRStyle();
    TFile *_file0 = TFile::Open(filename.c_str());

    TTree *gentree = (TTree*)gDirectory->Get("genevents");
    TTree *recotree = (TTree*)gDirectory->Get("events");

    TCanvas *c = new TCanvas("c","c",1200,800);
    c->Divide(2,2);
    int canvasctr(0);
    
    // Zerfallszeit
    c->cd(1);
    string hbin = "(40,0,20e-12)";
    gentree->Draw(("ctB0>>hgen_ct_pbins1"+hbin).c_str(),"pB0>4.52&&pB0<5.64");
    gentree->Draw(("ctB0>>hgen_ct_pbins2"+hbin).c_str(),"pB0>10.20&&pB0<10.63","same");
    gentree->Draw(("ctB0>>hgen_ct_pbins3"+hbin).c_str(),"pB0>24&&pB0<25","same");
    gentree->Draw(("ctB0>>hgen_ct_pbins4"+hbin).c_str(),"pB0>55&&pB0<67","same");

    TH1F *hgen_ct_pbins1 = (TH1F*)gDirectory->GetList()->FindObject("hgen_ct_pbins1");
    TH1F *hgen_ct_pbins2 = (TH1F*)gDirectory->GetList()->FindObject("hgen_ct_pbins2");
    TH1F *hgen_ct_pbins3 = (TH1F*)gDirectory->GetList()->FindObject("hgen_ct_pbins3");
    TH1F *hgen_ct_pbins4 = (TH1F*)gDirectory->GetList()->FindObject("hgen_ct_pbins4");

    Int_t n1 = hgen_ct_pbins1->GetEntries();
    hgen_ct_pbins2->Scale(n1/hgen_ct_pbins2->GetEntries());
    hgen_ct_pbins3->Scale(n1/hgen_ct_pbins3->GetEntries());
    hgen_ct_pbins4->Scale(n1/hgen_ct_pbins4->GetEntries());

    hgen_ct_pbins1->SetLineColor(1);
    hgen_ct_pbins2->SetLineColor(2);
    hgen_ct_pbins3->SetLineColor(3);
    hgen_ct_pbins4->SetLineColor(4);

    hgen_ct_pbins1->SetTitle("Lifetime MC truth generator level");
    hgen_ct_pbins1->GetXaxis()->SetTitle("Lifetime (s)");

    TLegend *leg_gen_ct_pbins = new TLegend(0.4,0.6,0.89,0.89);
    leg_gen_ct_pbins->AddEntry(hgen_ct_pbins1,"p #in (4.52,5.64)","lf");
    leg_gen_ct_pbins->AddEntry(hgen_ct_pbins2,"p #in (10.20,10.63)","lf");
    leg_gen_ct_pbins->AddEntry(hgen_ct_pbins3,"p #in (24.0,25.0)","lf");
    leg_gen_ct_pbins->AddEntry(hgen_ct_pbins4,"p #in (55.0,67.0)","lf");
    leg_gen_ct_pbins->Draw();

    gPad->SetLogy();


    // nun dsselbe fuer Flugstrecke
    c->cd(2);

    hbin = "(40,0,2.5)";
    gentree->Draw(("d3dB0>>hgen_d3_pbins1"+hbin).c_str(),"pB0>4.52&&pB0<5.64");
    gentree->Draw(("d3dB0>>hgen_d3_pbins2"+hbin).c_str(),"pB0>10.20&&pB0<10.63","same");
    gentree->Draw(("d3dB0>>hgen_d3_pbins3"+hbin).c_str(),"pB0>24&&pB0<25","same");
    gentree->Draw(("d3dB0>>hgen_d3_pbins4"+hbin).c_str(),"pB0>55&&pB0<67","same");

    TH1F *hgen_d3_pbins1 = (TH1F*)gDirectory->GetList()->FindObject("hgen_d3_pbins1");
    TH1F *hgen_d3_pbins2 = (TH1F*)gDirectory->GetList()->FindObject("hgen_d3_pbins2");
    TH1F *hgen_d3_pbins3 = (TH1F*)gDirectory->GetList()->FindObject("hgen_d3_pbins3");
    TH1F *hgen_d3_pbins4 = (TH1F*)gDirectory->GetList()->FindObject("hgen_d3_pbins4");

    n1 = hgen_d3_pbins1->GetEntries();
    hgen_d3_pbins2->Scale(n1/hgen_d3_pbins2->GetEntries());
    hgen_d3_pbins3->Scale(n1/hgen_d3_pbins3->GetEntries());
    hgen_d3_pbins4->Scale(n1/hgen_d3_pbins4->GetEntries());

    hgen_d3_pbins1->SetLineColor(1);
    hgen_d3_pbins2->SetLineColor(2);
    hgen_d3_pbins3->SetLineColor(3);
    hgen_d3_pbins4->SetLineColor(4);

    hgen_d3_pbins1->SetTitle("Lifetime MC truth generator level");
    hgen_d3_pbins1->GetXaxis()->SetTitle("Flight length (cm)");

    TLegend *leg_gen_d3_pbins = new TLegend(0.4,0.6,0.89,0.89);
    leg_gen_d3_pbins->AddEntry(hgen_d3_pbins1,"p #in (4.52,5.64)","lf");
    leg_gen_d3_pbins->AddEntry(hgen_d3_pbins2,"p #in (10.20,10.63)","lf");
    leg_gen_d3_pbins->AddEntry(hgen_d3_pbins3,"p #in (24.0,25.0)","lf");
    leg_gen_d3_pbins->AddEntry(hgen_d3_pbins4,"p #in (55.0,67.0)","lf");
    leg_gen_d3_pbins->Draw();

    gPad->SetLogy();

    // Zerfallszeit in reco
    /*
    c->cd(3);
    hbin = "(40,0,20e-12)";
    recotree->Draw(("ctB0truth>>hreco_ct_pbins1"+hbin).c_str(),"pB0truth>0.62&&pB0truth<10.92");
    recotree->Draw(("ctB0truth>>hreco_ct_pbins2"+hbin).c_str(),"pB0truth>10.92&&pB0truth<17.03","same");
    recotree->Draw(("ctB0truth>>hreco_ct_pbins3"+hbin).c_str(),"pB0truth>17.03&&pB0truth<26.47","same");
    recotree->Draw(("ctB0truth>>hreco_ct_pbins4"+hbin).c_str(),"pB0truth>26.47&&pB0truth<61.41","same");

    TH1F *hreco_ct_pbins1 = (TH1F*)gDirectory->GetList()->FindObject("hreco_ct_pbins1");
    TH1F *hreco_ct_pbins2 = (TH1F*)gDirectory->GetList()->FindObject("hreco_ct_pbins2");
    TH1F *hreco_ct_pbins3 = (TH1F*)gDirectory->GetList()->FindObject("hreco_ct_pbins3");
    TH1F *hreco_ct_pbins4 = (TH1F*)gDirectory->GetList()->FindObject("hreco_ct_pbins4");

    n1 = hreco_ct_pbins1->GetEntries();
    hreco_ct_pbins2->Scale(n1/hreco_ct_pbins2->GetEntries());
    hreco_ct_pbins3->Scale(n1/hreco_ct_pbins3->GetEntries());
    hreco_ct_pbins4->Scale(n1/hreco_ct_pbins4->GetEntries());

    hreco_ct_pbins1->SetLineColor(1);
    hreco_ct_pbins2->SetLineColor(2);
    hreco_ct_pbins3->SetLineColor(3);
    hreco_ct_pbins4->SetLineColor(4);

    hreco_ct_pbins1->SetTitle("Lifetime MC reco, MC truth");
    hreco_ct_pbins1->GetXaxis()->SetTitle("Lifetime (s)");

    TLegend *leg_reco_ct_pbins = new TLegend(0.4,0.6,0.89,0.89);
    leg_reco_ct_pbins->AddEntry(hreco_ct_pbins1,"p #in (0.62,10.97)","lf");
    leg_reco_ct_pbins->AddEntry(hreco_ct_pbins2,"p #in (10.97,17.03)","lf");
    leg_reco_ct_pbins->AddEntry(hreco_ct_pbins3,"p #in (17.03,26.47)","lf");
    leg_reco_ct_pbins->AddEntry(hreco_ct_pbins4,"p #in (26.47,61.41)","lf");
    leg_reco_ct_pbins->Draw();

    gPad->SetLogy();
*/
    // Zerfallszeit in reco
    c->cd(4);
    hbin = "(40,0,20e-12)";
    recotree->Draw(("ctB0truth>>hreco_ct_pbins1"+hbin).c_str(),"isMCmatch==1&&pB0truth>0.62&&pB0truth<10.92");
    recotree->Draw(("ctB0truth>>hreco_ct_pbins2"+hbin).c_str(),"isMCmatch==1&&pB0truth>10.92&&pB0truth<17.03","same");
    recotree->Draw(("ctB0truth>>hreco_ct_pbins3"+hbin).c_str(),"isMCmatch==1&&pB0truth>17.03&&pB0truth<26.47","same");
    recotree->Draw(("ctB0truth>>hreco_ct_pbins4"+hbin).c_str(),"isMCmatch==1&&pB0truth>26.47&&pB0truth<61.41","same");

    TH1F *hreco_ct_pbins1 = (TH1F*)gDirectory->GetList()->FindObject("hreco_ct_pbins1");
    TH1F *hreco_ct_pbins2 = (TH1F*)gDirectory->GetList()->FindObject("hreco_ct_pbins2");
    TH1F *hreco_ct_pbins3 = (TH1F*)gDirectory->GetList()->FindObject("hreco_ct_pbins3");
    TH1F *hreco_ct_pbins4 = (TH1F*)gDirectory->GetList()->FindObject("hreco_ct_pbins4");

    n1 = hreco_ct_pbins1->GetEntries();
    hreco_ct_pbins2->Scale(n1/hreco_ct_pbins2->GetEntries());
    hreco_ct_pbins3->Scale(n1/hreco_ct_pbins3->GetEntries());
    hreco_ct_pbins4->Scale(n1/hreco_ct_pbins4->GetEntries());

    hreco_ct_pbins1->SetLineColor(1);
    hreco_ct_pbins2->SetLineColor(2);
    hreco_ct_pbins3->SetLineColor(3);
    hreco_ct_pbins4->SetLineColor(4);

    hreco_ct_pbins1->SetTitle("Lifetime MC reco truthmatched, MC truth");
    hreco_ct_pbins1->GetXaxis()->SetTitle("Lifetime (s)");

    TLegend *leg_reco_ct_pbins = new TLegend(0.4,0.6,0.89,0.89);
    leg_reco_ct_pbins->AddEntry(hreco_ct_pbins1,"p #in (0.62,10.97)","lf");
    leg_reco_ct_pbins->AddEntry(hreco_ct_pbins2,"p #in (10.97,17.03)","lf");
    leg_reco_ct_pbins->AddEntry(hreco_ct_pbins3,"p #in (17.03,26.47)","lf");
    leg_reco_ct_pbins->AddEntry(hreco_ct_pbins4,"p #in (26.47,61.41)","lf");
    leg_reco_ct_pbins->Draw();

    gPad->SetLogy();
}

