#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLegend.h"

#include "Cuts.C"

void quickRunComparison()
{
    /*
    TFile *f1 = TFile::Open("../data/run380.root");
    TFile *f2 = TFile::Open("../data/run381.root");
    TFile *f3 = TFile::Open("../data/run382.root");
    TFile *f4 = TFile::Open("../data/run383.root");
    TFile *f5 = TFile::Open("../data/run384.root");
    */
    TFile *f1 = TFile::Open("../data/run393.root");
    TFile *f2 = TFile::Open("../data/run394.root");
    TFile *f3 = TFile::Open("../data/run395.root");
    TFile *f4 = TFile::Open("../data/run396.root");
    TFile *f5 = TFile::Open("../data/run397.root");

    TTree *t1 = (TTree*) f1->Get("events");
    TTree *t2 = (TTree*) f2->Get("events");
    TTree *t3 = (TTree*) f3->Get("events");
    TTree *t4 = (TTree*) f4->Get("events");
    TTree *t5 = (TTree*) f5->Get("events");
    
    Cuts cut;
    //cut.selectCut("acc05B0", "muSoft", "B004", "HLT_jpsiDispl", "HLT_matched");
    //cut.selectCut("acc05Lb", "muSoft", "lb11", "HLT_jpsiBarrel");
    cut.selectCut("acc05B0", "muSoft", "B004", "HLT_jpsiBarrel");

    //const string drawstring = "(25,-.5e-12,4.5e-12)";
    const string drawstring = "(100,-.5e-12,2e-12)";

    const string toDraw = "ct3dB0";

    t5->Draw((toDraw+">>h5"+drawstring).c_str(), cut.getCut().c_str());
    t1->Draw((toDraw+">>h1"+drawstring).c_str(), cut.getCut().c_str(),"same");
    t2->Draw((toDraw+">>h2"+drawstring).c_str(), cut.getCut().c_str(),"same");
    t3->Draw((toDraw+">>h3"+drawstring).c_str(), cut.getCut().c_str(),"same");
    t4->Draw((toDraw+">>h4"+drawstring).c_str(), cut.getCut().c_str(),"same");

    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("h1");
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject("h2");
    TH1F *h3 = (TH1F*)gDirectory->GetList()->FindObject("h3");
    TH1F *h4 = (TH1F*)gDirectory->GetList()->FindObject("h4");
    TH1F *h5 = (TH1F*)gDirectory->GetList()->FindObject("h5");

    /*
    h1->SetMinimum(0);
    h2->SetMinimum(0);
    h3->SetMinimum(0);
    h4->SetMinimum(0);
    h5->SetMinimum(0);
    */

    h1->SetLineColor(2);
    h2->SetLineColor(3);
    h3->SetLineColor(4);
    h4->SetLineColor(6);
    h5->SetLineColor(7);

    TLegend *leg = new TLegend(0.4,0.6,0.89,0.89);
    leg->AddEntry(h1,"Run2011A_MayReReco","l");
    leg->AddEntry(h2,"Run2011A_PromptReco_v4","l");
    leg->AddEntry(h3,"Run2011A_PromptReco_v5","l");
    leg->AddEntry(h4,"Run2011A_PromptReco_v6","l");
    leg->AddEntry(h5,"Run2011B_PromptReco_v1","l");
    leg->SetFillColor(0);
    leg->Draw();
}

