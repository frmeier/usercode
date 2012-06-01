#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"

#include "Cuts.C"

void printResult(TH2F *h, string title)
{
    double N = (double)h->GetEntries();
    cout << "Results for " << title << ":" << endl;
    cout << "          | HLTmatch" << endl;
    cout << "          | false    true" << endl;
    cout << "----------+-----------------------" << endl;
    cout << "Trg false |" << h->GetBinContent(1,1) << "  " << h->GetBinContent(2,1) << "        " << h->GetBinContent(1,1)/N << "  " << h->GetBinContent(2,1)/N << endl;
    cout << "    true  |" << h->GetBinContent(1,2) << "  " << h->GetBinContent(2,2) << "        " << h->GetBinContent(1,2)/N << "  " << h->GetBinContent(2,2)/N <<endl;
    cout << endl;
}

void quickHLTmatchValidator()
{
    TFile *f1 = TFile::Open("../data/run380.root");
    TFile *f2 = TFile::Open("../data/run429.root");
    TFile *f3 = TFile::Open("../data/run430.root");
    TFile *f4 = TFile::Open("../data/run431.root");
    TFile *f5 = TFile::Open("../data/run384.root");
    /*
    TFile *f1 = TFile::Open("../data/run393.root");
    TFile *f2 = TFile::Open("../data/run423.root");
    TFile *f3 = TFile::Open("../data/run424.root");
    TFile *f4 = TFile::Open("../data/run425.root");
    TFile *f5 = TFile::Open("../data/run397.root");
    */

    TTree *t1 = (TTree*) f1->Get("events");
    TTree *t2 = (TTree*) f2->Get("events");
    TTree *t3 = (TTree*) f3->Get("events");
    TTree *t4 = (TTree*) f4->Get("events");
    TTree *t5 = (TTree*) f5->Get("events");
    
    //const string drawstring = "(25,-.5e-12,4.5e-12)";
    const string drawstring = "";

    const string toDraw = "HLTokBarrelJpsi:HLTmatch";
    //const string toDraw = "HLTokDisplJpsi:HLTmatch";

    t5->Draw((toDraw+">>h5"+drawstring).c_str(), "","lego");
    t1->Draw((toDraw+">>h1"+drawstring).c_str(), "","lego");
    t2->Draw((toDraw+">>h2"+drawstring).c_str(), "","lego");
    t3->Draw((toDraw+">>h3"+drawstring).c_str(), "","lego");
    t4->Draw((toDraw+">>h4"+drawstring).c_str(), "","lego");

    TH2F *h1 = (TH2F*)gDirectory->GetList()->FindObject("h1");
    TH2F *h2 = (TH2F*)gDirectory->GetList()->FindObject("h2");
    TH2F *h3 = (TH2F*)gDirectory->GetList()->FindObject("h3");
    TH2F *h4 = (TH2F*)gDirectory->GetList()->FindObject("h4");
    TH2F *h5 = (TH2F*)gDirectory->GetList()->FindObject("h5");

    /*
    printResult(h1, "run393");
    printResult(h2, "run423");
    printResult(h3, "run424");
    printResult(h4, "run425");
    printResult(h5, "run397");
    */
    printResult(h1, "run380");
    printResult(h2, "run429");
    printResult(h3, "run430");
    printResult(h4, "run431");
    printResult(h5, "run384");
}

