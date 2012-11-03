#include "TROOT.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include <iostream>
#include "TMath.h"

#include "utils.h"
#include "Ran.h"

using std::cout;
using std::endl;

void toyCtRadial(unsigned int N)
{
    TCanvas *c1 = new TCanvas();
    const double ctauLb = 0.0427;
    const double ctauL0 = 7.89;

    TRandom3 *ran = new TRandom3();
    // init tree
    TTree *tree = new TTree("tree","tree with lifetimes");
    Double_t ctLb, ctL0, ctSum;
    tree->Branch("ctLb",&ctLb,"ctLb/D");
    tree->Branch("ctL0",&ctL0,"ctL0/D");
    tree->Branch("ctSum",&ctSum,"ctSum/D");

    // fill tree
    for (unsigned int i = 0; i!=N; i++)
    {
	ctLb = ran->Exp(ctauLb);
	ctL0 = ran->Exp(ctauL0);
	ctSum = ctLb+ctL0;

	tree->Fill();
    }

    const int nBins = 400;
    const double lo = 0;
    const double hi = 2;
    tree->Draw(makePlotsString("ctLb", "h1", nBins, lo, hi).c_str());
    tree->Draw(makePlotsString("ctL0", "h2", nBins, lo, hi).c_str(), "", "same");
    tree->Draw(makePlotsString("ctSum", "h3", nBins, lo, hi).c_str(), "", "same");
    tree->Draw(makePlotsString("ctLb", "h4", nBins, lo, hi).c_str(), "getProb(ctSum)>getrn()", "same");
    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("h1");
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject("h2");
    TH1F *h3 = (TH1F*)gDirectory->GetList()->FindObject("h3");
    TH1F *h4 = (TH1F*)gDirectory->GetList()->FindObject("h3");
    h1->SetLineColor(2);
    h2->SetLineColor(3);
    h3->SetLineColor(4);
    h1->SetLineColor(1);

    gPad->SetLogy();
}

