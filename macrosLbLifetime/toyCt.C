#include "TROOT.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include <iostream>
#include "TMath.h"

using std::cout;
using std::endl;

void toyCt(unsigned int N)
{
    TCanvas *c1 = new TCanvas();
    const double tau = 1.5e-12;
    TRandom3 *ran = new TRandom3();
    // init tree
    TTree *tree = new TTree("tree","tree with lifetimes");
    Double_t t;
    tree->Branch("t",&t,"t/D");
    // fill tree
    const double mean(0), sigma(1e-12);
    for (unsigned int i = 0; i!=N; i++)
    {
	//t = ran->Exp(tau) + ran->Gaus(mean,sigma) - 3e-12;
	double cursigma = TMath::Abs(ran->Gaus(2e-12,.5e-12));
	double curt = ran->Exp(tau) + ran->Gaus(0,cursigma);
	t = curt - 3*cursigma;
	tree->Fill();
    }
    for (unsigned int i = 0; i!=N; i++)
    {
	//t = ran->Exp(tau) + ran->Gaus(mean,sigma) - 3e-12;
	double cursigma = TMath::Abs(ran->Gaus(1e-12,.1e-12));
	double curt = ran->Exp(tau) + ran->Gaus(0,cursigma);
	t = curt*1.1 - 3*cursigma;
	tree->Fill();
    }
    tree->Draw("t>>h1(40,0,10e-12)");
    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("h1");
    //TF1 *f1 = new TF1("f1","exp([0]+x/[1])",0,10e-12);
    TF1 *f1 = new TF1("f1","expo(0)",0,1e-12);
    f1->SetParameter(0,0);
    f1->SetParameter(1,-1e+12);
    h1->Fit(f1);
    cout << -1./f1->GetParameter(1) << "+/-" << -1./f1->GetParError(1) << endl;
}

