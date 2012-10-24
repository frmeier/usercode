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

void quickplot(string filename)
{
    TH1F *h1 = new TH1F("h1","h1",40,-2e-12,18e-12);
    
    TFile *f = TFile::Open(filename.c_str());
    TTree *t = (TTree*)f->Get("fittree");

    double time, eff, mass;
    t->SetBranchAddress("mass", &mass);
    t->SetBranchAddress("t", &time);
    t->SetBranchAddress("eff", &eff);

    const int N = t->GetEntries();
    for (int i=0; i!=N; i++)
    {
	t->GetEntry(i);
	if (eff <= 0) continue;
	if (mass < 5.6) continue;
	if (mass > 5.64) continue;
	h1->Fill(time, 1./eff);
    }
    h1->Draw();
}

