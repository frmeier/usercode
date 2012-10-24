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

void quickMultiplot()
{
    setTDRStyle();
    for(int i=0; i<=12; i++)
    {
	TFile *f = TFile::Open(("../data/vrt_r460_472_MC_B0_acc_barrel_ptbins_B0_"+toString(i)+".root").c_str());
	TTree *t = (TTree*)f->Get("fittree");
	t->Draw(("t>>h"+toString(i)+"(500,-2e-12,18e-12)").c_str(),"",i!=0?"same":"");
    }
}

