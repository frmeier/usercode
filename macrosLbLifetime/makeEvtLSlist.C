#include<iostream>
#include<string>

#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"

using std::cout;
using std::endl;
using std::string;

void makeEvtLSlist(string filename)
{
    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
    }
    TTree *tree = (TTree*)f->Get("genevents");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
    }
    tree->Draw(">>lst", "hasCand==1");
    TEventList *lst = (TEventList*)gDirectory->Get("lst");
    cout << "Selection has " << lst->GetN() << " entries" << endl;

    Int_t evt, ls;
    tree->SetBranchAddress("LS", &ls);
    tree->SetBranchAddress("event", &evt);
    for(Long64_t i = 0; i!=lst->GetN(); i++)
    {
	tree->GetEntry(lst->GetEntry(i));
	cout << ls << " " << evt << endl;
    }
}

