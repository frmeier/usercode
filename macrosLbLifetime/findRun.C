#include "TFile.h"
#include "TTree.h"
#include <string>
#include <iostream>
#include <map>

using std::string;
using std::cout;
using std::endl;

void getMaxRun(string filename)
{
    TFile *f=new TFile(filename.c_str());
    TTree *t = (TTree*)f->Get("events");
    const Long64_t N = t->GetEntries();
    std::map<Int_t,Int_t> runmap;
    Int_t currun;
    t->SetBranchAddress("run",&currun);
    for (Int_t i=0; i!=N; i++)
    {
	t->GetEntry(i);
	runmap[currun]++;
    }
    cout << "No. of runs in Tree:" << runmap.size() << endl;
    for (std::map<Int_t,Int_t>::const_iterator it=runmap.begin(); it!=runmap.end(); it++)
    {
	cout << it->first << " " << it->second << endl;
    }
}

