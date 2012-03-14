#include <iostream>
#include <string>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"

#include "utils.h"
#include "setTDRStyle_modified.C"

using std::cout;
using std::endl;
using std::string;

struct runlsevent
{
    runlsevent(int r, int l, int e) : run(r), ls(l), event(e) {};
    int run, ls, event;
    bool operator<(const runlsevent &rhs) const
    {
	if (run<rhs.run) return true;
	if (ls<rhs.ls) return true;
	if (event<rhs.event) return true;
	return false;
    };
};

void checkDoublecounts(string filename)
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
    int run, ls, event;
    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("LS",&ls);
    tree->SetBranchAddress("event",&event);
    const int nEvents = tree->GetEntries();
    std::map<runlsevent,int> eventmap;
    for (int i = 0; i!=nEvents; i++)
    {
	tree->GetEvent(i);
	eventmap[runlsevent(run,ls,event)]++;
    }
    cout << "nEvents: " << nEvents << " Mapsize: " << eventmap.size() << endl;
    if (nEvents == eventmap.size()) cout << "No doublecounts found" << endl;
    else cout << "OOOOPS!!!! There must be some doublecounts." << endl;
}


