#include<set>
#include<iostream>
#include<fstream>
#include<string>
#include<utility>
#include<map>

#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TEventList.h"
#include "TH1F.h"

using std::string;
using std::cout;
using std::endl;

class Listentry
{
    public:
	Listentry() {};
	Listentry(Int_t a, Int_t b)
	{
	    entry_.first = a;
	    entry_.second = b;
	};
	Int_t first() const {return entry_.first; }
	Int_t second() const {return entry_.second; }
	bool operator()(const Listentry &a, const Listentry &b)
	{
	    if (a.first() == b.first()) return (a.second() < b.second());
	    else return (a.first() < b.first());
	};
    private:
	std::pair<Int_t, Int_t> entry_; // first: LS, second: evt
};


void plotGodBad(string filename, string listfilename, const std::string toPlot, const int nBins, const double lo, const double hi, const std::string divideBy = "")
{
    // read file with badlist
    ifstream infile(listfilename.c_str());
    Int_t ls, evt;
    std::set<Listentry, Listentry> badlist;
    while (!infile.eof())
    {
	infile >> ls >> evt;
	badlist.insert(Listentry(ls,evt));
    }
    cout << "badlist.size(): " << badlist.size() << endl;

    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
    }
    TTree *tree = (TTree*)f->Get("events");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
    }
    tree->Draw(">>lst", "");
    TEventList *lst = (TEventList*)gDirectory->Get("lst");
    cout << "Selection has " << lst->GetN() << " entries" << endl;

    tree->SetBranchAddress("LS", &ls);
    tree->SetBranchAddress("event", &evt);
    Double_t val, divisor;
    tree->SetBranchAddress(toPlot.c_str(), &val);
    const bool doDivision = (divideBy.size() != 0);
    if (doDivision) tree->SetBranchAddress(divideBy.c_str(), &divisor);

    std::string title = toPlot;
    if (doDivision) title = toPlot +"/"+ divideBy;
    TH1F *hgood = new TH1F("hgood", ("hgood "+title).c_str(), nBins, lo, hi);
    TH1F *hbad = new TH1F("hbad", "hbad", nBins, lo, hi);

    for(Long64_t i = 0; i!=lst->GetN(); i++)
    {
	tree->GetEntry(lst->GetEntry(i));
	if (badlist.find(Listentry(ls,evt)) != badlist.end())
	    hbad->Fill(doDivision ? val/divisor : val);
	else
	    hgood->Fill(doDivision ? val/divisor : val);
    }
    hgood->Draw();
    hbad->Draw("same");
    hgood->SetLineColor(3);
    hbad->SetLineColor(2);
    hbad->Scale(hgood->GetEntries()/hbad->GetEntries());
    if (hbad->GetMaximum() > hgood->GetMaximum()) hbad->SetMaximum(hbad->GetMaximum());
    double max = (hbad->GetMaximum() > hgood->GetMaximum()) ? hbad->GetMaximum() : hgood->GetMaximum();
    hgood->SetMaximum(max);
    hbad->SetMaximum(max);
    cout << "hbad->GetMaximum(): " << hbad->GetMaximum() << endl;
    cout << "hgood->GetMaximum(): " << hgood->GetMaximum() << endl;
    gPad->Update();
}

