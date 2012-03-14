#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"

#include "Cut.h"
#include "Cuts.h"

using std::string;
using std::cout;
using std::endl;

void dumpSome(string infile, string field, string obj = "lb", string addCut = "")
{
    TFile* tfIn = TFile::Open(infile.c_str(), "READ");
    if (tfIn == 0)
    {
	cout << "Problem opening infile \"" << infile << "\" - exting";
	return;
    }
    TTree* intree = (TTree*) tfIn->Get("events");
    // cut selection
    Cuts cut;
    if (obj == "lb") cut.selectCut("acc03","lb07");
    if (obj == "B0") cut.selectCut("acc03B0","B001");
    //if (obj == "B0") cut.selectCut("acc03B0","B001exp");

    if (addCut.size() > 0)
	intree->Draw(">>lst", (cut.getCut()+"&&"+addCut).c_str());
    else
	intree->Draw(">>lst", cut.getCut().c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");
    intree->SetEventList(lst);

    intree->SetScanField(-1);
    intree->Scan(field.c_str());
}

