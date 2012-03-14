#include <iostream>
#include <string>
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEventList.h"
#include "cut.h"
#include "cuts.C"
#include "setTDRStyle_modified.C"

TH1F *h1, *h2, *h3, *h4;
TCanvas *c1, *c2;

typedef cutBase::valueType valueType;
typedef std::vector<cutBase::valueType> parVecType;

void cutRates(std::string fullPath, bool savePdf = false, std::string pdfnamebase = "")
{
    // Set verbosity to stdout
    const int fVerbose(0);
    // Style issues
    setTDRStyle();
    c1 = new TCanvas("c1","c1",1000,600);
    unsigned int noOfPads = 4;
    c2 = new TCanvas;
    c2->Divide(2,2);
    TPad* pad1= (TPad*)c1->cd();
    pad1->cd();
    pad1->SetLeftMargin(0.20);
    pad1->SetTopMargin(0.10);
    TPad* pad2;
    for(unsigned int i=0; i!=noOfPads; i++)
    {
	pad2 = (TPad*)c2->cd(i+1);
        pad2->cd();
	pad2->SetLeftMargin(0.20);
	pad2->SetTopMargin(0.10);
    }

    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;

    // get cuts
    Cuts cutConst, cutVary;
    //cutConst.selectCut("acc03","HLT_matched_01");
    cutVary.selectCut("lb07");

    // Cutstring
    cutSet cs; parVecType parvec;
    cs = cutVary.cs;
    parvec = cutVary.parvec;

    if(fVerbose>0) cout << "Current full cut string is:" << endl;
    if(fVerbose>0) cout << cs.getCut(parvec) << endl;

    // Get TTree
    TTree* tfull = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << tfull->GetEntries() << " entries" << endl;
    TEventList* lst;

    // check if we can make a preselection
    TFile* tmpfile = new TFile("tmp.root","RECREATE");
    if (tmpfile == 0)
    {
	cout << "temporary file creation failed - exiting" << endl;
	return;
    }
    cout << "tmpfile: " << tmpfile << endl;
    TTree* t = 0;
    bool usePreseltree(false);
    if(cutConst.getCut() != "")
    {
        cout << "Preselecting events..." << endl;
        t = tfull->CopyTree(cutConst.getCut().c_str());
        cout << "now " << t->GetEntries() << " of " << tfull->GetEntries()<< endl;
        delete tfull;
        usePreseltree = true; // needed to detach the preseltree from the original one
    }
    else
	t = tfull;

    // ---------------------------------------------------------------------
    // Single cut histo
    {
	cout << "Single cut histo" << endl;
	const int bincount = cs.size()+2;
	h1 = new TH1F("h1","Single cut rate",bincount,0,bincount);
	h1->GetYaxis()->SetTitle("Cut rate wrt. no cuts");
	// get entries without any cut
	h1->SetBinContent(1,t->GetEntries());
	h1->GetXaxis()->SetBinLabel(1,"no cut at all");
	// get entries for individual cuts
	for (unsigned int i=0; i!= cs.size(); i++)
	{
	    const std::string lstname = "lst" + toString(i);
	    const std::string curcut = cs.getOneCut(i,parvec[i]);
	    t->Draw((">>"+lstname).c_str(),curcut.c_str());
	    lst = (TEventList*)gDirectory->Get(lstname.c_str());
	    h1->SetBinContent(i+2,lst->GetN());
	    h1->GetXaxis()->SetBinLabel(i+2,curcut.c_str());
	    cout << "  " << i << "/" << cs.size() << ": " << lst->GetN() << " events" << endl;
	}
	// get entries with all cuts
	t->Draw(">>lst",cs.getCut(parvec).c_str());
	lst = (TEventList*)gDirectory->Get("lst");
	if(fVerbose>0) cout << "Got TEventList with " << lst->GetN() << " events" << endl;
	h1->SetBinContent(bincount,lst->GetN());
	h1->GetXaxis()->SetBinLabel(bincount,"all cuts");
	// normalize to max
	h1->Scale(1./h1->GetBinContent(1));
	// finally draw the histo
	bool setLog = (h1->GetBinContent(1)-h1->GetBinContent(bincount))>.9;
	c1->cd();
	pad1->SetLogx(setLog);
	h1->Draw("hbar2");
	if (savePdf) c1->SaveAs((pdfnamebase+"_h1.pdf").c_str());
	c2->cd(1);
	pad2 = (TPad*)c2->cd(1);
	pad2->SetLogx(setLog);
	h1->Draw("hbar2");
	cout << "Single cut histo finished ----" << endl;
    }

    // ---------------------------------------------------------------------
    // Cut histo except one cut
    {
	cout << "Cut histo except one cut" << endl;
	const int bincount = cs.size()+1;
	h2 = new TH1F("h2","Cut rate, one cut excluded",bincount,0,bincount);
	h2->GetYaxis()->SetTitle("Cut rate wrt. full cut set");
	// get entries for individual cuts
	for (unsigned int i=0; i!= cs.size(); i++)
	{
	    const std::string lstname = "lstexcl" + toString(i);
	    const std::string curcut = cs.getCutExceptOne(parvec,i);
	    if(fVerbose>2) cout << i << ": " << curcut << endl;
	    const std::string leftcut = cs.getOneCut(i,parvec[i]);
	    t->Draw((">>"+lstname).c_str(),curcut.c_str());
	    lst = (TEventList*)gDirectory->Get(lstname.c_str());
	    h2->SetBinContent(i+1,lst->GetN());
	    h2->GetXaxis()->SetBinLabel(i+1,leftcut.c_str());
	    cout << "  " << i << "/" << cs.size() << ": " << lst->GetN() << " events" << endl;
	}
	// get entries with all cuts
	t->Draw(">>lstall",cs.getCut(parvec).c_str());
	lst = (TEventList*)gDirectory->Get("lstall");
	h2->SetBinContent(bincount,lst->GetN());
	h2->GetXaxis()->SetBinLabel(bincount,"all cuts");
	// normalize to all cuts
	h2->Scale(1./h2->GetBinContent(bincount));
	// finally draw the histo
	bool setLog = (1./h2->GetBinContent(bincount))>5;
	c1->cd();
	pad1->SetLogx(setLog);
	h2->Draw("hbar2");
	if (savePdf) c1->SaveAs((pdfnamebase+"_h2.pdf").c_str());
	c2->cd(2);
	pad2 = (TPad*)c2->cd(2);
	pad2->SetLogx(setLog);
	h2->Draw("hbar2");
	cout << "Cut histo except one cut finished -------" << endl;
    }

    // ---------------------------------------------------------------------
    // Cut histo sequential application of cuts
    {
	cout << "Cut histo sequential application of cuts" << endl;
	const int bincount = cs.size()+1;
	h3 = new TH1F("h3","Cut rate, sequential application",bincount,0,bincount);
	h3->GetYaxis()->SetTitle("Cut rate wrt. no cuts");
	// now start the work
	cutSet workcs=cs;
	std::vector<cutBase::valueType> workparvec = parvec;
	unsigned int lstcounter(0); // counter to enumerate lists for probable later use
	std::string constcuts("");
	// get number of events for empty cutset
	t->Draw(">>lstseqnocut",constcuts.c_str());
	lst = (TEventList*)gDirectory->Get("lstseqnocut");
	unsigned int bincounter(bincount); // Start from top to down
	h3->SetBinContent(bincounter, lst->GetN());
	h3->GetXaxis()->SetBinLabel(bincounter,"no cuts");
	while (workcs.size() != 0)
	{
	    cout << "   " << workcs.size() << endl;
	    // using workcs.begin() ensures that at the end at least one cut will be removed
	    // even if no minimum has been found
	    double minValue = h3->GetBinContent(bincounter);
	    cutSet::iterator minIt = workcs.begin();
	    unsigned int minInt = 0;
	    std::string minCut = workcs.getOneCut(0,workparvec[0]);
	    unsigned int loopcounter = 0;
	    if (constcuts.size()!=0) constcuts+="&&";
	    // now loop over all cuts, add them and check the outcome
	    for (cutSet::iterator it = workcs.begin(); it!=workcs.end(); it++)
	    {
		const std::string lstname = "lstseq" + toString(lstcounter);
		const std::string curcut = workcs.getOneCut(loopcounter,workparvec[loopcounter]);
		if (fVerbose>1) cout << curcut << endl;
		t->Draw((">>"+lstname).c_str(),(constcuts+curcut).c_str());
		lst = (TEventList*)gDirectory->Get(lstname.c_str());
		if(lst->GetN() < minValue)
		{
		    minValue = lst->GetN();
		    minIt = it;
		    minInt = loopcounter;
		    minCut = curcut;
		}
		lstcounter++;
		loopcounter++;
	    }
	    if (fVerbose>1) cout << endl;
	    workcs.erase(minIt);
	    workparvec.erase(workparvec.begin()+minInt);
	    bincounter--;
	    h3->SetBinContent(bincounter, minValue);
	    h3->GetXaxis()->SetBinLabel(bincounter,("+" + minCut).c_str());
	    constcuts += minCut;
	}
	// normalize to all cuts
	h3->Scale(1./h3->GetBinContent(bincount));
	// finally draw the histo
	bool setLog = (h3->GetBinContent(bincount)-h3->GetBinContent(1))>.9;
	c1->cd();
	pad1->SetLogx(setLog);
	h3->Draw("hbar2");
	if (savePdf) c1->SaveAs((pdfnamebase+"_h3.pdf").c_str());
	pad2 = (TPad*)c2->cd(3);
	pad2->SetLogx(setLog);
	h3->Draw("hbar2");
	cout << "Cut histo sequential application of cuts finished -------" << endl;
    }

    // ---------------------------------------------------------------------
    // Cut histo sequential removal of cuts
    {
	cout << "Cut histo sequential removal of cuts" << endl;
	const int bincount = cs.size()+1;
	h4 = new TH1F("h4","Cut rate, sequential removal",bincount,0,bincount);
	h4->GetYaxis()->SetTitle("Cut rate wrt. no cuts");
	// now start the work
	cutSet workcs=cs;
	std::vector<cutBase::valueType> workparvec = parvec;
	unsigned int lstcounter(0); // counter to enumerate lists for probable later use
	// get number of events for empty cutset
	t->Draw(">>lstseqremallcut",workcs.getCut(workparvec).c_str());
	lst = (TEventList*)gDirectory->Get("lstseqremallcut");
	unsigned int bincounter(bincount); // Start from top to down
	h4->SetBinContent(bincounter, lst->GetN());
	h4->GetXaxis()->SetBinLabel(bincounter,"all cuts");
	while (workcs.size() != 0)
	{
	    cout << "   " << workcs.size() << endl;
	    // using workcs.begin() ensures that at the end at least one cut will be removed
	    // even if no maximum has been found
	    double maxValue = h4->GetBinContent(bincounter);
	    cutSet::iterator maxIt = workcs.begin();
	    unsigned int maxInt = 0;
	    std::string maxCut = workcs.getOneCut(0,workparvec[0]);
	    unsigned int loopcounter = 0;
	    // now loop over all cuts, add them and check the outcome
	    for (cutSet::iterator it = workcs.begin(); it!=workcs.end(); it++)
	    {
		const std::string lstname = "lstseqrem" + toString(lstcounter);
		const std::string curcut = workcs.getCutExceptOne(workparvec,loopcounter);
		if (fVerbose>1) cout << curcut << endl;
		t->Draw((">>"+lstname).c_str(),curcut.c_str());
		lst = (TEventList*)gDirectory->Get(lstname.c_str());
		if(lst->GetN() > maxValue)
		{
		    maxValue = lst->GetN();
		    maxIt = it;
		    maxInt = loopcounter;
		    maxCut = workcs.getOneCut(loopcounter,workparvec[loopcounter]);
		}
		lstcounter++;
		loopcounter++;
	    }
	    if (fVerbose>1) cout << endl;
	    workcs.erase(maxIt);
	    workparvec.erase(workparvec.begin()+maxInt);
	    bincounter--;
	    h4->SetBinContent(bincounter, maxValue);
	    h4->GetXaxis()->SetBinLabel(bincounter,("-" + maxCut).c_str());
	}
	// normalize to all cuts
	h4->Scale(1./h4->GetBinContent(1));
	// finally draw the histo
	bool setLog = (h4->GetBinContent(1)-h4->GetBinContent(bincount))>.9;
	c1->cd();
	pad1->SetLogx(setLog);
	h4->Draw("hbar2");
	if (savePdf) c1->SaveAs((pdfnamebase+"_h4.pdf").c_str());
	pad2 = (TPad*)c2->cd(4);
	pad2->SetLogx(setLog);
	h4->Draw("hbar2");
	cout << "Cut histo sequential removal of cuts -- finished" << endl;
    }
}

void cutRates()
{
    // Path to datafile
    std::string strPath = "/home/frank/psi/lambda/data/";
    std::string strFile = "../data/run041.root";
    cutRates(strPath+strFile);
}
