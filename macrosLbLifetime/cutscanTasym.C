//#include "cut.C"
#include "cuts.C"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TMath.h"
#include "TEventList.h"

#include "setTDRStyle_modified.C"
#include "utils.h"

using std::cout;
using std::endl;

// Set verbosity to stdout
const int fVerbose(0);
const bool isB0(false);

TCanvas *c1, *c2;

typedef cutBase::valueType valueType;
typedef std::vector<cutBase::valueType> parVecType;

// a struct to collect the comprehensive results from a cutscan
struct cutscanresult
{
    // even thoufg a vect of some variables would be nicer, this is easier to plug into TGraps whithout transposing the vector<struct>
    std::vector<valueType> par;
    std::vector<valueType> asym;
    std::string cutname;
    std::string cuttitle;
    void reset() { par.clear(); asym.clear(); }
};

// a struct to hold the characteristic settings for one cutscan
struct cutscanitem
{
    cutscanitem(unsigned int cutno, unsigned int nsteps, valueType min, valueType max) :
	cutno_(cutno), nsteps_(nsteps), min_(min), max_(max) {};
    unsigned int cutno_;
    unsigned int nsteps_;
    valueType min_;
    valueType max_;
};

// type for vector of cutscanitems
typedef std::vector<cutscanitem> cutscanItemVecType;

// perform one cutscan, i.e. scans one variable in the defined range. All other cuts stay on their default value
cutscanresult doScan(TTree *tree, TCanvas *c, std::string hname, const cutSet &cs, parVecType pv,
	const unsigned int cutToVary, const unsigned int Ncuts, const valueType min, const valueType max)
    // pv is an intentional copy as this method changes it
{
    if (fVerbose > 0) cout << "cutToVary: " << cutToVary << endl;
    cutscanresult csr;
    const valueType cutslice = Ncuts>1 ? (max-min)/(Ncuts-1) : 0;
    // calculate an 'optimal' size of the canvas
    const int nHistos = Ncuts;
    const int cSizeY = TMath::Ceil(sqrt(nHistos));
    const int cSizeX = TMath::Ceil((double)nHistos/(double)cSizeY);
    c->Divide(cSizeX,cSizeY);

    Double_t sigLo, sigHi, bgrLo, nBgr;
    std::string massstring, timestring;
    if (isB0)
	{ sigLo = 5.24897; sigHi = 5.30833; bgrLo = 5.32323; nBgr = 710.5; massstring = "mB0"; timestring = "ct3dB0";}
    else
	{ sigLo = 5.57042; sigHi = 5.66637; bgrLo = 5.47359; nBgr = 685.264; massstring = "mlb"; timestring = "ct3dlb";}
    string sigCut = massstring + ">" + toString(sigLo) + "&&" + massstring + "<" + toString(sigHi);
    string bgrCut = massstring + ">" + toString(bgrLo);
    string hstring = "(80,-5e-12,15e-12)";
    TEventList* lst;
    tree->Draw(">>lst",cs.getCutExceptOne(pv,cutToVary).c_str());
    lst = (TEventList*)gDirectory->Get("lst");
    tree->SetEventList(lst);
    for(unsigned int i=0; i!=Ncuts; i++)
    {
	c->cd(i+1);
	// calculate current cut
	const valueType cutval = min+i*cutslice;
	pv[cutToVary] = cutval;
	if (fVerbose > 2) cout << cs.getCut(pv) << endl;
	// draw histo
	std::string hsig_name = "hsig"+hname+toString(i);
	std::string hbgr_name = "hbgr"+hname+toString(i);
	std::string hbgri_name = "hbgri"+hname+toString(i);
	cout << ("ct3dlb>>"+hsig_name+hstring) << endl;

	tree->Draw((timestring+">>"+hsig_name+hstring).c_str(), (cs.getCut(pv) + "&&" + sigCut).c_str());
	tree->Draw((timestring+">>"+hbgr_name+hstring).c_str(), (cs.getCut(pv) + "&&" + bgrCut).c_str(), "same");
	tree->Draw(("-"+timestring+">>"+hbgri_name+hstring).c_str(), (cs.getCut(pv) + "&&" + bgrCut).c_str(), "same");

	TH1F *hsig = (TH1F*)gDirectory->GetList()->FindObject(hsig_name.c_str());
	TH1F *hbgr = (TH1F*)gDirectory->GetList()->FindObject(hbgr_name.c_str());
	TH1F *hbgri= (TH1F*)gDirectory->GetList()->FindObject(hbgri_name.c_str());

	hsig->SetLineColor(8);
	hsig->SetLineWidth(1);
	hbgr->SetLineColor(9);
	hbgr->SetLineWidth(1);
	hbgri->SetLineColor(9);
	hbgri->SetLineStyle(2);
	hbgri->SetLineWidth(1);

	hbgr->Scale(nBgr/hbgr->GetEntries());
	hbgri->Scale(nBgr/hbgr->GetEntries());

	hsig->SetTitle(((isB0 ? "B^{0} " : "#Lambda_{b} ") + cs.getCutTitle(cutToVary) + " " + toString(cutval)).c_str());
	hsig->GetXaxis()->SetTitle("Lifetime (ps)");
	hsig->GetYaxis()->SetTitle("Entries");
	gPad->SetLogy();
	gPad->SetRightMargin(0.1);
	gPad->SetTopMargin(0.1);

	// calculate asymmetry
	const double right = hbgr->Integral(hbgr->FindBin(0),hbgr->GetNbinsX());
	const double left  = hbgri->Integral(hbgr->FindBin(0),hbgr->GetNbinsX());
	const double asym = (right-left)/(right+left);
	cout << "Asym: " << asym << endl;
	if (!isnan(asym))
	{
	    csr.par.push_back(cutval);
	    csr.asym.push_back(asym);
	}
    }
    tree->SetEventList(0);
    csr.cutname = cs.getCutName(cutToVary);
    csr.cuttitle = cs.getCutTitle(cutToVary);
    return csr;
}

// wrapper to invoke using a cutscanitem struct
cutscanresult doScan(TTree *t, TCanvas *c, std::string hname, const cutSet &cs, parVecType pv, cutscanitem csl)
{
    return doScan(t, c, hname, cs, pv, csl.cutno_, csl.nsteps_, csl.min_, csl.max_);
}

// draws a TGraph of a vector with data
void drawGraph(TPad* pad, parVecType x, parVecType y, std::string title, std::string xtitle, std::string ytitle, valueType ymin, valueType ymax)
{
    pad->cd();
    pad->SetLeftMargin(0.15);
    pad->SetTopMargin(0.10);
    TGraph *gr = new TGraph(x.size(),&x[0],&y[0]);
    gr->Draw("AP");
    gr->SetTitle(title.c_str());
    gr->GetXaxis()->SetTitle(xtitle.c_str());
    gr->GetYaxis()->SetTitle(ytitle.c_str());
    gr->SetMinimum(ymin);
    gr->SetMaximum(ymax);
}
void drawGraph(TPad* pad, parVecType x, parVecType y, std::string title, std::string xtitle, std::string ytitle)
{
    drawGraph(pad,x,y,title,xtitle,ytitle,0,-1111); // -1111 is the value to force root to autoscale....
}

// draws a TGraph of a vector with data
void drawGraphCutline(TPad* pad, parVecType x, parVecType y, std::string title, std::string xtitle, std::string ytitle, valueType xcut, valueType ymin, valueType ymax)
{
    pad->cd();
    pad->SetLeftMargin(0.15);
    pad->SetTopMargin(0.10);
    TGraph *gr = new TGraph(x.size(),&x[0],&y[0]);
    gr->Draw("AP");
    gr->SetTitle(title.c_str());
    gr->GetXaxis()->SetTitle(xtitle.c_str());
    gr->GetYaxis()->SetTitle(ytitle.c_str());
    gr->SetMinimum(ymin);
    gr->SetMaximum(ymax);

    // draw a vertical line at the given x position
    const valueType fractionOfHeight = .96;
    const valueType curMax = *(std::max_element(y.begin(), y.end()));
    const valueType curMin = ymin;
    const valueType height = curMax-curMin ;
    const valueType lineMin = curMin + (1-fractionOfHeight) * height;
    const valueType lineMax = curMax - (1-fractionOfHeight) * height;

    cout << "curMax: " << curMax << endl;
    TLine *l = new TLine(xcut,lineMin,xcut,lineMax);
    l->SetLineColor(2);
    l->Draw();
}
void drawGraphCutline(TPad* pad, parVecType x, parVecType y, std::string title, std::string xtitle, std::string ytitle, valueType xcut)
{
    drawGraphCutline(pad,x,y,title,xtitle,ytitle,xcut,0,-1111); // -1111 is the value to force root to autoscale....
}

// draws a TGraph of a vector with data
void drawGraphErrors(TPad* pad, parVecType x, parVecType y, parVecType yerr, std::string title, std::string xtitle, std::string ytitle, valueType ymin, valueType ymax)
{
    pad->cd();
    pad->SetLeftMargin(0.15);
    pad->SetTopMargin(0.10);
    TGraph *gr = new TGraphErrors(x.size(),&x[0],&y[0],0,&yerr[0]);
    gr->Draw("AP");
    gr->SetTitle(title.c_str());
    gr->GetXaxis()->SetTitle(xtitle.c_str());
    gr->GetYaxis()->SetTitle(ytitle.c_str());
    gr->SetMinimum(ymin);
    gr->SetMaximum(ymax);
}
void drawGraphErrors(TPad* pad, parVecType x, parVecType y, parVecType yerr, std::string title, std::string xtitle, std::string ytitle)
{
    drawGraphErrors(pad,x,y,yerr,title,xtitle,ytitle,0,-1111); // -1111 is the value to force root to autoscale....
}

// main function
void cutscan(std::string fullPath)
{
    // Style issues
    setTDRStyle();
    // initialize the canvases
    c1 = new TCanvas("c1","c1",1000,800);
    c2 = new TCanvas("c2","c2",8000,5000);

    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;

    // set what variable from the TTree to plot
    const std::string plotstring = "mlb";

    // Select cuts to use
    Cuts cuts;
    cuts.selectCut("acc04Lb","HLT_jpsiBarrel","lb07","iso01");

    // Cutstring
    cutSet cs = cuts.cs;
    parVecType parvec = cuts.parvec;

    // define the list of cutscans to perform
    cutscanItemVecType cutlist;
    cutlist.push_back(cutscanitem(cs.getCutPos("mjp"),12,0.060,0.300));
    cutlist.push_back(cutscanitem(cs.getCutPos("rpt1m"),16,2.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("rpt2m"),16,2.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptjp"),11,5,15));
    cutlist.push_back(cutscanitem(cs.getCutPos("ml0"),20,0.0001,0.010));
    cutlist.push_back(cutscanitem(cs.getCutPos("Kshypo"),10,0.002,0.020));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptl0"),13,2.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("rptpr"),11,1.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("rptpi"),12,0.1,3.3));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptgangDRl0"),26,0.005,0.060));
    cutlist.push_back(cutscanitem(cs.getCutPos("d3l0"),12,0.5,10));
    cutlist.push_back(cutscanitem(cs.getCutPos("d3l0/d3El0"),12,1.0,16));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptlb"),16,5.0,15.0));

    cutlist.push_back(cutscanitem(cs.getCutPos("isoClostrk"),21,0,22));
    cutlist.push_back(cutscanitem(cs.getCutPos("isoDocatrk"),21,0,0.2));

    if(fVerbose>-1) cout << "Current full cut string is:" << endl;
    if(fVerbose>-1) cout << cs.getCut(parvec) << endl;

    // Get TTree
    TTree* t = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;

    // set up the canvas for the summary graphs
    unsigned int nGraphs(1); // no. of TGraphs per scan
    unsigned int nGrDrawn(0); // counter for the graphs drawn so far
    c1->Clear();
    const double c1w = nGraphs*200;
    const double c1h = cutlist.size()*200;
    c1->SetWindowSize(c1w, c1h);
    //c1->SetCanvasSize(c1w, c1h);
    //c1->Divide(nGraphs,cutlist.size());
    const int nHistos = nGraphs*cutlist.size();
    const int cSizeY = TMath::Ceil(sqrt(nHistos));
    const int cSizeX = TMath::Ceil((double)nHistos/(double)cSizeY);
    c1->Divide(cSizeX,cSizeY);
    c1->Resize();

    // Initialize printing.
    const std::string pdfNameHistos = "cutscanTasymHistos";

    unsigned int canvascounter(0);
    // loop and do all cutscans in the list
    int itemno(0);
    for(cutscanItemVecType::iterator it = cutlist.begin(); it!=cutlist.end(); it++)
    {
	c2->Clear();
	const int c2w = 2000;
	const int c2h = (int)(c2w / sqrt(2));
	c2->SetCanvasSize(c2w, c2h);
	// actually do the cutscan
	cutscanresult csr = doScan(t, c2, "cs"+toString(itemno), cs, parvec, (*it));
	c2->SaveAs((pdfNameHistos+toString(canvascounter)+".pdf").c_str());
	const std::string title = "Cutscan results for " + csr.cutname;
	// plot graphs
	if(csr.par.size()>0)
	{
	    const valueType stdCutVal = parvec[(*it).cutno_];
	    cout << "stdCutVal: " << stdCutVal << endl;
	    drawGraphCutline((TPad*)c1->cd(nGrDrawn+1), csr.par, csr.asym, title, "Value for " + csr.cuttitle, "asymmetry",stdCutVal);
	}
	nGrDrawn+=nGraphs;
	canvascounter++;
	itemno++;
    }
    c1->SaveAs((pdfNameHistos+".pdf").c_str());
}

