//#include "cut.C"
#include "cuts.C"
#include "fitfunc.C"
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

using std::cout;
using std::endl;

// Set verbosity to stdout
const int fVerbose(0);

TCanvas *c1, *c2;

typedef cutBase::valueType valueType;
typedef std::vector<cutBase::valueType> parVecType;

// a struct to collect the comprehensive results from a cutscan
struct cutscanresult
{
    // even thoufg a vect of some variables would be nicer, this is easier to plug into TGraps whithout transposing the vector<struct>
    std::vector<valueType> par;
    std::vector<valueType> SoverSqrtSB;
    std::vector<valueType> sig;
    std::vector<valueType> chi2;
    std::vector<valueType> mean;
    std::vector<valueType> sigma;
    std::string cutname;
    std::string cuttitle;
    void reset() { par.clear(); SoverSqrtSB.clear(); }
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
cutscanresult doScan(TTree *t, TCanvas *c, TH1F *h, fitfuncBase *ff, const std::string toDraw, const cutSet &cs, parVecType pv, 
	const unsigned int cutToVary, const unsigned int Ncuts, const valueType min, const valueType max)
    // pv is an intentional copy as this method changes it
{
    if (fVerbose > 0) cout << "cutToVary: " << cutToVary << endl;
    cutscanresult csr;
    const valueType cutslice = Ncuts>1 ? (max-min)/(Ncuts-1) : 0;
    const std::string hname = toString(h->GetName());
    const std::string drawstring = toDraw + ">>" + hname;
    // calculate an 'optimal' size of the canvas
    const int nHistos = Ncuts;
    const int cSizeY = TMath::Ceil(sqrt(nHistos));
    const int cSizeX = TMath::Ceil((double)nHistos/(double)cSizeY);
    c->Divide(cSizeX,cSizeY);

    TEventList* lst;
    t->Draw(">>lst",cs.getCutExceptOne(pv,cutToVary).c_str());
    lst = (TEventList*)gDirectory->Get("lst");
    t->SetEventList(lst);
    for(unsigned int i=0; i!=Ncuts; i++)
    {
	c->cd(i+1);
	// calculate current cut
	const valueType cutval = min+i*cutslice;
	pv[cutToVary] = cutval;
	if (fVerbose > 2) cout << cs.getCut(pv) << endl;
	// draw histo
	t->Draw(drawstring.c_str(),cs.getCut(pv).c_str()); 
	h->SetTitle(cs.getOneCut(cutToVary,cutval).c_str());
	TH1F *newh = (TH1F*)h->Clone((drawstring+toString(i)).c_str());
	newh->SetMinimum(0);
	newh->Draw();
	// perform fit
	ff->reset();
	ff->doFit(newh);
	// fill result vectors
	if(ff->getSig()>=0 && ff->getBgr()>=0 && ff->getChi2()<70) // omit rubish entries
	{
	    csr.par.push_back(cutval);
	    csr.SoverSqrtSB.push_back(ff->getSoverSqrtSB());
	    csr.sig.push_back(ff->getSig());
	    csr.chi2.push_back(ff->getChi2());
	    csr.mean.push_back(ff->getParameter(3));
	    csr.sigma.push_back(ff->getParameter(4));
	}
	if (fVerbose > 2) cout << "Chi2 of this fit: " << ff->getChi2() << endl;
    }
    t->SetEventList(0);
    csr.cutname = cs.getCutName(cutToVary);
    csr.cuttitle = cs.getCutTitle(cutToVary);
    return csr;
}

// wrapper to invoke using a cutscanitem struct
cutscanresult doScan(TTree *t, TCanvas *c, TH1F *h, fitfuncBase *ff, const std::string toDraw, const cutSet &cs, parVecType pv, 
	cutscanitem csl)
{
    return doScan(t,c,h,ff,toDraw,cs,pv,csl.cutno_,csl.nsteps_,csl.min_,csl.max_);
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
    c2 = new TCanvas("c2","c2",1000,600);

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
    cuts.selectCut("acc04Lb","HLT_jpsiBarrel","lb07");

    // Cutstring
    cutSet cs = cuts.cs;
    parVecType parvec = cuts.parvec;

    // define the list of cutscans to perform
    cutscanItemVecType cutlist;
    
    /*
    cutlist.push_back(cutscanitem(cs.getCutPos("mjp"),12,0.115,0.445));
    cutlist.push_back(cutscanitem(cs.getCutPos("prob1m"),12,0.000,0.33));
    cutlist.push_back(cutscanitem(cs.getCutPos("prob2m"),12,0.000,0.33));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptjp"),12,0,11));
    cutlist.push_back(cutscanitem(cs.getCutPos("probjp"),12,0.000,0.11));

    cutlist.push_back(cutscanitem(cs.getCutPos("ml0"),12,0.005,0.027));
    cutlist.push_back(cutscanitem(cs.getCutPos("probpr"),12,0.00,0.033));
    cutlist.push_back(cutscanitem(cs.getCutPos("probpi"),12,0.00,0.033));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptl0"),12,0,5.5));
    cutlist.push_back(cutscanitem(cs.getCutPos("rptpr"),12,0,2.75));
    cutlist.push_back(cutscanitem(cs.getCutPos("rptpi"),12,0,1.375));

    */
    /* 25.5.11
    cutlist.push_back(cutscanitem(cs.getCutPos("mjp"),12,0.060,0.300));
    cutlist.push_back(cutscanitem(cs.getCutPos("rpt1m"),16,2.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("rpt2m"),16,2.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptjp"),11,5,15));
    cutlist.push_back(cutscanitem(cs.getCutPos("ml0"),12,0.005,0.027));
    cutlist.push_back(cutscanitem(cs.getCutPos("Kshypo"),10,0.002,0.020));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptl0"),13,2.0,3.2));
    cutlist.push_back(cutscanitem(cs.getCutPos("rptpr"),11,1.0,3.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("rptpi"),12,0.1,1.3));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptgangDRl0"),26,0.005,0.030));
    cutlist.push_back(cutscanitem(cs.getCutPos("d3l0"),12,3.0,36));
    cutlist.push_back(cutscanitem(cs.getCutPos("d3l0/d3El0"),12,3.0,36));
    */
    //cutlist.push_back(cutscanitem(cs.getCutPos("probjp"),21,0.0,0.10));
    //cutlist.push_back(cutscanitem(cs.getCutPos("probl0"),21,0.0,0.10));
    //cutlist.push_back(cutscanitem(cs.getCutPos("problb"),21,0.0,0.10));
    //cutlist.push_back(cutscanitem(cs.getCutPos("ptgangDRl0"),12,0.001,0.013));
    //cutlist.push_back(cutscanitem(cs.getCutPos("ptgangDRlb"),12,0.05,0.60));
    cutlist.push_back(cutscanitem(cs.getCutPos("mjp"),12,0.060,0.300));
    cutlist.push_back(cutscanitem(cs.getCutPos("rpt1m"),16,2.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("rpt2m"),16,2.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptjp"),11,5,15));
    cutlist.push_back(cutscanitem(cs.getCutPos("ml0"),12,0.002,0.027));
    cutlist.push_back(cutscanitem(cs.getCutPos("Kshypo"),10,0.002,0.020));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptl0"),13,2.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("rptpr"),11,1.0,5.0));
    cutlist.push_back(cutscanitem(cs.getCutPos("rptpi"),12,0.1,3.3));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptgangDRl0"),26,0.005,0.060));
    cutlist.push_back(cutscanitem(cs.getCutPos("d3l0"),12,0.5,10));
    cutlist.push_back(cutscanitem(cs.getCutPos("d3l0/d3El0"),12,1.0,16));
    cutlist.push_back(cutscanitem(cs.getCutPos("ptlb"),16,5.0,15.0));

    /*
    cutlist.push_back(cutscanitem(cs.getCutPos("iprpr"),12,1.1,4.7));
    cutlist.push_back(cutscanitem(cs.getCutPos("iprpi"),12,1.1,4.7));
    cutlist.push_back(cutscanitem(cs.getCutPos("ipr1m"),12,1.1,4.7));
    cutlist.push_back(cutscanitem(cs.getCutPos("ipr2m"),12,1.1,4.7));
    */
    //cutlist.push_back(cutscanitem(cs.getCutPos("probl0"),12,0.00,0.11));
    //cutlist.push_back(cutscanitem(cs.getCutPos("ctl0"),12,0.10,1.20));
    //cutlist.push_back(cutscanitem(cs.getCutPos("ctlb"),12,0.002,0.035));
    //cutlist.push_back(cutscanitem(cs.getCutPos("problb"),12,0.000,0.033));
    //cutlist.push_back(cutscanitem(cs.getCutPos("iprpr"),12,1.5,6.5));

    //cutlist.push_back(cutscanitem(cs.getCutPos("d3l0"),12,0.0,3.3));
    //cutlist.push_back(cutscanitem(cs.getCutPos("d3l0/d3El0"),12,3.0,36));
    //cutlist.push_back(cutscanitem(cs.getCutPos("d3lb/d3Elb"),12,0.0,33));

    //alte Cuts - pro memoria
    //cutlist.push_back(cutscanitem(2,12,0.05,0.07)); // mjp
    //cutlist.push_back(cutscanitem(3,12,0.06,0.08));// ml0
    //cutlist.push_back(cutscanitem(5,12,0.0016,0.0019));// alpha
    //cutlist.push_back(cutscanitem(6,12,2.3,2.7));// ptl0
    //cutlist.push_back(cutscanitem(7,12,2.4,2.8));// d3l0
    //cutlist.push_back(cutscanitem(8,12,27,31));// d3l0/d3El0
    //cutlist.push_back(cutscanitem(9,12,1.8,2.3));// d3lb/d3Elb

    if(fVerbose>-1) cout << "Current full cut string is:" << endl;
    if(fVerbose>-1) cout << cs.getCut(parvec) << endl;

    // Get TTree
    TTree* t = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;

    // define the fitfunction
    fitfuncGausBgrExp *fitfunc = new fitfuncGausBgrExp();
    //fitfuncGausBgrExpConst *fitfunc = new fitfuncGausBgrExpConst();
    //fitfunc->setInitValues(3,6,0,0,1,5.62,.1,2); // for Gaus + slope
    //fitfunc->setInitValues(5.23,6.03,14,-1.8,10,5.62,.02,2); // for Gaus + slope
    fitfunc->setInitValues(5.40,6.30,14,-1.8,10,5.62,.02,2); // for Gaus + slope, Lb, moved to right
    //fitfunc->setInitValues(5.0,6.4,120000,-1.95,1,5.62,.1,2); // for Gaus + Exp
    //fitfunc->setInitValues(5.1,6.4,120000,-1.5,1,5.62,.1,1,2); // for Gaus + Exp + Const

    // Set the histo for ALL fits
    const valueType histoXmin(4.5), histoXmax(6.5);
    TH1F *h = new TH1F("h1","h1",80,histoXmin,histoXmax);

    // limits for the mass summary plot
    const valueType histoSummaryXmin(5.5), histoSummaryXmax(5.7);

    // set up the canvas for the summary graphs
    unsigned int nGraphs(4); // no. of TGraphs per scan
    unsigned int nGrDrawn(0); // counter for the graphs drawn so far
    c1->Clear();
    c1->SetWindowSize(nGraphs*200,cutlist.size()*200);
    c1->SetCanvasSize(nGraphs*200,cutlist.size()*200);
    c1->Divide(nGraphs,cutlist.size());

    // Initialize printing.
    const std::string pdfNameHistos = "cutscanHistos";

    unsigned int canvascounter(0);
    // loop and do all cutscans in the list
    for(cutscanItemVecType::iterator it = cutlist.begin(); it!=cutlist.end(); it++)
    {
	c2->Clear();
	// actually do the cutscan
	cutscanresult csr = doScan(t,c2,h,fitfunc,plotstring,cs,parvec,(*it));
	c2->SaveAs((pdfNameHistos+toString(canvascounter)+".pdf").c_str());
	const std::string title = "Cutscan results for " + csr.cutname;
	// plot graphs
	if(csr.par.size()>0)
	{
	    const valueType stdCutVal = parvec[(*it).cutno_];
	    cout << "stdCutVal: " << stdCutVal << endl;
	    drawGraphCutline((TPad*)c1->cd(nGrDrawn+1), csr.par, csr.SoverSqrtSB, title, "Value for " + csr.cuttitle, "#frac{S}{#sqrt{S+B}}",stdCutVal);
	    drawGraphCutline((TPad*)c1->cd(nGrDrawn+2), csr.par, csr.sig, title, "Value for " + csr.cuttitle, "S",stdCutVal);
	    drawGraphCutline((TPad*)c1->cd(nGrDrawn+3), csr.par, csr.chi2, title, "Value for " + csr.cuttitle, "#chi^{2}",stdCutVal);
	    drawGraphErrors((TPad*)c1->cd(nGrDrawn+4), csr.par, csr.mean, csr.sigma, title, "Value for " + csr.cuttitle, "mean, i.e. mass",histoSummaryXmin,histoSummaryXmax);
	}
	nGrDrawn+=nGraphs;
	canvascounter++;
    }
    c1->SaveAs((pdfNameHistos+".pdf").c_str());
}

