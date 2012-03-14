#include "TROOT.h"
#include "TTree.h"
#include "TH2F.h"
#include "TPad.h"
#include "TGraph.h"
#include "TStyle.h"
#include "utils.h"
#include <string>
#include "Cuts.C"

using std::string;

void do2dPlot(TTree *t, string hname, string todraw, string cut,
	int nBinsX, double minX, double maxX, int nBinsY, double minY, double maxY,
	string title, string titleX, string unitX, string titleY, string unitY)
{
    gStyle->SetPalette(1);
    const std::string plotstring = todraw + ">>" + hname + "(" +
	toString(nBinsX) + "," + toString(minX) + "," + toString(maxX) + "," +
	toString(nBinsY) + "," + toString(minY) + "," + toString(maxY) + ")";
    t->Draw(plotstring.c_str(), cut.c_str(),"COLZ");
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(hname.c_str());
    gPad->SetLeftMargin(.15);
    gPad->SetRightMargin(.2);
    gPad->SetTopMargin(.1);
    gPad->SetBottomMargin(.12);
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitle(unitY.size()>0 ? (titleY+" / "+unitY).c_str() : titleY.c_str());
    return; 
}

void do2dPlot(TTree *t, string hname, string todraw, const Cuts &cut,
	int nBinsX, double minX, double maxX, int nBinsY, double minY, double maxY,
	string title, string titleX, string unitX, string titleY, string unitY)
{
    do2dPlot(t, hname, todraw, cut.getCut(), nBinsX, minX, maxX, nBinsY, minY, maxY, title, titleX, unitX, titleY, unitY);
}

// a version without bins
void do2dPlot(TTree *t, string hname, string todraw, string cut,
	double minX, double maxX, double minY, double maxY,
	string title, string titleX, string unitX, string titleY, string unitY)
{   // TODO for some strange reasons the titles are not adjusted...
    (void)minX;
    (void)maxX; // this plot is not mission critical and I was unable to find out how to set the xrange... Suppresses the annoying compiler warning of unused variables
    gStyle->SetPalette(1);
    t->Draw(todraw.c_str(), cut.c_str());
    TGraph *gr = (TGraph*)gPad->GetPrimitive("Graph");
    gr->SetName(hname.c_str());
    gr->SetTitle(title.c_str());
    gr->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    gr->GetYaxis()->SetTitle(unitY.size()>0 ? (titleY+" / "+unitY).c_str() : titleY.c_str());
    gr->SetMinimum(minY);
    gr->SetMaximum(maxY);
    gr->Draw("P");
    return; 
}

