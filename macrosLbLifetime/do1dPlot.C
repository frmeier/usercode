#include "TROOT.h"
#include "TTree.h"
#include "TH1F.h"
#include "TPad.h"
#include "TGraph.h"
#include "TStyle.h"
#include "utils.h"
#include <string>
//#include "Cuts.C"

using std::string;

void do1dPlot(TTree *t, string hname, string todraw, string cut,
	int nBinsX, double minX, double maxX, double minY, double maxY,
	string title, string titleX, string unitX)
{
    gStyle->SetPalette(1);
    const std::string plotstring = todraw + ">>" + hname + "(" +
	toString(nBinsX) + "," + toString(minX) + "," + toString(maxX) + ")";
    t->Draw(plotstring.c_str(), cut.c_str());
    TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(hname.c_str());
    gPad->SetLeftMargin(.15);
    gPad->SetRightMargin(.2);
    gPad->SetTopMargin(.1);
    gPad->SetBottomMargin(.12);
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h->GetXaxis()->SetNdivisions(505);
    const string titleY = "Entries per bin (width: " + toString(h->GetBinWidth(1)) + (unitX.size()>0 ? unitX : "") + ")";
    h->GetYaxis()->SetTitle(titleY.c_str());
    h->SetMinimum(minY);
    h->SetMaximum(maxY);
    return; 
}

void do1dPlot(TTree *t, string hname, string todraw, const Cuts &cut,
	int nBinsX, double minX, double maxX, double minY, double maxY,
	string title, string titleX, string unitX)
{
    do1dPlot(t, hname, todraw, cut.getCut(), nBinsX, minX, maxX, minY, maxY, title, titleX, unitX);
}

void do1dPlot(TTree *t, string hname, string todraw, const Cuts &cut,
	int nBinsX, double minX, double maxX,
	string title, string titleX, string unitX)
{
    do1dPlot(t, hname, todraw, cut.getCut(), nBinsX, minX, maxX, -1111, -1111, title, titleX, unitX);
}

void do1dPlot(TTree *t, string hname, string todraw, const string cut,
	int nBinsX, double minX, double maxX,
	string title, string titleX, string unitX)
{
    do1dPlot(t, hname, todraw, cut, nBinsX, minX, maxX, -1111, -1111, title, titleX, unitX);
}

