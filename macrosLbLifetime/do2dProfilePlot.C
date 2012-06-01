#include "TROOT.h"
#include "TTree.h"
#include "TProfile.h"
#include "TPad.h"
#include "TGraph.h"
#include "TStyle.h"
#include "utils.h"
#include <string>
//#include "Cuts.C"

using std::string;

void do2dProfilePlot(TTree *t, string hname, string todraw, string cut,
	int nBinsX, double minX, double maxX, double minY, double maxY,
	string title, string titleX, string unitX, string titleY, string unitY)
{
    TProfile *h = new TProfile(hname.c_str(), title.c_str(), nBinsX, minX, maxX);
    const std::string plotstring = todraw + ">>" + hname;
    t->Draw(plotstring.c_str(), cut.c_str());
    gPad->SetLeftMargin(.15);
    gPad->SetRightMargin(.1);
    gPad->SetTopMargin(.1);
    gPad->SetBottomMargin(.12);
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(valueWithUnit(titleX, unitX).c_str());
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitle(valueWithUnit(titleY, unitY).c_str());
    h->SetMaximum(maxY);
    h->SetMinimum(minY);
    return; 
}

void do2dProfilePlot(TTree *t, string hname, string todraw, string cut,
	int nBinsX, double minX, double maxX,
	string title, string titleX, string unitX, string titleY, string unitY)
{
    do2dProfilePlot(t, hname, todraw, cut, nBinsX, minX, maxX, -1111, -1111, title, titleX, unitX, titleY, unitY);
}

void do2dProfilePlot(TTree *t, string hname, string todraw, const Cuts &cut,
	int nBinsX, double minX, double maxX,
	string title, string titleX, string unitX, string titleY, string unitY)
{
    do2dProfilePlot(t, hname, todraw, cut.getCut(), nBinsX, minX, maxX, -1111, -1111, title, titleX, unitX, titleY, unitY);
}

void do2dProfilePlot(TTree *t, string hname, string todraw, const Cuts &cut,
	int nBinsX, double minX, double maxX, double minY, double maxY,
	string title, string titleX, string unitX, string titleY, string unitY)
{
    do2dProfilePlot(t, hname, todraw, cut.getCut(), nBinsX, minX, maxX, minY, maxY, title, titleX, unitX, titleY, unitY);
}

