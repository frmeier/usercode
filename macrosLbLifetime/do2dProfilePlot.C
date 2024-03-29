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
    const double maxscaler(.1);
    h->SetMaximum(maxY*maxscaler);
    h->SetMinimum(minY*maxscaler);
    return; 
}

void do2dProfilePlotVariableBins(TTree *t, string hname, string todrawX, string todrawY, string cut,
	int nBinsX, double minX, double maxX, double minY, double maxY,
	string title, string titleX, string unitX, string titleY, string unitY, double scaleX = 1.0, double scaleY = 1.0)
{
    std::vector<double> bins = makeDynamicBins(t, todrawX, cut, nBinsX, minX, maxX, scaleX);
    TProfile *h = new TProfile(hname.c_str(), title.c_str(), nBinsX, &bins[0]);
    //const std::string plotstring = todrawY+(scaleY!=1.0?"*"+toString(scaleY):"") +":"+ todrawX+(scaleX!=1.0?"*"+toString(scaleX):"") + ">>" + hname;
    std::string plotstring;
    if (scaleY!=1.0)
	plotstring += "("+todrawY+")*"+toString(scaleY);
    else
	plotstring += todrawY;
    plotstring += ":";
    if (scaleX!=1.0)
	plotstring += "("+todrawX+")*"+toString(scaleX);
    else
	plotstring += todrawX;
    plotstring += ">>" + hname;

    cout << "plotstring: " << plotstring << endl;
    t->Draw(plotstring.c_str(), cut.c_str());
    gPad->SetLeftMargin(.15);
    gPad->SetRightMargin(.1);
    gPad->SetTopMargin(.1);
    gPad->SetBottomMargin(.12);
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(valueWithUnit(titleX, unitX).c_str());
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitle(valueWithUnit(titleY, unitY).c_str());
    const double maxscaler(.1);
    h->SetMaximum(maxY*maxscaler);
    h->SetMinimum(minY*maxscaler);
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

