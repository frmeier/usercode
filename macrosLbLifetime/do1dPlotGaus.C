#include "TROOT.h"
#include "TTree.h"
#include "TH1F.h"
#include "TPad.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLatex.h"
#include "utils.h"
#include <string>
//#include "Cuts.C"

using std::string;

void do1dPlotGaus(TTree *t, string hname, string todraw, string cut,
	int nBinsX, double minX, double maxX, double minY, double maxY,
	string title, string titleX, string unitX)
{
    gStyle->SetPalette(1);
    const std::string plotstring = todraw + ">>" + hname + "(" +
	toString(nBinsX) + "," + toString(minX) + "," + toString(maxX) + ")";
    cout << plotstring << endl;
    t->Draw(plotstring.c_str(), cut.c_str());
    TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(hname.c_str());
    gPad->SetLeftMargin(.15);
    gPad->SetRightMargin(.1);
    gPad->SetTopMargin(.08);
    gPad->SetBottomMargin(.14);
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h->GetXaxis()->SetNdivisions(505);
    const string titleY = "Entries per bin (width: " + toString(h->GetBinWidth(1)) + (unitX.size()>0 ? unitX : "") + ")";
    h->GetYaxis()->SetTitle(titleY.c_str());
    h->SetMinimum(minY);
    h->SetMaximum(maxY);

    gStyle->SetOptFit(0);
    h->Fit("gaus");
    TF1 *fit = h->GetFunction("gaus");
    const double xpos = 0.6;
    const double ypos = 0.7;
    const double yspacing = 0.05;
    int yline = 0;
    TLatex *tl0 = writeTLatex("Gaussian fit:", xpos, ypos-yspacing*yline++);
    TLatex *tl1 = writeTLatex("mean: "+roundToString(fit->GetParameter(1),3)+" #pm "+roundToString(fit->GetParError(1),3)+" "+unitX, xpos, ypos-yspacing*yline++);
    TLatex *tl2 = writeTLatex("#sigma: "+roundToString(fit->GetParameter(2),3)+" #pm "+roundToString(fit->GetParError(2),3)+" "+unitX, xpos, ypos-yspacing*yline++);
    tl0->Draw();
    tl1->Draw();
    tl2->Draw();
    gPad->Update();
    return; 
}

void do1dPlotGaus(TTree *t, string hname, string todraw, const Cuts &cut,
	int nBinsX, double minX, double maxX, double minY, double maxY,
	string title, string titleX, string unitX)
{
    do1dPlotGaus(t, hname, todraw, cut.getCut(), nBinsX, minX, maxX, minY, maxY, title, titleX, unitX);
}

void do1dPlotGaus(TTree *t, string hname, string todraw, const Cuts &cut,
	int nBinsX, double minX, double maxX,
	string title, string titleX, string unitX)
{
    do1dPlotGaus(t, hname, todraw, cut.getCut(), nBinsX, minX, maxX, -1111, -1111, title, titleX, unitX);
}

void do1dPlotGaus(TTree *t, string hname, string todraw, const string cut,
	int nBinsX, double minX, double maxX,
	string title, string titleX, string unitX)
{
    do1dPlotGaus(t, hname, todraw, cut, nBinsX, minX, maxX, -1111, -1111, title, titleX, unitX);
}

