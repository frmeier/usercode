#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPaletteAxis.h"

#include "utils.h"
#include "setTDRStyle_modified.C"

using std::string;

void doMassT3dPlot(TTree *tree, string title = "default title", double lumiPb = 0, bool prelim = true)
{
    setTDRStyle();
    TCanvas *c1 = new TCanvas("mtplane","mass lifetime plot", 800,600);
    gStyle->SetPalette(1);
    // define plot range
    const double massLo(4.6), massHi(6.3);
    const double tLo(-5e-12), tHi(25e-12);
    string cut = "mass>" + toString(massLo) + "&&mass<" + toString(massHi) + "&&t>" + toString(tLo) + "&&t<" + toString(tHi);
    // define binning
    const int nBinsMass(34), nBinsT(20);
    const string histoName = "histo";
    string toDraw = "mass:t>>" + histoName +
	"(" + toString(nBinsT) + "," + toString(tLo) + "," + toString(tHi) +
	"," + toString(nBinsMass) + "," + toString(massLo) + "," + toString(massHi) + ")";
    // draw histo
    tree->Draw(toDraw.c_str(),cut.c_str(),"surf3z");
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(histoName.c_str());
    // now make it nicer
    h->SetTitle(title.c_str());
    h->GetYaxis()->SetTitle("mass [GeV/c^{2}]");
    h->GetYaxis()->SetTitleOffset(1.38);
    h->GetYaxis()->CenterTitle(true);
    h->GetXaxis()->SetTitle("lifetime [ps]");
    h->GetXaxis()->SetTitleOffset(1.38);
    h->GetXaxis()->CenterTitle(true);
    gPad->SetPhi(240);
    gPad->SetTheta(15);
    gPad->SetLogz();
    gPad->SetRightMargin(0.2);
    c1->Update();
    TPaletteAxis *pal = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
    pal->SetX1NDC(0.86); pal->SetY1NDC(0.24); pal->SetX2NDC(0.90); pal->SetY2NDC(0.70);
    pal->Draw();
    c1->Update();
    // draw lumi if provided
    if (lumiPb >0)
    {
	TLatex * l = new TLatex;
	l->SetTextSize(0.04);
	l->SetNDC(true);
	double lumi;
	string lumiUnit;
	if (lumiPb > 1100)
	{
	    lumi = lumiPb/1000;
	    lumiUnit = " fb^{-1}";
	}
	else
	{
	    lumi = lumiPb;
	    lumiUnit = " pb^{-1}";
	}
	string lumistring = roundToString(lumi,2-int(TMath::Log10(lumi))) + lumiUnit;
	l->DrawLatex(0.80,0.05,("#int Ldt = " + lumistring).c_str());
    }
    if (prelim)
    {
	TLatex * l = new TLatex;
	l->SetTextSize(0.05);
	l->SetNDC(true);
	l->DrawLatex(0.65,0.95,"Work in progress");
    }
}

void doMassT3dPlot(string filename, string title = "default title", double lumiPb = 0, bool prelim = true)
{
    TFile *_file0 = TFile::Open(filename.c_str());
    TTree *tree = (TTree*)gDirectory->Get("fittree");
    doMassT3dPlot(tree, title, lumiPb, prelim);
}

void test(double l)
{
    TFile *_file0 = TFile::Open("veryreducedtrees/fit_B0.root");
    TTree *tree = (TTree*)gDirectory->Get("fittree");
    doMassT3dPlot(tree, "Titel", l);
}

