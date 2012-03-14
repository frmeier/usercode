#include <string>
#include <iostream>

#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"

#include "utils.h"
#include "setTDRStyle_modified.C"

using std::string;
using std::cout;
using std::endl;

Double_t computeIntegralRatio(TH1F *h)
{
    Int_t zeroBin = 0;
    for (zeroBin=1; zeroBin!=h->GetNbinsX(); zeroBin++)
	if (h->GetBinCenter(zeroBin) > 0) break;
    zeroBin--;
    Double_t left = h->Integral(1,zeroBin);
    Double_t right = h->Integral(zeroBin+1,h->GetNbinsX());
    return right/left;
}

void quickTsigBgrPlot(string filename, bool isB0)
{
    setTDRStyle();

    TFile *_file0 = TFile::Open(filename.c_str());

    TTree *tree = (TTree*)gDirectory->Get("fittree");

    TCanvas *c = new TCanvas("c","c",600,500);
    int canvasctr(0);

    Double_t sigLo, sigHi, bgrLo, nBgr;
    if (isB0)
    {
	sigLo = 5.24897; sigHi = 5.30833; bgrLo = 5.32323; nBgr = 710.5;
    }
    else
    {
	sigLo = 5.57042; sigHi = 5.66637; bgrLo = 5.47359; nBgr = 685.264;
    }
    string sigCut = "mass>" + toString(sigLo) + "&&mass<" + toString(sigHi);
    string bgrCut = "mass>" + toString(bgrLo);
    string hstring = "(80,-5e-12,15e-12)";
    tree->Draw(("t>>hsig"+hstring).c_str(),sigCut.c_str());
    tree->Draw(("t>>hbgr"+hstring).c_str(),bgrCut.c_str(),"same");
    tree->Draw(("-t>>hbgri"+hstring).c_str(),bgrCut.c_str(),"same");

    TH1F *hsig = (TH1F*)gDirectory->GetList()->FindObject("hsig");
    TH1F *hbgr = (TH1F*)gDirectory->GetList()->FindObject("hbgr");
    TH1F *hbgri= (TH1F*)gDirectory->GetList()->FindObject("hbgri");

    hsig->SetLineColor(8);
    hsig->SetLineWidth(2);
    hbgr->SetLineColor(9);
    hbgr->SetLineWidth(2);
    hbgri->SetLineColor(9);
    hbgri->SetLineStyle(2);
    hbgri->SetLineWidth(2);

    /*
    cout << "Mean sig: " << hsig->GetMean(1) << " bgr: " << hbgr->GetMean(1) << endl;
    cout << "RMS sig: " << hsig->GetRMS(1) << " bgr: " << hbgr->GetRMS(1) << endl;
    cout << "Skewness sig: " << hsig->GetSkewness(1) << " bgr: " << hbgr->GetSkewness(1) << endl;
    cout << "Maximum bin: " << hsig->GetMaximumBin() << " bgr: " << hbgr->GetMaximumBin() << endl;
    cout << "Bin center of max: " << hsig->GetBinCenter(hsig->GetMaximumBin()) << " bgr: " << hbgr->GetBinCenter(hbgr->GetMaximumBin()) << endl;
    cout << "Pearson mode: " << (hsig->GetMean()-hsig->GetBinCenter(hsig->GetMaximumBin()))/hsig->GetRMS()
         << " bgr: " << (hbgr->GetMean()-hbgr->GetBinCenter(hbgr->GetMaximumBin()))/hbgr->GetRMS() << endl;
    */
    cout << "right/left ratio around zero:" << endl;
    cout << "hsig: " << computeIntegralRatio(hsig) << endl;
    cout << "hbgr: " << computeIntegralRatio(hbgr) << endl;
    cout << "hbgr flipped: " << computeIntegralRatio(hbgri) << endl;

    hbgr->Scale(nBgr/hbgr->GetEntries());
    hbgri->Scale(nBgr/hbgr->GetEntries());

    hsig->SetTitle(isB0 ? "B^{0}" : "#Lambda_{b}");
    hsig->GetXaxis()->SetTitle("Lifetime (ps)");
    hsig->GetYaxis()->SetTitle("Entries");
    gPad->SetLogy();
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);
}

