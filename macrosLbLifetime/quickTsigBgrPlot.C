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

void quickTsigBgrPlot(string filename, bool isB0, string title = "", string saveAs = "")
{
    setTDRStyle();

    TFile *_file0 = TFile::Open(filename.c_str());

    TTree *tree = (TTree*)gDirectory->Get("fittree");

    TCanvas *c = new TCanvas("c","c",600,500);
    int canvasctr(0);

    Double_t sigLo, sigHi, bgrLo, nSig, nBgr;
    if (isB0)
    {
	//sigLo = 5.24897; sigHi = 5.30833; bgrLo = 5.32323; nBgr = 710.5;
	sigLo = 5.24976; sigHi = 5.30849; bgrLo = 5.32571; nSig = 5697; nBgr = 615.4;
	// 2σ range: 5.24976-5.30849 3σ range: 5.23254-5.32571
    }
    else
    {
	// sigLo = 5.57042; sigHi = 5.66637; bgrLo = 5.47359; nBgr = 685.264; old
	sigLo = 5.57884; sigHi = 5.65855; bgrLo = 5.68451; nSig = 752.2; nBgr = 332.4;
	// 2σ range: 5.57884-5.65855 3σ range: 5.55288-5.68451
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
    hsig->SetFillColor(0);
    hbgr->SetLineColor(9);
    hbgr->SetLineWidth(2);
    hbgr->SetFillColor(0);
    hbgri->SetLineColor(9);
    hbgri->SetLineStyle(2);
    hbgri->SetLineWidth(2);
    hbgri->SetFillColor(0);

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

    //const double scale = nBgr/hbgr->GetEntries(); // old
    const double scale = hsig->GetEntries()/hbgr->GetEntries()*nBgr/(nBgr+nSig); // same as in sidebandsubstracted plots
    hbgr->Scale(scale);
    hbgri->Scale(scale);

    if (title.size() == 0)
	hsig->SetTitle(isB0 ? "B^{0}" : "#Lambda_{b}");
    else
	hsig->SetTitle(title.c_str());
    hsig->GetXaxis()->SetTitle("Lifetime (ps)");
    hsig->GetYaxis()->SetTitle("Entries");
    gPad->SetLogy();
    gPad->SetRightMargin(0.1);
    gPad->SetTopMargin(0.1);

    if (saveAs.size() != 0) c->SaveAs(saveAs.c_str());
}

void doSomePlots()
{
    quickTsigBgrPlot("../data/vrt_r333_lb_data_barrelMatch.root", false, "#Lambda_{b} in data r333", "r333_quickTsigBgrPlot_lb_data.pdf");
    quickTsigBgrPlot("../data/vrt_r332_lb_MC_barrelMatch.root", false, "#Lambda_{b} in MC r332", "r332_quickTsigBgrPlot_lb_mc.pdf");
    quickTsigBgrPlot("../data/vrt_r335_B0_data_barrelMatch.root", true, "B^{0} in data", "r335_quickTsigBgrPlot_B0_data.pdf");
    quickTsigBgrPlot("../data/vrt_r334_B0_MC_barrelMatch.root", true, "B^{0} in MC", "r334_quickTsigBgrPlot_B0_mc.pdf");
}

/*
../data/vrt_r306_bgr_JPsiToMuMu.root
../data/vrt_r307_bgr_BsToPsiMuMu.root
../data/vrt_r308_bgr_BpToPsiMuMu.root
../data/vrt_r309_bgr_B0ToPsiMuMu.root
../data/vrt_r310_bgr_XibToJpsiXi.root
../data/vrt_r311_bgr_OmegabToJpsiOmega.root
../data/vrt_r313_lb_data.root
../data/vrt_r314_lb_mc.root
../data/vrt_r314_lb_mc_match.root
../data/vrt_r315_B0_data.root
../data/vrt_r316_B0_mc_match.root
*/

