#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"


using std::cout;
using std::endl;
using std::string;
using namespace RooFit;

int doCtFitTest(string filename, bool isB0 = false)
{
    const int nBins(50);

    // Observables in TTree
    //RooRealVar t("t","t",-5e-12,20e-12);
    RooRealVar t("t","t",-4e-12,20e-12);
    RooRealVar tE("tE","t error", -10e-12, 10e-12);
    const double massLb(5.620), massLbWindow(0.2);
    const double massB0(5.27950), massB0Window(0.2);
    const double mass(isB0 ? 5.27950 : 5.620);
    const double massWindow = 0.2;
    RooRealVar m("mass","m", mass-massWindow, mass+massWindow);
    //RooRealVar m("mass","m", mass-0.05, mass+massWindow);

    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
	return -1;
    }
    TTree *tree = (TTree*)f->Get("fittree");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
	return -2;
    }
    RooDataSet data("data", "dataset", RooArgSet(t, tE, m), Import(*tree));
    int curEntries = data.numEntries();
    cout << "Dataset contains " << curEntries << " entries." << endl;

    // Fit parameters
    const double tauLb(1.391e-12), tauB0(1.525e-12);
    RooRealVar tau("tau","tau", (isB0 ? tauB0 : tauLb), 1e-12, 2e-12);

    RooRealVar tEbias("tEbias","tEbias",0,-1e-13,1e-12);
    RooRealVar tEsigma("tEsigma","per-event error scale factor",1,0.1,10);

    RooRealVar mMean("mMean","mass mean", mass, mass-massWindow, mass+massWindow);
    RooRealVar mSigma("mSigma","mass sigma", 0.014, 0.01, 0.03);
    // Background
    RooRealVar tEbias_bgr("tEbias_bgr","tEbias_bgr",0,-1e-13,1e-12);
    RooRealVar tEsigma_bgr("tEsigma_bgr","per-event error scale factor",1,0.1,10);
    RooRealVar tau_bgr("tau_bgr","tau_bgr", 1.5e-12, .1e-12, 10e-12);
    RooRealVar c0_bgr("c0_bgr","coefficient c0 bgr", 1.0,-1,2) ;
    RooRealVar c1_bgr("c1_bgr","coefficient c1 bgr", 0.1,-1,1) ;
    RooChebychev mBgr("mBgr","background mass",m,RooArgSet(c0_bgr)) ;

    RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries*1.4) ;
    RooRealVar nbkg("nbkg","Background fraction", .5* curEntries,0.,curEntries*1.4) ;

    // set some parameters constant
    tEbias.setConstant(true);
    tEbias_bgr.setConstant(true);

    // Model building
    // -- decay
    RooGaussModel gm("gm1","gauss model scaled by per-event error",t,tEbias,tEsigma,tE);
    RooDecay decay_gm("decay_gm", "decay",t,tau,gm,RooDecay::SingleSided);
    // -- mass fit
    RooGaussian mGauss("mGauss", "Gauss model for signal", m, mMean, mSigma);
    // combined mass and lifetime
    RooProdPdf massTimeModel("masstimemodel", "Model for Mass and Lifetime", RooArgSet(decay_gm, mGauss));

    // now background: flat plus decay
    RooGaussModel gm_bgr("gm_bgr","gauss model scaled by per-event error",t,tEbias_bgr,tEsigma_bgr,tE);
    RooDecay decay_gm_bgr("decay_gm_bgr", "decay background",t,tau_bgr,gm_bgr,RooDecay::SingleSided);
    RooProdPdf bgrModel("bgrmodel", "Model for Background", RooArgSet(decay_gm_bgr, mBgr));

    // final model
    RooAddPdf model("model", "model", RooArgSet(massTimeModel, bgrModel), RooArgList(nsig,nbkg));

    // do the fit
    model.fitTo(data,ConditionalObservables(tE));

    // some plotting
    RooPlot * frame_t = t.frame(Title("Lifetime"));
    data.plotOn(frame_t, Binning(nBins));
    model.plotOn(frame_t,ProjWData(data));

    RooPlot * frame_m = m.frame(Title("Mass"));
    data.plotOn(frame_m, Binning(nBins));
    model.plotOn(frame_m,ProjWData(data));

    // draw the plots
    TCanvas *canvas = new TCanvas("fitres","fit results",800,800);
    canvas->Divide(1,2);
    canvas->cd(1);
    frame_t->Draw();
    gPad->SetLogy();
    canvas->cd(2);
    frame_m->Draw();
    gPad->SetLogy(false);

    // we reached the end without an error
    return 0;
}

