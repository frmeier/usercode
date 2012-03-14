#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooAddModel.h"
#include "RooChebychev.h"
#include "RooTruthModel.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"

#include "utils.h"
#include "setTDRStyle_modified.C"

// Fits resolution function, lifetime and mass to MC "data"

using std::cout;
using std::endl;
using std::string;
using namespace RooFit;

int fVerbose(5);

// doResoTMassFitGenReco.C
int doResoTMassFitGenReco(TTree *tree, string title, bool isB0, double tLo = -5e-12, double tHi = 5e-12)
{
    setTDRStyle();
    const int nBins(100);

    // Observables in TTree
    RooRealVar tDiff("tDiff", "tDiff", -8e-12, 8e-12);
    RooRealVar t("t", "t", -5e-12, 20e-12);
    RooRealVar mass("mass", "mass", 4.6, 6.0);

    // Data
    RooDataSet data("data", "dataset", RooArgSet(tDiff, t, mass), Import(*tree));
    RooDataHist* binnedData = data.binnedClone();
    int curEntries = data.numEntries();
    if (fVerbose > 0) cout << "Dataset contains " << curEntries << " entries." << endl;

    RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries*1.4) ;
    RooRealVar nbgrpr("nbgrpr","prompt bgr fraction", .1*curEntries,0.,curEntries*1.0) ;

    // ================================================================================================================
    // Fitfunktion
    RooRealVar resoMean("resoMean", "resoMean", 0, tLo, tHi);
    RooRealVar resoSigma1("resoSigma1","resoSigma center", 0.15e-12, 0.001e-12, 0.5e-12);
    RooRealVar resoSigma2("resoSigma2","resoSigma tails",  0.30e-12, 0.1e-12,   1.7e-12);
    RooRealVar resoSigma3("resoSigma3","resoSigma tails",  2.0e-12, 0.5e-12,   50e-12);
    RooRealVar reso1frac("reso1frac","gauss1 fraction", .4, .01, .99);
    RooRealVar reso2frac("reso2frac","gauss1 fraction", .4, .01, .99);
    RooGaussModel resoGm1("resoGm1","gauss model center",tDiff,resoMean,resoSigma1) ;
    RooGaussModel resoGm2("resoGm2","gauss model tail",tDiff,resoMean,resoSigma2) ;
    RooGaussModel resoGm3("resoGm3","gauss model extended tail",tDiff,resoMean,resoSigma3) ;
    RooAddModel resoGm("resoGm","double gaussian",RooArgList(resoGm1,resoGm2,resoGm3),RooArgList(reso1frac,reso2frac));

    // do the fit
#define doNOTResoTMassFitGenReco_DOFIT
#ifdef doResoTMassFitGenReco_DOFIT
    resoSigma3.setConstant(true);
    resoGm.fitTo(*binnedData, NumCPU(4));
    if (fVerbose > 0) cout << "---------------------------------------------------------------------- this was with s3 const" << endl;
    resoSigma1.setConstant(true);
    resoSigma2.setConstant(true);
    resoSigma3.setConstant(false);
    resoGm.fitTo(*binnedData, NumCPU(4));
    if (fVerbose > 0) cout << "---------------------------------------------------------------------- this was with s1,s2 const" << endl;
    resoSigma1.setConstant(false);
    resoSigma2.setConstant(false);
    resoSigma3.setConstant(true);
    resoGm.fitTo(*binnedData, NumCPU(4));
    if (fVerbose > 0) cout << "---------------------------------------------------------------------- this was with s3 const" << endl;
    resoSigma1.setConstant(true);
    resoSigma2.setConstant(true);
    resoSigma3.setConstant(false);
    resoGm.fitTo(*binnedData, NumCPU(4));
    if (fVerbose > 0) cout << "---------------------------------------------------------------------- this was with s1,s2 const" << endl;
    resoSigma1.setConstant(false);
    resoSigma2.setConstant(false);
    resoSigma3.setConstant(true);
    resoGm.fitTo(*binnedData, NumCPU(4));
    if (fVerbose > 0) cout << "---------------------------------------------------------------------- this was with s3 const" << endl;
#else
    // NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
    // 1  reso1frac    4.90695e-01   9.16452e-03   1.26007e-04  -1.89905e-02
    // 2  reso2frac    3.63939e-01   8.98458e-03   1.55707e-04  -2.81373e-01
    // 3  resoMean     2.47458e-15   2.76896e-15   2.36725e-05  -5.23313e-01
    // 4  resoSigma1   1.60428e-13   3.41247e-15   1.65076e-04  -3.69352e-01
    // 5  resoSigma2   8.60333e-13   2.38555e-14   3.11409e-04  -4.96037e-02

    reso1frac.setVal(4.90695e-01); reso1frac.setError(9.16452e-03);
    reso2frac.setVal(3.63939e-01); reso2frac.setError(8.98458e-03);
    resoMean.setVal(2.47458e-15); resoMean.setError(2.76896e-15);
    resoSigma1.setVal(1.60428e-13); resoSigma1.setError(3.41247e-15);
    resoSigma2.setVal(8.60333e-13); resoSigma2.setError(2.38555e-14);
    resoSigma3.setVal(3.20110e-12); resoSigma3.setError(8.50472e-14);
#endif

    // Plotten
    RooPlot* frame = tDiff.frame(Title((title+" - Resolution function before full fit").c_str()));
    data.plotOn(frame, Binning(nBins));
    resoGm.plotOn(frame,ProjWData(data));
    resoGm.plotOn(frame,Components("resoGm1"),LineStyle(kDashed),LineColor(kGreen));
    resoGm.plotOn(frame,Components("resoGm2"),LineStyle(kDashed),LineColor(kCyan));
    resoGm.plotOn(frame,Components("resoGm3"),LineStyle(kDashed),LineColor(kMagenta));

    // write the results to the canvas
    {
	const double txtPosLeft   = .60;
	const double txtPosTop    = .80;
	const double txtLineSpace = .05;
	const double textsize     = .04;
	int i(0);

	frame->addObject(writeTLatex("mean: "             + roundToString(resoMean.getVal()*1e12  , 3) + " #pm " + roundToString(resoMean.getError()*1e12,   3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame->addObject(writeTLatex("center #sigma: "    + roundToString(resoSigma1.getVal()*1e12, 3) + " #pm " + roundToString(resoSigma1.getError()*1e12, 3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame->addObject(writeTLatex("tail #sigma: "      + roundToString(resoSigma2.getVal()*1e12, 3) + " #pm " + roundToString(resoSigma2.getError()*1e12, 3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame->addObject(writeTLatex("ext. tail #sigma: " + roundToString(resoSigma3.getVal()*1e12, 3) + " #pm " + roundToString(resoSigma3.getError()*1e12, 3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame->addObject(writeTLatex("(errors statistical only)", txtPosLeft, txtPosTop-4*txtLineSpace, textsize));
    }

    TCanvas *canvas = new TCanvas("fitres","fit results",1200,800);
    int canvascounter(0);
    canvas->Divide(3,2);
    canvas->cd(++canvascounter);
    frame->Draw();

    // Now we freeze some parameters of the resolution function
    resoMean.setConstant(true);
    resoSigma1.setConstant(false);
    resoSigma2.setConstant(false);
    resoSigma3.setConstant(true);
    reso1frac.setConstant(false);
    reso2frac.setConstant(false);

    // ================================================================================================================
    // Fit parameters
    //const double tauLb(1.391e-12), tauB0(1.525e-12);
    //const double tauLb(1.229e-12), tauB0(1.536e-12); // corresponds to the evt.pdl table used in production
    const double tauLb(1.424e-12), tauB0(1.536e-12); // corresponds to the evt.pdl table used in production
    const double tauTruthVal = (isB0 ? tauB0 : tauLb);
    RooRealVar tau("tau","tau", tauTruthVal, 0e-12, (tHi < 2*tauTruthVal ? 2*tauTruthVal : tHi));

    // model for signal:
    // - lifetime: decay using triple Gaussian as resolution model
    // - mass: double Gaussian

    // lifetime fit setup - we use values from fit above
    RooGaussModel resmodelGm1("resmodelGm1","gauss model center",t,resoMean,resoSigma1) ;
    RooGaussModel resmodelGm2("resmodelGm2","gauss model tail",t,resoMean,resoSigma2) ;
    RooGaussModel resmodelGm3("resmodelGm3","gauss model extended tail",t,resoMean,resoSigma3) ;
    RooAddModel resmodel("resmodel","multigaussian",RooArgList(resmodelGm1,resmodelGm2,resmodelGm3),RooArgList(reso1frac,reso2frac));

    RooDecay sig_t_model("sig_t_model", "signal model for time", t, tau, resmodel, RooDecay::SingleSided);

    // mass fit setup for signal
    RooRealVar sig_m_mean("sig_m_mean","sig_m_mean",5.28,5.26,5.30);
    RooRealVar sig_m_sigma1("sig_m_sigma1","sig_m_sigma1",0.010,0.005,0.020);
    RooRealVar sig_m_sigma2("sig_m_sigma2","sig_m_sigma2",0.020,0.012,0.040);
    RooRealVar sig_m_gauss1frac("sig_m_gauss1frac","gauss1 fraction", .3, .1, .9);

    RooGaussian sig_m_gauss1("sig_m_gauss1","signal Gaussian 1", mass, sig_m_mean, sig_m_sigma1);
    RooGaussian sig_m_gauss2("sig_m_gauss2","signal Gaussian 2", mass, sig_m_mean, sig_m_sigma2);
    RooAddPdf sig_m_model("sig_m_model","double gaussian", RooArgList(sig_m_gauss1, sig_m_gauss2), RooArgList(sig_m_gauss1frac));

    RooProdPdf sig_model("sig_model", "signal model", RooArgSet(sig_t_model, sig_m_model));

    // model for prompt background (mainly combinatorics)
    RooRealVar bgrpr_m_mean("bgrpr_m_mean", "prompt bgr mass mean", 5.28,5.0, 5.6);
    RooRealVar bgrpr_m_sigma("bgrpr_m_sigma", "prompt bgr mass sigma", 2.0, 1.0, 8.0);
    RooGaussian bgrpr_m("bgrpr_m", "prompt bgr mass Gaussian", mass, bgrpr_m_mean, bgrpr_m_sigma);

    RooRealVar bgrpr_t_mean("bgrpr_t_mean", "prompt bgr time mean", 0, -2e-12, 2e-12);
    RooRealVar bgrpr_t_sigma("bgrpr_t_sigma", "prompt bgr time sigma", 1e-12, 0.1e-12, 5e-12);
    RooGaussian bgrpr_t("bgrpr_t", "prompt bgr mass Gaussian", mass, bgrpr_t_mean, bgrpr_t_sigma);

    RooProdPdf bgrpr_model("bgrpr_model", "model for prompt bgr", RooArgList(bgrpr_m, bgrpr_t));

    // final model
    RooAddPdf model("model", "model", RooArgSet(sig_model, bgrpr_model), RooArgList(nsig, nbgrpr));

    // now select what we want to float or fix
    resoMean.setConstant(true);
    resoSigma1.setConstant(true);
    resoSigma2.setConstant(true);
    resoSigma3.setConstant(true);
    reso1frac.setConstant(true);
    reso2frac.setConstant(true);

    sig_m_mean.setConstant(false);
    sig_m_sigma1.setConstant(false);
    sig_m_sigma2.setConstant(false);
    sig_m_gauss1frac.setConstant(false);

    bgrpr_m_mean.setConstant(true);
    bgrpr_m_sigma.setConstant(true);
    bgrpr_t_mean.setConstant(true);
    bgrpr_t_sigma.setConstant(true);

    tau.setConstant(true);

    // do the fit
    model.fitTo(data, NumCPU(4));
    if (fVerbose > 0) cout << "model.fitTo(data) finished. -------------------------------" << endl;

    // some plotting
    canvas->cd(++canvascounter);
    RooPlot * frame_t = t.frame(Title((title+" - Lifetime").c_str()));
    data.plotOn(frame_t, Binning(nBins));
    if (fVerbose > 0) cout << "data.plotOn finished. -------------------------------" << endl;
    //modelTruth.plotOn(frame_t,ProjWData(data),LineStyle(kDashed),LineColor(kGreen));
    model.plotOn(frame_t,ProjWData(data));
    //model.plotOn(frame_t);
    if (fVerbose > 0) cout << "model.plotOn finished. -------------------------------" << endl;
    //model.plotOn(frame_t,Components("t"),LineStyle(kDashed),LineColor(kGreen));

    // write the results to the canvas
    {
	const double txtPosLeft   = .60;
	const double txtPosTop    = .80;
	const double txtLineSpace = .05;
	const double textsize     = .04;
	int i(0);

	frame_t->addObject(writeTLatex("#tau: " + roundToString(tau.getVal()*1e12, 3) + " #pm " + roundToString(tau.getError()*1e12, 3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame_t->addObject(writeTLatex("#tau expected: " + roundToString(tauTruthVal*1e12, 3)  + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame_t->addObject(writeTLatex("this is " + roundToString(tau.getVal()/tauTruthVal*100, 1) + " #pm " + roundToString(tau.getError()/tauTruthVal*100, 1) + " %", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame_t->addObject(writeTLatex("n Signal: " + roundToString(nsig.getVal(), 0) + " #pm " + roundToString(nsig.getError(), 0), txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame_t->addObject(writeTLatex("(errors statistical only)", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
    }
    frame_t->Draw();
    gPad->SetLogy();

    // plot mass
    canvas->cd(++canvascounter);
    RooPlot * frame_m = mass.frame(Title((title+" - Mass").c_str()));
    data.plotOn(frame_m, Binning(nBins));
    model.plotOn(frame_m,ProjWData(data));
    if (fVerbose > 0) cout << "mass plot finished. -------------------------------" << endl;
    frame_m->Draw();

    // draw fit of resofunction again
    // Plotten
    RooPlot* frame2 = tDiff.frame(Title((title+" - Resolution function after full fit").c_str()));
    data.plotOn(frame2, Binning(nBins));
    resoGm.plotOn(frame2,ProjWData(data));
    resoGm.plotOn(frame2,Components("resoGm1"),LineStyle(kDashed),LineColor(kGreen));
    resoGm.plotOn(frame2,Components("resoGm2"),LineStyle(kDashed),LineColor(kCyan));
    resoGm.plotOn(frame2,Components("resoGm3"),LineStyle(kDashed),LineColor(kMagenta));

    // write the results to the canvas
    {
	const double txtPosLeft   = .60;
	const double txtPosTop    = .80;
	const double txtLineSpace = .05;
	const double textsize     = .04;
	int i(0);

	frame2->addObject(writeTLatex("mean: "             + roundToString(resoMean.getVal()*1e12  , 3) + " #pm " + roundToString(resoMean.getError()*1e12,   3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame2->addObject(writeTLatex("center #sigma: "    + roundToString(resoSigma1.getVal()*1e12, 3) + " #pm " + roundToString(resoSigma1.getError()*1e12, 3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame2->addObject(writeTLatex("tail #sigma: "      + roundToString(resoSigma2.getVal()*1e12, 3) + " #pm " + roundToString(resoSigma2.getError()*1e12, 3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame2->addObject(writeTLatex("ext. tail #sigma: " + roundToString(resoSigma3.getVal()*1e12, 3) + " #pm " + roundToString(resoSigma3.getError()*1e12, 3) + " ps", txtPosLeft, txtPosTop-(i++)*txtLineSpace, textsize));
	frame2->addObject(writeTLatex("(errors statistical only)", txtPosLeft, txtPosTop-4*txtLineSpace, textsize));
    }

    canvas->cd(++canvascounter);
    frame2->Draw();
    // we reached the end without an error
    return 0;
}

int doResoTMassFitGenReco(string filename, string title, bool isB0, double tLo = -5e-12, double tHi = 15e-12)
{
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
	cout << "This script requires a very reduced tree called fittree" << endl;
	return -2;
    }
    return doResoTMassFitGenReco(tree, title, isB0, tLo, tHi);
}

