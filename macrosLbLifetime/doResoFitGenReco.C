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


using std::cout;
using std::endl;
using std::string;
using namespace RooFit;

int doResoFitGenReco(TTree *tree, string title, bool isB0, double tLo = -5e-12, double tHi = 5e-12)
{
    setTDRStyle();
    const int nBins(100);

    // Observables in TTree
    RooRealVar tDiff("tDiff", "tDiff", -8e-12, 8e-12);
    RooRealVar t("t", "t", -5e-12, 20e-12);

    // Data
    RooDataSet data("data", "dataset", RooArgSet(tDiff,t), Import(*tree));
    RooDataHist* binnedData = data.binnedClone();
    int curEntries = data.numEntries();
    cout << "Dataset contains " << curEntries << " entries." << endl;

    RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries*1.4) ;

    // Fitfunktion
    RooRealVar resoMean("resoMean", "resoMean", 0, tLo, tHi);
    RooRealVar resoSigma1("resoSigma1","resoSigma center", 0.15e-12, 0.001e-12, 0.5e-12);
    RooRealVar resoSigma2("resoSigma2","resoSigma tails",  0.30e-12, 0.1e-12,   1.7e-12);
    RooRealVar resoSigma3("resoSigma3","resoSigma tails",  2.0e-12, 0.5e-12,   50e-12);
    RooRealVar reso1frac("reso1frac","gauss1 fraction", .4, .01, .99);
    RooRealVar reso2frac("reso2frac","gauss1 fraction", .4, .01, .99);
    RooRealVar reso3frac("reso3frac","gauss1 fraction", .2, .01, .99);
    RooGaussModel resoGm1("resoGm1","gauss model center",tDiff,resoMean,resoSigma1) ;
    RooGaussModel resoGm2("resoGm2","gauss model tail",tDiff,resoMean,resoSigma2) ;
    RooGaussModel resoGm3("resoGm3","gauss model extended tail",tDiff,resoMean,resoSigma3) ;
    //RooRealVar frac_core("frac_core","core fraction",0.9);
    //RooAddModel resoGm("resoGm","double gaussian",RooArgList(resoGm1,resoGm2), frac_core);
    RooAddModel resoGm("resoGm","double gaussian",RooArgList(resoGm1,resoGm2,resoGm3),RooArgList(reso1frac,reso2frac));

    // do the fit
    resoSigma3.setConstant(true);
    resoGm.fitTo(*binnedData);
    cout << "---------------------------------------------------------------------- this was with s3 const" << endl;
    resoSigma1.setConstant(true);
    resoSigma2.setConstant(true);
    resoSigma3.setConstant(false);
    resoGm.fitTo(*binnedData);
    cout << "---------------------------------------------------------------------- this was with s1,s2 const" << endl;
    resoSigma1.setConstant(false);
    resoSigma2.setConstant(false);
    resoSigma3.setConstant(true);
    resoGm.fitTo(*binnedData);
    cout << "---------------------------------------------------------------------- this was with s3 const" << endl;
    resoSigma1.setConstant(true);
    resoSigma2.setConstant(true);
    resoSigma3.setConstant(false);
    resoGm.fitTo(*binnedData);
    cout << "---------------------------------------------------------------------- this was with s1,s2 const" << endl;
    resoSigma1.setConstant(false);
    resoSigma2.setConstant(false);
    resoSigma3.setConstant(true);
    resoGm.fitTo(*binnedData);
    cout << "---------------------------------------------------------------------- this was with s3 const" << endl;

    // Plotten
    RooPlot* frame = tDiff.frame(Title((title+" - Resolution function").c_str()));
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

    TCanvas *canvas = new TCanvas("fitres","fit results",800,800);
    canvas->Divide(2,2);
    canvas->cd(1);
    frame->Draw();

    // Now we freeze all parameters of the resolution function
    resoMean.setConstant(true);
    resoSigma1.setConstant(false);
    resoSigma2.setConstant(false);
    resoSigma3.setConstant(true);
    reso1frac.setConstant(false);
    reso2frac.setConstant(false);
    //reso3frac.setConstant(true);

    // Fit parameters
    //const double tauLb(1.391e-12), tauB0(1.525e-12);
    //const double tauLb(1.229e-12), tauB0(1.536e-12); // corresponds to the evt.pdl table used in production
    const double tauLb(1.424e-12), tauB0(1.536e-12); // corresponds to the evt.pdl table used in production
    const double tauTruthVal = (isB0 ? tauB0 : tauLb);
    RooRealVar tau("tau","tau", tauTruthVal, 0e-12, (tHi < 2*tauTruthVal ? 2*tauTruthVal : tHi));

    // lifetime fit setup
    RooGaussModel resmodelGm1("resmodelGm1","gauss model center",t,resoMean,resoSigma1) ;
    RooGaussModel resmodelGm2("resmodelGm2","gauss model tail",t,resoMean,resoSigma2) ;
    RooGaussModel resmodelGm3("resmodelGm3","gauss model extended tail",t,resoMean,resoSigma3) ;
    RooAddModel resmodel("resmodel","multigaussian",RooArgList(resmodelGm1,resmodelGm2,resmodelGm3),RooArgList(reso1frac,reso2frac));
    RooDecay decay("decay", "decay", t, tau, resmodel, RooDecay::SingleSided);

    // final model
    RooAddPdf model("model", "model", RooArgSet(decay), RooArgList(nsig));

    // do the fit
    model.fitTo(data);

    // some plotting
    canvas->cd(2);
    RooPlot * frame_t = t.frame(Title((title+" - Lifetime").c_str()));
    data.plotOn(frame_t, Binning(nBins));
    //modelTruth.plotOn(frame_t,ProjWData(data),LineStyle(kDashed),LineColor(kGreen));
    model.plotOn(frame_t,ProjWData(data));
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

    // draw fit of resofunction again
    // Plotten
    RooPlot* frame2 = tDiff.frame(Title((title+" - Resolution function").c_str()));
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

    canvas->cd(3);
    frame2->Draw();
    // we reached the end without an error
    return 0;
}

int doResoFitGenReco(string filename, string title, bool isB0, double tLo = -5e-12, double tHi = 15e-12)
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
    return doResoFitGenReco(tree, title, isB0, tLo, tHi);
}

