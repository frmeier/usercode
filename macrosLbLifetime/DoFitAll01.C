#include <iostream>
#include <string>
#include <memory>

#ifndef __CINT__
#include "boost/shared_ptr.hpp"
#endif /* __CINT __ */

#include "TROOT.h"
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
#include "RooVoigtian.h"
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

#include "DoFitAll01.h"

using std::cout;
using std::endl;
using std::string;
using std::auto_ptr;

/* DoFitAll01
 * -------------------------------------------------
 * Do something like this:
 .L DoFitAll01.C+
 dfa = new DoFitAll01;
 dfa->doFitDataB0mass();

 Dataset files are in init*datasets()
 */

DoFitAll01::DoFitAll01()
{
    Init();
}

DoFitAll01::~DoFitAll01()
{
}

void DoFitAll01::Init()
{
    fnCPU_ = 6; // number of CPU used for fitting and plotting
    fVerbose_ = 3;
    fnBins_t_ = 100; // choose an even number of bins here
    fnBins_m_ = 200; // choose an even number of bins here

    //rrvMCB0_m_ = new RooRealVar("rrvMCB0_m", "mass", 4.6, 6.0);
    //rrvMCB0_t_ = new RooRealVar("rrvMCB0_t", "decay time", -5e-12, 20e-12);
    //rrvMCB0_dt_ = new RooRealVar("rrvMCB0_dt", "time difference truth reco", -8e-12, 8e-12);

    //rrvMCB0_m_ = new RooRealVar("mass", "mass", 4.6, 6.0, "GeV/c^{2}");
    rrvMCB0_m_ = new RooRealVar("mass", "mass", 5.14, 5.50, "GeV/c^{2}"); // TODO/WORKAROUND: Bgr wird nicht richtig berechnet, siehe Journal 13, p.90
    rrvMCB0_t_ = new RooRealVar("t", "decay time", -5e-12, 20e-12, "s");
    rrvMCB0_dt_ = new RooRealVar("tDiff", "time difference truth reco", -8e-12, 8e-12, "s");

    rrvDatLb_m_ = new RooRealVar("mass", "mass", 5.35, 6.10, "GeV/c^{2}");
    rrvDatLb_t_ = new RooRealVar("t", "decay time", -5e-12, 20e-12, "s");

    setTDRStyle();

    isPrelim = false;
    isWorkInProgress = true;
    noTitle = false;
}

void DoFitAll01::initMCB0datasets()
{
    /*
    if (!ds_MCB0_.isInitialised)
    {
	ds_MCB0_.filename = "../data/vrt_r276_nocut.root";
	ds_MCB0_.treename = "fittree";
	openFileAndGetTree(ds_MCB0_);
	//openFileAndGetTree(ds_MCB0_, "vrtMCB0");
	ds_MCB0_.rds = new RooDataSet("rds_dataMCB0", "dataset MC B0 all", RooArgSet(*rrvMCB0_m_, *rrvMCB0_t_, *rrvMCB0_dt_),
		RooFit::Import(*ds_MCB0_.tree));
	ds_MCB0_.isInitialised = true;
    }*/

    if (!ds_MCB0match_.isInitialised)
    {
	ds_MCB0match_.filename = "../data/vrt_r316_B0_mc_match_barrel.root";
	ds_MCB0match_.treename = "fittree";
	openFileAndGetTree(ds_MCB0match_);
	//openFileAndGetTree(ds_MCB0match_, "vrtMCB0match");
	ds_MCB0match_.rds = new RooDataSet("rds_dataMCB0match", "dataset MC B0 truth matched events",
		RooArgSet(*rrvMCB0_m_, *rrvMCB0_t_, *rrvMCB0_dt_), RooFit::Import(*ds_MCB0match_.tree));
	ds_MCB0match_.isInitialised = true;
    }

    /*
    if (!ds_MCB0nomatch_.isInitialised)
    {
	ds_MCB0nomatch_.filename = "../data/vrt_r276_noMCmatch.root";
	ds_MCB0nomatch_.treename = "fittree";
	openFileAndGetTree(ds_MCB0nomatch_);
	//openFileAndGetTree(ds_MCB0nomatch_, "vrtMCB0nomatch");
	ds_MCB0nomatch_.rds = new RooDataSet("rds_dataMCB0nomatch", "dataset MC B0 not truth matched events",
		RooArgSet(*rrvMCB0_m_, *rrvMCB0_t_, *rrvMCB0_dt_), RooFit::Import(*ds_MCB0nomatch_.tree));
	ds_MCB0nomatch_.isInitialised = true;
    }*/
}

void DoFitAll01::initDataB0datasets()
{
    if (!ds_dataB0barrel.isInitialised)
    {
	ds_dataB0barrel.filename = "../data/vrt_r315_B0_data_barrel.root";
	ds_dataB0barrel.treename = "fittree";
	openFileAndGetTree(ds_dataB0barrel);
	ds_dataB0barrel.rds = new RooDataSet("rds_dataDataB0barrel", "dataset data B0 barrel", RooArgSet(*rrvMCB0_m_, *rrvMCB0_t_),
		RooFit::Import(*ds_dataB0barrel.tree));
	ds_dataB0barrel.isInitialised = true;
    }
}

void DoFitAll01::initDataLbdatasets()
{
    if (!ds_dataLbbarrel.isInitialised)
    {
	ds_dataLbbarrel.filename = "../data/vrt_r313_lb_data_barrel.root";
	ds_dataLbbarrel.treename = "fittree";
	openFileAndGetTree(ds_dataLbbarrel);
	ds_dataLbbarrel.rds = new RooDataSet("rds_dataDataLbbarrel", "dataset data Lb barrel", RooArgSet(*rrvDatLb_m_, *rrvDatLb_t_),
		RooFit::Import(*ds_dataLbbarrel.tree));
	ds_dataLbbarrel.isInitialised = true;
    }
}

void DoFitAll01::initMCLbdatasets()
{
    if (!ds_MCLbbarrel.isInitialised)
    {
	ds_MCLbbarrel.filename = "../data/vrt_r313_lb_data_barrel.root";
	ds_MCLbbarrel.treename = "fittree";
	openFileAndGetTree(ds_MCLbbarrel);
	ds_MCLbbarrel.rds = new RooDataSet("rds_dataDataLbbarrel", "dataset data Lb barrel", RooArgSet(*rrvDatLb_m_, *rrvDatLb_t_),
		RooFit::Import(*ds_MCLbbarrel.tree));
	ds_MCLbbarrel.isInitialised = true;
    }
}

// fit the resolution function. Assumes clean, signal only events, i.e. truth matched MC
void DoFitAll01::doFitMCB0reso()
{
    initMCB0datasets();
    mtgMCB0reso.name   = "MCB0reso";
    mtgMCB0reso.title  = "MC B0 resolution";
    mtgMCB0reso.initMean(0.0, -1e-12, 1e-12);
    mtgMCB0reso.initSigma1(0.15e-12, 0.001e-12, 0.5e-12);
    mtgMCB0reso.initSigma2(0.30e-12, 0.1e-12,   1.7e-12);
    mtgMCB0reso.initSigma3(1.80e-12, 0.5e-12,   50e-12);
    mtgMCB0reso.initFrac1(.4, .01, .99);
    mtgMCB0reso.initFrac2(.4, .01, .99);
    mtgMCB0reso.setVar(*rrvMCB0_dt_);
    mtgMCB0reso.initModel();

    mtgMCB0reso.setConstant(false, false, false, true, false, false);
    mtgMCB0reso.model->fitTo(*ds_MCB0match_.rds, RooFit::NumCPU(fnCPU_));
    /*
    mtgMCB0reso.setConstant(false, true, true, false, false, false);
    mtgMCB0reso.model->fitTo(*ds_MCB0match_.rds, RooFit::NumCPU(fnCPU_));
    mtgMCB0reso.setConstant(false, false, false, true, false, false);
    mtgMCB0reso.model->fitTo(*ds_MCB0match_.rds, RooFit::NumCPU(fnCPU_));
    mtgMCB0reso.setConstant(false, true, true, false, false, false);
    mtgMCB0reso.model->fitTo(*ds_MCB0match_.rds, RooFit::NumCPU(fnCPU_));
    mtgMCB0reso.setConstant(false, false, false, true, false, false);
    mtgMCB0reso.model->fitTo(*ds_MCB0match_.rds, RooFit::NumCPU(fnCPU_));
    */

    curFrame = auto_ptr<RooPlot>(rrvMCB0_dt_->frame(RooFit::Title((mtgMCB0reso.title+" - "+ds_MCB0match_.rds->GetTitle()).c_str())));
    ds_MCB0match_.rds->plotOn(curFrame.get(), RooFit::Binning(fnBins_t_-1)); // make an odd number of bins
    mtgMCB0reso.plotOnFrame(curFrame.get(), fnCPU_);
    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresReso","fit results",1200,800));
    mtgMCB0reso.writeResultsOnFrame(curFrame.get(), 1e12, "ps");
    curFrame->Draw();
    curCanvas->SaveAs(("DoFitAll01plot_"+mdtgMCB0match.name+".pdf").c_str());
}

// fit the resolution function. Assumes clean, signal only events, i.e. truth matched MC
void DoFitAll01::doFitMCB0matchTau(bool floatingTau, bool floatingReso)
{
    initMCB0datasets();

    mdtgMCB0match.name   = "MCB0matchTau";
    mdtgMCB0match.title  = "MC B0 tau fit on matched MC";
    mdtgMCB0match.initTau(1.535e-12, 0e-12, 3e-12);
    mdtgMCB0match.initMean(0.0, -1e-12, 1e-12);
    mdtgMCB0match.initSigma1(0.15e-12, 0.001e-12, 0.5e-12);
    mdtgMCB0match.initSigma2(0.30e-12, 0.1e-12,   1.7e-12);
    mdtgMCB0match.initSigma3(1.80e-12, 0.5e-12,   50e-12);
    mdtgMCB0match.initFrac1(.4, .01, .99);
    mdtgMCB0match.initFrac2(.4, .01, .99);
    mdtgMCB0match.setVar(*rrvMCB0_t_);
    mdtgMCB0match.initModel();

    doFitMCB0matchTauWorker(floatingTau, floatingReso);
}

void DoFitAll01::doFitMCB0matchTau(const modelTripleGauss &mtg, bool floatingTau, bool floatingReso)
{
    initMCB0datasets();

    mdtgMCB0match.name   = "MCB0matchTau";
    mdtgMCB0match.title  = "MC B0 tau fit on matched MC";
    mdtgMCB0match.initTau(1.535e-12, 0e-12, 3e-12);
    mdtgMCB0match.copyFromModel(mtg);
    mdtgMCB0match.setVar(*rrvMCB0_t_);
    mdtgMCB0match.initModel();

    doFitMCB0matchTauWorker(floatingTau, floatingReso);
}

void DoFitAll01::doFitMCB0matchTauWorker(const bool &floatingTau, const bool &floatingReso)
{
    //int i(0);
    //cout << "ok " << ++i << endl;
    if (floatingTau && floatingReso) mdtgMCB0match.setConstant(false, false, false, false, true, false, false);
    if (!floatingTau && floatingReso) mdtgMCB0match.setConstant(true, false, false, false, true, false, false);
    if (floatingTau && !floatingReso) mdtgMCB0match.setConstant(false, true, true, true, true, true, true);
    if (!(floatingTau || floatingReso))
    {
	mdtgMCB0match.setConstant(true, true, true, true, true, true, true);
	cout << "WARNING (DoFitAll01::doFitMCB0matchTau): "
	     << "neither tau nor resolution function requested as floating - performing a fit makes no sense. No guarantee for results..." << endl;
    }

    mdtgMCB0match.model->fitTo(*ds_MCB0match_.rds, RooFit::NumCPU(fnCPU_));

    curFrame = auto_ptr<RooPlot>(rrvMCB0_t_->frame(RooFit::Title((mdtgMCB0match.title+" - "+ds_MCB0match_.rds->GetTitle()).c_str())));
    ds_MCB0match_.rds->plotOn(curFrame.get(), RooFit::Binning(fnBins_t_));
    mdtgMCB0match.plotOnFrame(curFrame.get(), fnCPU_);
    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresTau","fit results",1200,800));
    mdtgMCB0match.writeResultsOnFrame(curFrame.get(), 1e12, "ps");
    curFrame->Draw();
    curCanvas->SaveAs(("DoFitAll01plot_"+mdtgMCB0match.name+".pdf").c_str());
}

// fit the mass. Assumes clean, signal only events, i.e. truth matched MC
void DoFitAll01::doFitMCB0mass()
{
    initMCB0datasets();
    mdgMCB0mass.name   = "MCB0mass";
    mdgMCB0mass.title  = "MC B0 mass";
    mdgMCB0mass.initMean(5.28, 5.14, 5.50);
    mdgMCB0mass.initSigma1(0.008, 0.002, 0.020);
    mdgMCB0mass.initSigma2(0.018, 0.012, 0.080);
    mdgMCB0mass.initFrac1(.4, .01, .99);
    mdgMCB0mass.setVar(*rrvMCB0_m_);
    mdgMCB0mass.initModel();
    rrvMCB0_m_->setRange("B0massrange", 5.14, 5.50);

    mdgMCB0mass.setConstant(false, false, false, false);
    mdgMCB0mass.model->fitTo(*ds_MCB0match_.rds, RooFit::NumCPU(fnCPU_), RooFit::Range("B0massrange"));
    mdgMCB0mass.calculateRange(2.0);

    curFrame = auto_ptr<RooPlot>(rrvMCB0_m_->frame(RooFit::Title((mdgMCB0mass.title+" - "+ds_MCB0match_.rds->GetTitle()).c_str()),RooFit::Range("B0massrange")));
    ds_MCB0match_.rds->plotOn(curFrame.get(), RooFit::Binning(fnBins_m_-1), RooFit::Range("B0massrange")); // make an odd number of bins
    mdgMCB0mass.plotOnFrame(curFrame.get(), fnCPU_);
    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresMass","fit results",1200,800));
    mdgMCB0mass.writeResultsOnFrame(curFrame.get(), 1, "GeV/c^{2}", 1000, "MeV/c^{2}");
    curFrame->Draw();
    curCanvas->SaveAs(("DoFitAll01plot_"+mdgMCB0mass.name+".pdf").c_str());
}

// fit the mass. Currently assumes clean, signal only events, i.e. truth matched MC. But intended as basis to extend for data
void DoFitAll01::doFitDataB0mass()
{
    //initMCB0datasets();
    initDataB0datasets();
    //DoFitAll01_dataset *dset = &ds_MCB0match_;
    DoFitAll01_dataset *dset = &ds_dataB0barrel;
    msgDataB0mass.name   = "DataB0mass";
    msgDataB0mass.title  = "Data B^{0} mass";
    msgDataB0mass.initMean(5.28, 5.14, 5.50);
    msgDataB0mass.initSigma1(0.010, 0.005, 0.020);
    msgDataB0mass.initSigma2(0.020, 0.012, 0.080);
    msgDataB0mass.initFrac1(0.4, 0.0, 1.0);
    msgDataB0mass.initc0(0.0, -1.0, 1.0);
    const int nentries = dset->rds->numEntries();
    msgDataB0mass.initNsig(.5*nentries, 0, nentries);
    msgDataB0mass.initNbgr(.5*nentries, 0, nentries);
    msgDataB0mass.setVar(*rrvMCB0_m_);
    msgDataB0mass.initModel();
    rrvMCB0_m_->setRange("B0massrange", 5.14, 5.50);
    //msgDataB0mass.model->fixCoefRange("B0massrange");

    msgDataB0mass.setConstant(false, false, false, false, false, false, false);
    //msgDataB0mass.setConstant(false, false, false);
    msgDataB0mass.model->fitTo( *(dset->rds), RooFit::NumCPU(fnCPU_), RooFit::Range("B0massrange"));
    msgDataB0mass.calculateRange(2.0);
    msgDataB0mass.calculateSigInRange(2.0);
    cout << "2 sigma equivalent range calculated: " << msgDataB0mass.sigLo << " to " << msgDataB0mass.sigHi << endl;
    const Double_t signalLowerBound = msgDataB0mass.sigLo;
    const Double_t signalUpperBound = msgDataB0mass.sigHi;
    cout << "sig in " << msgDataB0mass.Nsig_sigmas << " sigmas: " << msgDataB0mass.Nsig_val << " +/- " << msgDataB0mass.Nsig_err << endl;
    cout << "bgr in " << msgDataB0mass.Nsig_sigmas << " sigmas: " << msgDataB0mass.Nbgr_val << " +/- " << msgDataB0mass.Nbgr_err << endl;
    msgDataB0mass.calculateRange(3.0);
    cout << "3 sigma equivalent range calculated for bgr selection: " << msgDataB0mass.sigLo << " to " << msgDataB0mass.sigHi << endl;
    const Double_t upperSidebandLowerBound = msgDataB0mass.sigHi;
    rrvMCB0_m_->setRange("B0_bgr_range", msgDataB0mass.sigHi, 5.50); // TODO: obere Limite ist hart codiert - sinnvoll anpassen

    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresMass","fit results",1200,800));
    curFrame = auto_ptr<RooPlot>(rrvMCB0_m_->frame(RooFit::Title((msgDataB0mass.title+" - "+dset->rds->GetTitle()).c_str()),RooFit::Range("B0massrange")));
    dset->rds->plotOn(curFrame.get(), RooFit::Binning(fnBins_m_-1), RooFit::Range("B0massrange")); // make an odd number of bins
    msgDataB0mass.plotOnFrame(curFrame.get(), fnCPU_);
    msgDataB0mass.writeResultsOnFrame(curFrame.get(), 1, "GeV/c^{2}", 1000, "MeV/c^{2}");
    curFrame->Draw();
    gPad->SetLogy(0);
    curCanvas->SaveAs(("DoFitAll01plot_"+msgDataB0mass.name+".pdf").c_str());

    //double upperSidebandLowerBound = 5.3295; // DEBUG

    mdgDataB0BgrT.name= "dataB0_t_bgr";
    mdgDataB0BgrT.title = "Data B^{0} lifetime background";
    mdgDataB0BgrT.initMean(0, -0.5e-12, 0.5e-12);
    mdgDataB0BgrT.initSigma1(.2e-12, .01e-12, 1e-12);
    mdgDataB0BgrT.initSigma2(1e-12, .5e-12, 3e-12);
    mdgDataB0BgrT.initFrac1(.5, 0,1);
    mdgDataB0BgrT.setVar(*rrvMCB0_t_);
    mdgDataB0BgrT.initModel();

    mdgDataB0BgrT.setConstant(false, false, false, false);
    //mdgDataB0BgrT.model->fitTo(*ds_MCB0match_.rds, RooFit::NumCPU(fnCPU_), RooFit::Range("B0_bgr_range"));
    string cutBgr = "mass>" + toString(upperSidebandLowerBound) + "&&mass<5.5"; // TODO: 8tung, obere Grenze hart codiert
    auto_ptr<RooDataSet> rds_bgr ((RooDataSet*) dset->rds->reduce(RooFit::Cut(cutBgr.c_str())));
    mdgDataB0BgrT.model->fitTo(*rds_bgr.get(), RooFit::NumCPU(fnCPU_), RooFit::Cut(cutBgr.c_str()));

    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresTime","fit results",1200,800));
    //curFrame = auto_ptr<RooPlot>(rrvMCB0_t_->frame(RooFit::Title((mdgDataB0BgrT.title+" - "+dset->rds->GetTitle()).c_str()),RooFit::Range("B0_bgr_range")));
    curFrame = auto_ptr<RooPlot>(rrvMCB0_t_->frame(RooFit::Title((mdgDataB0BgrT.title+" - "+dset->rds->GetTitle()).c_str()),RooFit::Cut(cutBgr.c_str())));
    //dset->rds->plotOn(curFrame.get(), RooFit::Binning(fnBins_t_-1), RooFit::Range("B0_bgr_range")); // make an odd number of bins
    rds_bgr->plotOn(curFrame.get(), RooFit::Binning(fnBins_t_-1), RooFit::Cut(cutBgr.c_str())); // make an odd number of bins
    mdgDataB0BgrT.plotOnFrame(curFrame.get(), fnCPU_);
    mdgDataB0BgrT.writeResultsOnFrame(curFrame.get(), 1e12, "ps", 1e12, "ps");
    curFrame->Draw();
    gPad->SetLogy(1);
    curCanvas->SaveAs(("DoFitAll01plot_"+mdgDataB0BgrT.name+".pdf").c_str());

    // -------------------------------------------------------------------
    // ------ and now we fit the lifetime with bgr
    mdtgbdgDataB0lifetime.name = "DataB0t_lifetime";
    mdtgbdgDataB0lifetime.title = "Data B0 lifetime fit with bgr";

    mdtgbdgDataB0lifetime.initTau(1.5e-12, 1e-12, 2e-12);

    // 1  MCB0reso_frac1   4.00929e-01   1.39919e-02   1.57613e-04  -2.03590e-01
    // 2  MCB0reso_frac2   2.82830e-01   1.30211e-02   2.19086e-04  -4.59169e-01
    // 3  MCB0reso_mean   3.12890e-15   2.54505e-15   1.87699e-04   3.12891e-03
    // 4  MCB0reso_sigma1   1.29598e-13   3.83095e-15   1.63254e-04  -5.05880e-01
    // 5  MCB0reso_sigma2   4.67371e-13   2.14021e-14   2.79640e-04  -5.71372e-01

    mdtgbdgDataB0lifetime.initResoMean(0e-12, -0.1e-12, 0.1e-12);
    mdtgbdgDataB0lifetime.initResoSigma1(1.29598e-13, 0.001e-12, 0.5e-12);
    mdtgbdgDataB0lifetime.initResoSigma2(4.67371e-13, 0.1e-12, 1.7e-12);
    mdtgbdgDataB0lifetime.initResoSigma3(1.8e-12, 0.5e-12, 50e-12);
    mdtgbdgDataB0lifetime.initResoFrac1(4.00929e-01, .01, .99);
    mdtgbdgDataB0lifetime.initResoFrac2(2.82830e-01, .01, .99);

    mdtgbdgDataB0lifetime.initBgrMean(mdgDataB0BgrT.mean->getVal(), 1e-12, 2e-12);
    mdtgbdgDataB0lifetime.initBgrSigma1(mdgDataB0BgrT.sigma1->getVal(), 0.001e-12, 0.5e-12);
    mdtgbdgDataB0lifetime.initBgrSigma2(mdgDataB0BgrT.sigma2->getVal(), 0.1e-12, 1.7e-12);
    mdtgbdgDataB0lifetime.initBgrGaussFrac1(mdgDataB0BgrT.frac1->getVal(), 0, 1);
    
    const Double_t sigfrac = msgDataB0mass.Nsig_val / (msgDataB0mass.Nsig_val + msgDataB0mass.Nbgr_val);
    cout << "The signal fraction is: " << sigfrac << endl;
    mdtgbdgDataB0lifetime.initSigFrac(sigfrac, 0, 1);

    mdtgbdgDataB0lifetime.setVar(*rrvMCB0_t_);
    mdtgbdgDataB0lifetime.initModel();

    mdtgbdgDataB0lifetime.setConstant(false, false, false, false, false, false, false);
    mdtgbdgDataB0lifetime.setConstantBgr(false, false, false, false); // keep bgr shape fixed
    mdtgbdgDataB0lifetime.setConstantSigFrac(false); // keep sig/bgr ratio fixed

    // now get the data from signal region
    string cutSig = "mass>" + toString(signalLowerBound) + "&&mass<" + toString(signalUpperBound); // TODO: 8tung, obere Grenze hart codiert
    auto_ptr<RooDataSet> rds_sig ((RooDataSet*) dset->rds->reduce(RooFit::Cut(cutSig.c_str())));

    mdtgbdgDataB0lifetime.full_model->fitTo(*rds_sig.get(), RooFit::NumCPU(fnCPU_));
    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresTime","fit results",1200,800));
    curFrame = auto_ptr<RooPlot>(rrvMCB0_t_->frame(RooFit::Title((noTitle ? " " : mdtgbdgDataB0lifetime.title+" - "+rds_sig->GetTitle()).c_str()),RooFit::Cut(cutBgr.c_str())));
    rds_sig->plotOn(curFrame.get(), RooFit::Binning(fnBins_t_));
    mdtgbdgDataB0lifetime.plotOnFrame(curFrame.get(), fnCPU_);
    mdtgbdgDataB0lifetime.writeResultsOnFrame(curFrame.get(), 1e12, "ps");
    if (isWorkInProgress) mdtgbdgDataB0lifetime.markWorkInProgress(curFrame.get());
    curFrame->Draw();
    gPad->SetLogy(1);
    gPad->SetRightMargin(0.09);
    curCanvas->SaveAs(("DoFitAll01plot_"+mdtgbdgDataB0lifetime.name+".pdf").c_str());
}

// fit the mass. Currently assumes clean, signal only events, i.e. truth matched MC. But intended as basis to extend for data
void DoFitAll01::doFitDataLbmass()
{
    initDataLbdatasets();
    DoFitAll01_dataset *dset = &ds_dataLbbarrel;
    msgDataLbmass.name   = "DataLbmass";
    msgDataLbmass.title  = "Data Lb mass";
    const Double_t massUpperBound = 6.10;
    msgDataLbmass.initMean(5.62, 5.35, massUpperBound);
    msgDataLbmass.initSigma1(0.010, 0.005, 0.020);
    msgDataLbmass.initSigma2(0.020, 0.012, 0.080);
    msgDataLbmass.initFrac1(0.4, 0.0, 1.0);
    msgDataLbmass.initc0(0.0, -1.0, 1.0);
    const int nentries = dset->rds->numEntries();
    msgDataLbmass.initNsig(.5*nentries, 0, nentries);
    msgDataLbmass.initNbgr(.5*nentries, 0, nentries);
    msgDataLbmass.setVar(*rrvDatLb_m_);
    msgDataLbmass.initModel();
    rrvDatLb_m_->setRange("Lbmassrange", 5.35, massUpperBound);
    //msgDataLbmass.model->fixCoefRange("Lbmassrange");

    msgDataLbmass.setConstant(false, false, false, false, false, false, false);
    //msgDataLbmass.setConstant(false, false, false);
    msgDataLbmass.model->fitTo( *(dset->rds), RooFit::NumCPU(fnCPU_), RooFit::Range("Lbmassrange"));
    msgDataLbmass.calculateRange(2.0);
    msgDataLbmass.calculateSigInRange(2.0);
    cout << "2 sigma equivalent range calculated: " << msgDataLbmass.sigLo << " to " << msgDataLbmass.sigHi << endl;
    const Double_t signalLowerBound = msgDataLbmass.sigLo;
    const Double_t signalUpperBound = msgDataLbmass.sigHi;
    cout << "sig in " << msgDataLbmass.Nsig_sigmas << " sigmas: " << msgDataLbmass.Nsig_val << " +/- " << msgDataLbmass.Nsig_err << endl;
    cout << "bgr in " << msgDataLbmass.Nsig_sigmas << " sigmas: " << msgDataLbmass.Nbgr_val << " +/- " << msgDataLbmass.Nbgr_err << endl;
    msgDataLbmass.calculateRange(3.0);
    cout << "3 sigma equivalent range calculated for bgr selection: " << msgDataLbmass.sigLo << " to " << msgDataLbmass.sigHi << endl;
    const Double_t upperSidebandLowerBound = msgDataLbmass.sigHi;
    rrvDatLb_m_->setRange("Lb_bgr_range", msgDataLbmass.sigHi, massUpperBound); // TODO: obere Limite ist hart codiert - sinnvoll anpassen

    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresMass","fit results",1200,800));
    curFrame = auto_ptr<RooPlot>(rrvDatLb_m_->frame(RooFit::Title((msgDataLbmass.title+" - "+dset->rds->GetTitle()).c_str()),RooFit::Range("Lbmassrange")));
    dset->rds->plotOn(curFrame.get(), RooFit::Binning(fnBins_m_-1), RooFit::Range("Lbmassrange")); // make an odd number of bins
    msgDataLbmass.plotOnFrame(curFrame.get(), fnCPU_);
    msgDataLbmass.writeResultsOnFrame(curFrame.get(), 1, "GeV/c^{2}", 1000, "MeV/c^{2}");
    curFrame->Draw();
    curCanvas->SaveAs(("DoFitAll01plot_"+msgDataLbmass.name+".pdf").c_str());

    //double upperSidebandLowerBound = 5.3295; // DEBUG

    mdgDataLbBgrT.title = "Data Lb lifetime background";
    mdgDataLbBgrT.initMean(0, -0.5e-12, 0.5e-12);
    mdgDataLbBgrT.initSigma1(.2e-12, .01e-12, 1e-12);
    mdgDataLbBgrT.initSigma2(1e-12, .5e-12, 1.1e-12);
    mdgDataLbBgrT.initFrac1(.5, 0,1);
    mdgDataLbBgrT.setVar(*rrvDatLb_t_);
    mdgDataLbBgrT.initModel();

    mdgDataLbBgrT.setConstant(false, false, false, false);
    //mdgDataLbBgrT.model->fitTo(*ds_MCLbmatch_.rds, RooFit::NumCPU(fnCPU_), RooFit::Range("Lb_bgr_range"));
    string cutBgr = "mass>" + toString(upperSidebandLowerBound) + "&&mass<" + toString(massUpperBound); // TODO: 8tung, obere Grenze hart codiert
    auto_ptr<RooDataSet> rds_bgr ((RooDataSet*) dset->rds->reduce(RooFit::Cut(cutBgr.c_str())));
    mdgDataLbBgrT.model->fitTo(*rds_bgr.get(), RooFit::NumCPU(fnCPU_), RooFit::Cut(cutBgr.c_str()));

    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresTime","fit results",1200,800));
    //curFrame = auto_ptr<RooPlot>(rrvDatLb_t_->frame(RooFit::Title((mdgDataLbBgrT.title+" - "+dset->rds->GetTitle()).c_str()),RooFit::Range("Lb_bgr_range")));
    curFrame = auto_ptr<RooPlot>(rrvDatLb_t_->frame(RooFit::Title((mdgDataLbBgrT.title+" - "+dset->rds->GetTitle()).c_str()),RooFit::Cut(cutBgr.c_str())));
    //dset->rds->plotOn(curFrame.get(), RooFit::Binning(fnBins_t_-1), RooFit::Range("Lb_bgr_range")); // make an odd number of bins
    rds_bgr->plotOn(curFrame.get(), RooFit::Binning(fnBins_t_-1), RooFit::Cut(cutBgr.c_str())); // make an odd number of bins
    mdgDataLbBgrT.plotOnFrame(curFrame.get(), fnCPU_);
    mdgDataLbBgrT.writeResultsOnFrame(curFrame.get(), 1e12, "ps", 1e12, "ps");
    curFrame->Draw();
    curCanvas->SaveAs(("DoFitAll01plot_"+mdgDataLbBgrT.name+".pdf").c_str());

    // -------------------------------------------------------------------
    // ------ and now we fit the lifetime with bgr
    mdtgbdgDataLblifetime.name = "DataLbt_lifetime";
    mdtgbdgDataLblifetime.title = "Data Lb lifetime fit with bgr";

    mdtgbdgDataLblifetime.initTau(1.5e-12, 1e-12, 2e-12);

    // 1  MCLbreso_frac1   4.00929e-01   1.39919e-02   1.57613e-04  -2.03590e-01
    // 2  MCLbreso_frac2   2.82830e-01   1.30211e-02   2.19086e-04  -4.59169e-01
    // 3  MCLbreso_mean   3.12890e-15   2.54505e-15   1.87699e-04   3.12891e-03
    // 4  MCLbreso_sigma1   1.29598e-13   3.83095e-15   1.63254e-04  -5.05880e-01
    // 5  MCLbreso_sigma2   4.67371e-13   2.14021e-14   2.79640e-04  -5.71372e-01

    mdtgbdgDataLblifetime.initResoMean(0e-12, -0.1e-12, 0.1e-12);
    mdtgbdgDataLblifetime.initResoSigma1(1.29598e-13, 0.001e-12, 0.5e-12);
    mdtgbdgDataLblifetime.initResoSigma2(4.67371e-13, 0.1e-12, 1.7e-12);
    mdtgbdgDataLblifetime.initResoSigma3(1.8e-12, 0.5e-12, 50e-12);
    mdtgbdgDataLblifetime.initResoFrac1(4.00929e-01, .01, .99);
    mdtgbdgDataLblifetime.initResoFrac2(2.82830e-01, .01, .99);

    mdtgbdgDataLblifetime.initBgrMean(mdgDataLbBgrT.mean->getVal(), 1e-12, 2e-12);
    mdtgbdgDataLblifetime.initBgrSigma1(mdgDataLbBgrT.sigma1->getVal(), 0.001e-12, 0.5e-12);
    mdtgbdgDataLblifetime.initBgrSigma2(mdgDataLbBgrT.sigma2->getVal(), 0.1e-12, 1.7e-12);
    mdtgbdgDataLblifetime.initBgrGaussFrac1(mdgDataLbBgrT.frac1->getVal(), 0, 1);
    
    const Double_t sigfrac = msgDataLbmass.Nsig_val / (msgDataLbmass.Nsig_val + msgDataLbmass.Nbgr_val);
    cout << "The signal fraction is: " << sigfrac << endl;
    mdtgbdgDataLblifetime.initSigFrac(sigfrac, 0, 1);

    mdtgbdgDataLblifetime.setVar(*rrvDatLb_t_);
    mdtgbdgDataLblifetime.initModel();

    mdtgbdgDataLblifetime.setConstant(false, false, false, false, false, false, false);
    mdtgbdgDataLblifetime.setConstantBgr(false, false, false, false); // keep bgr shape fixed
    mdtgbdgDataLblifetime.setConstantSigFrac(false); // keep sig/bgr ratio fixed

    // now get the data from signal region
    string cutSig = "mass>" + toString(signalLowerBound) + "&&mass<" + toString(signalUpperBound); // TODO: 8tung, obere Grenze hart codiert
    auto_ptr<RooDataSet> rds_sig ((RooDataSet*) dset->rds->reduce(RooFit::Cut(cutSig.c_str())));

    mdtgbdgDataLblifetime.full_model->fitTo(*rds_sig.get(), RooFit::NumCPU(fnCPU_));
    if (curCanvas.get() != 0) curCanvas.reset();
    curCanvas = auto_ptr<TCanvas>(new TCanvas("fitresTime","fit results",1200,800));
    curFrame = auto_ptr<RooPlot>(rrvDatLb_t_->frame(RooFit::Title((mdtgbdgDataLblifetime.title+" - "+rds_sig->GetTitle()).c_str()),RooFit::Cut(cutBgr.c_str())));
    rds_sig->plotOn(curFrame.get(), RooFit::Binning(fnBins_t_-1)); // make an odd number of bins
    mdtgbdgDataLblifetime.plotOnFrame(curFrame.get(), fnCPU_);
    mdtgbdgDataLblifetime.writeResultsOnFrame(curFrame.get(), 1e12, "ps");
    curFrame->Draw();
    gPad->SetLogy(1);
    curCanvas->SaveAs(("DoFitAll01plot_"+mdtgbdgDataLblifetime.name+".pdf").c_str());
}

// -------------------------------------------------------------------
TFile* DoFitAll01::openFile(std::string filename)
{
    if (fVerbose_ > 4) cout << "Attempting to open file " << filename << "..." << endl;
    TFile* f = (TFile::Open(filename.c_str()));
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
	throw;
    }
    if (fVerbose_ > 2) cout << "File " << filename << " opened successfully." << endl;
    return f;
}

TTree* DoFitAll01::getTreeFromFile(std::string treename, TFile* f)
{
    TTree* t = (TTree*)f->Get(treename.c_str());
    if (t==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
	cout << "This script requires a very reduced tree called fittree" << endl;
    }
    if (fVerbose_ > 2) cout << "TTree " << treename << " accessed successfully." << endl;
    return t;
}

TTree* DoFitAll01::getTreeFromFile(std::string treename, std::string rename, TFile* f)
{
    TTree* t = (TTree*)f->Get(treename.c_str())->Clone(rename.c_str());
    if (t==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
	cout << "This script requires a very reduced tree called fittree" << endl;
    }
    if (fVerbose_ > 2) cout << "TTree " << treename << " accessed successfully." << endl;
    return t;
}

TTree* DoFitAll01::openFileAndGetTree(std::string filename, std::string treename, TFile* &file)
{
    file = openFile(filename);
    return getTreeFromFile(treename, file);
}

void DoFitAll01::openFileAndGetTree(DoFitAll01_dataset &ds)
{
    ds.tree = openFileAndGetTree(ds.filename, ds.treename, ds.file);
}

void DoFitAll01::openFileAndGetTree(DoFitAll01_dataset &ds, std::string rename)
{
    ds.file = openFile(ds.filename);
    ds.tree = getTreeFromFile(ds.treename, rename, ds.file);
}

