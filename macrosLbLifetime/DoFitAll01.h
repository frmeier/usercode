#ifndef GUARD_doFitAll01_H
#define GUARD_doFitAll01_H

#include <memory>
#include <string>
#include <iostream>

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooPlot.h"

#include "utils.h"

// Class to fit all stuff
class DoFitAll01
{
    public:
	DoFitAll01();
	~DoFitAll01();
	void Init();

	struct DoFitAll01_dataset
	{
	    DoFitAll01_dataset() : isInitialised(false) {};
	    bool isInitialised;
	    std::string filename;
	    std::string treename;
	    TFile* file;
	    TTree* tree;
	    RooDataSet* rds;
	};

	// ======================================================================================================= some containers for models
	struct modelTemplate
	{
	    std::string name;
	    std::string title;
	    RooRealVar* var;
	    Double_t sigSigmas, sigLo, sigHi;
	    std::auto_ptr<RooAddModel> model;
	    void setVar(RooRealVar& v) { var = &v; } ;
	    void calculateRange(Double_t sigmas)
	    {
		sigSigmas = sigmas;
		const Double_t area = .5*(TMath::Erf(-sigSigmas/TMath::Sqrt2())+1);
		std::auto_ptr<RooAbsReal> cdf (model->createCdf(*var));
		sigLo = cdf->findRoot(*var, var->getMin(), var->getMax(), area);
		sigHi = cdf->findRoot(*var, var->getMin(), var->getMax(), 1-area);
		cout << "sigLo: " << sigLo << " sigHi: " << sigHi << " \"sigma\": " << (sigHi-sigLo)/(2*sigSigmas) << " center: " << .5*(sigHi+sigLo) << endl;
	    };
	};

	struct modelTripleGauss : public modelTemplate
	{
	    std::auto_ptr<RooRealVar> mean, sigma1, sigma2, sigma3, frac1, frac2;
	    std::auto_ptr<RooGaussModel> gm1, gm2, gm3;
	    void initMean(double val, double lo, double hi)
		{ mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_mean").c_str(), (title+" mean").c_str(), val, lo, hi)); };
	    void initSigma1(double val, double lo, double hi)
		{ sigma1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma1").c_str(), (title+" sigma1").c_str(), val, lo, hi)); };
	    void initSigma2(double val, double lo, double hi)
		{ sigma2 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma2").c_str(), (title+" sigma2").c_str(), val, lo, hi)); };
	    void initSigma3(double val, double lo, double hi)
		{ sigma3 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma3").c_str(), (title+" sigma3").c_str(), val, lo, hi)); };
	    void initFrac1(double val, double lo, double hi)
		{ frac1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_frac1").c_str(), (title+" frac1").c_str(), val, lo, hi)); };
	    void initFrac2(double val, double lo, double hi)
		{ frac2 = auto_ptr<RooRealVar>(new RooRealVar((name+"_frac2").c_str(), (title+" frac2").c_str(), val, lo, hi)); };
	    void copyFromModel(const modelTripleGauss &mtg)
	    {
		mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_mean").c_str(), (title+" mean").c_str(),
			    mtg.mean->getVal(), mtg.mean->getMin(), mtg.mean->getMax()));
		sigma1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma1").c_str(), (title+" sigma1").c_str(),
			    mtg.sigma1->getVal(), mtg.sigma1->getMin(), mtg.sigma1->getMax()));
		sigma2 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma2").c_str(), (title+" sigma2").c_str(),
			    mtg.sigma2->getVal(), mtg.sigma2->getMin(), mtg.sigma2->getMax()));
		sigma3 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma3").c_str(), (title+" sigma3").c_str(),
			    mtg.sigma3->getVal(), mtg.sigma3->getMin(), mtg.sigma3->getMax()));
		frac1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_frac1").c_str(), (title+" frac1").c_str(),
			    mtg.frac1->getVal(), mtg.frac1->getMin(), mtg.frac1->getMax()));
		frac2 = auto_ptr<RooRealVar>(new RooRealVar((name+"_frac2").c_str(), (title+" frac2").c_str(),
			    mtg.frac2->getVal(), mtg.frac2->getMin(), mtg.frac2->getMax()));
		mean->setError(mtg.mean->getError());
		sigma1->setError(mtg.sigma1->getError());
		sigma2->setError(mtg.sigma2->getError());
		sigma3->setError(mtg.sigma3->getError());
		frac1->setError(mtg.frac1->getError());
		frac2->setError(mtg.frac2->getError());
	    };
	    void initModel()
	    {
		gm1 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm1").c_str(), (title+" Gauss 1").c_str(), *var, *mean.get(), *sigma1.get()));
		gm2 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm2").c_str(), (title+" Gauss 2").c_str(), *var, *mean.get(), *sigma2.get()));
		gm3 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm3").c_str(), (title+" Gauss 3").c_str(), *var, *mean.get(), *sigma3.get()));
		model = auto_ptr<RooAddModel>(new RooAddModel((name+"_model_trigauss").c_str(), (title+" model (triple Gauss)").c_str(),
						RooArgList(*gm1.get(), *gm2.get(), *gm3.get()), RooArgList(*frac1, *frac2)));
	    }
	    void setConstant(bool m, bool s1, bool s2, bool s3, bool f1, bool f2)
	    {
		mean->setConstant(m);
		sigma1->setConstant(s1);
		sigma2->setConstant(s2);
		sigma3->setConstant(s3);
		frac1->setConstant(f1);
		frac2->setConstant(f2);
	    };
	    void plotOnFrame(RooPlot* frame, int nCPU)
	    {	// make sure that the full model is plotted as last function so that chiSquare reports the correct value afterwards
		model->plotOn(frame, RooFit::Components((name+"_gm1").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
		model->plotOn(frame, RooFit::Components((name+"_gm2").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kCyan));
		model->plotOn(frame, RooFit::Components((name+"_gm3").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
		model->plotOn(frame, RooFit::NumCPU(nCPU));
	    };
	    void writeResultsOnFrame(RooPlot* frame, double scale = 1, std::string unit = "")
	    {
		frameWriter fw;
		fw.setFrame(frame);
		fw.addLine("Triple Gaussian fit:");
		fw.addLine("mean: " + roundToString(mean->getVal()*scale,3) + " #pm " + roundToString(mean->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 1: " + roundToString(sigma1->getVal()*scale,3) + " #pm " + roundToString(sigma1->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 2: " + roundToString(sigma2->getVal()*scale,3) + " #pm " + roundToString(sigma2->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 3: " + roundToString(sigma3->getVal()*scale,3) + " #pm " + roundToString(sigma3->getError()*scale,3) + " " + unit);
		fw.addLine("frac 1: " + roundToString(frac1->getVal(),3) + " #pm " + roundToString(frac1->getError(),3));
		fw.addLine("frac 2: " + roundToString(frac2->getVal(),3) + " #pm " + roundToString(frac2->getError(),3));
		const double restVal = 1 - (frac1->getVal()+frac2->getVal());
		const double restErr = sqrt(frac1->getError()*frac1->getError()+frac2->getError()*frac2->getError());
		fw.addLine("frac 3 (rest): " + roundToString(restVal,3) + " #pm " + roundToString(restErr,3));
		fw.addLine("#chi^{2}/ndof: " + roundToString(frame->chiSquare(),3));
	    };
	};

	struct modelDoubleGauss : public modelTemplate
	{
	    std::auto_ptr<RooRealVar> mean, sigma1, sigma2, frac1;
	    std::auto_ptr<RooGaussModel> gm1, gm2;
	    void initMean(double val, double lo, double hi)
		{ mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_mean").c_str(), (title+" mean").c_str(), val, lo, hi)); };
	    void initMean(RooRealVar& rrv) { mean = auto_ptr<RooRealVar>(new RooRealVar(rrv, (name+"_mean").c_str())); };
	    void initSigma1(double val, double lo, double hi)
		{ sigma1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma1").c_str(), (title+" sigma1").c_str(), val, lo, hi)); };
	    void initSigma2(double val, double lo, double hi)
		{ sigma2 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma2").c_str(), (title+" sigma2").c_str(), val, lo, hi)); };
	    void initFrac1(double val, double lo, double hi)
		{ frac1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_frac1").c_str(), (title+" frac1").c_str(), val, lo, hi)); };
	    void initModel()
	    {
		gm1 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm1").c_str(), (title+" Gauss 1").c_str(), *var, *mean.get(), *sigma1.get()));
		gm2 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm2").c_str(), (title+" Gauss 2").c_str(), *var, *mean.get(), *sigma2.get()));
		model = auto_ptr<RooAddModel>(new RooAddModel((name+"_model_dblgauss").c_str(), (title+" model (double Gauss)").c_str(), RooArgList(*gm1.get(), *gm2.get()),
						RooArgList(*frac1)));
	    }
	    void setConstant(bool m, bool s1, bool s2, bool f1)
	    {
		mean->setConstant(m);
		sigma1->setConstant(s1);
		sigma2->setConstant(s2);
		frac1->setConstant(f1);
	    };
	    void plotOnFrame(RooPlot* frame, int nCPU)
	    {	// make sure that the full model is plotted as last function so that chiSquare reports the correct value afterwards
		model->plotOn(frame, RooFit::Components((name+"_gm1").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
		model->plotOn(frame, RooFit::Components((name+"_gm2").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kCyan));
		model->plotOn(frame, RooFit::NumCPU(nCPU));
	    };
	    void writeResultsOnFrame(RooPlot* frame, double scale = 1, std::string unit = "", double sigmaScale = -1, std::string sigmaUnit = "")
	    {
		if (sigmaScale < 0)
		{
		    sigmaScale = scale;
		    sigmaUnit = unit;
		}
		frameWriter fw;
		fw.setFrame(frame);
		fw.addLine("Double Gaussian fit:");
		fw.addLine("mean: " + roundToString(mean->getVal()*scale,3) + " #pm " + roundToString(mean->getError()*scale,3) + " " + unit);
		fw.textsize *= 0.8; fw.txtLineSpace *= 0.8;
		fw.addLine("");
		fw.addLine("sigma 1: " + roundToString(sigma1->getVal()*sigmaScale,3) + " #pm " + roundToString(sigma1->getError()*sigmaScale,3) + " " + sigmaUnit);
		fw.addLine("sigma 2: " + roundToString(sigma2->getVal()*sigmaScale,3) + " #pm " + roundToString(sigma2->getError()*sigmaScale,3) + " " + sigmaUnit);
		fw.addLine("frac 1: " + roundToString(frac1->getVal(),3) + " #pm " + roundToString(frac1->getError(),3));
		const double restVal = 1 - frac1->getVal();
		const double restErr = frac1->getError();
		fw.addLine("frac 2 (rest): " + roundToString(restVal,3) + " #pm " + roundToString(restErr,3));
		fw.addLine("#chi^{2}/ndof: " + roundToString(frame->chiSquare(),3));
	    };
	};

	struct modelSingleGauss : public modelTemplate
	{
	    std::auto_ptr<RooRealVar> mean, sigma1;
	    std::auto_ptr<RooGaussModel> gm1;
	    void initMean(double val, double lo, double hi)
		{ mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_mean").c_str(), (title+" mean").c_str(), val, lo, hi)); };
	    void initSigma1(double val, double lo, double hi)
		{ sigma1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma1").c_str(), (title+" sigma1").c_str(), val, lo, hi)); };
	    void copyFromModel(const modelSingleGauss &mtg)
	    {
		mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_mean").c_str(), (title+" mean").c_str(),
			    mtg.mean->getVal(), mtg.mean->getMin(), mtg.mean->getMax()));
		sigma1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma1").c_str(), (title+" sigma1").c_str(),
			    mtg.sigma1->getVal(), mtg.sigma1->getMin(), mtg.sigma1->getMax()));
		mean->setError(mtg.mean->getError());
		sigma1->setError(mtg.sigma1->getError());
	    };
	    void initModel()
	    {
		gm1 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm1").c_str(), (title+" Gauss 1").c_str(), *var, *mean.get(), *sigma1.get()));
		model = auto_ptr<RooAddModel>(new RooAddModel((name+"_model_sglgauss").c_str(), (title+" model (single Gauss)").c_str(),
						RooArgList(*gm1.get()), RooArgList()));
	    }
	    void setConstant(bool m, bool s1)
	    {
		mean->setConstant(m);
		sigma1->setConstant(s1);
	    };
	    void plotOnFrame(RooPlot* frame, int nCPU)
	    {	// make sure that the full model is plotted as last function so that chiSquare reports the correct value afterwards
		model->plotOn(frame, RooFit::NumCPU(nCPU));
	    };
	    void writeResultsOnFrame(RooPlot* frame, double scale = 1, std::string unit = "", double sigmaScale = -1, std::string sigmaUnit = "")
	    {
		if (sigmaScale < 0)
		{
		    sigmaScale = scale;
		    sigmaUnit = unit;
		}
		frameWriter fw;
		fw.setFrame(frame);
		fw.addLine("Single Gaussian fit:");
		fw.addLine("mean: " + roundToString(mean->getVal()*scale,3) + " #pm " + roundToString(mean->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 1: " + roundToString(sigma1->getVal()*sigmaScale,3) + " #pm " + roundToString(sigma1->getError()*sigmaScale,3) + " " + unit);
		fw.addLine("#chi^{2}/ndof: " + roundToString(frame->chiSquare(),3));
	    };
	    void calculateRange(Double_t sigmas)
	    {
		sigSigmas = sigmas;
		// calculating the range is trivial for a simple Gaussian....
		sigLo = mean->getVal() - sigSigmas * sigma1->getVal();
		sigHi = mean->getVal() + sigSigmas * sigma1->getVal();
	    };
	};

	struct modelVoigtian : public modelTemplate
	{
	    std::auto_ptr<RooRealVar> mean, sigma, width;
	    std::auto_ptr<RooVoigtian> vgt;
	    std::auto_ptr<RooAddPdf> model; // Voigtian not available as a model, using pdf unstead
	    void initMean(double val, double lo, double hi)
		{ mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_mean").c_str(), (title+" mean").c_str(), val, lo, hi)); };
	    void initSigma(double val, double lo, double hi)
		{ sigma = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma").c_str(), (title+" sigma").c_str(), val, lo, hi)); };
	    void initWidth(double val, double lo, double hi)
		{ width = auto_ptr<RooRealVar>(new RooRealVar((name+"_width").c_str(), (title+" width").c_str(), val, lo, hi)); };
	    void copyFromModel(const modelVoigtian &mvgt)
	    {
		mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_mean").c_str(), (title+" mean").c_str(),
			    mvgt.mean->getVal(), mvgt.mean->getMin(), mvgt.mean->getMax()));
		sigma = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigma").c_str(), (title+" sigma").c_str(),
			    mvgt.sigma->getVal(), mvgt.sigma->getMin(), mvgt.sigma->getMax()));
		width = auto_ptr<RooRealVar>(new RooRealVar((name+"_width").c_str(), (title+" width").c_str(),
			    mvgt.width->getVal(), mvgt.width->getMin(), mvgt.width->getMax()));
		mean->setError(mvgt.mean->getError());
		sigma->setError(mvgt.sigma->getError());
		width->setError(mvgt.width->getError());
	    };
	    void initModel()
	    {
		// last arg of RooVoigtian: fast algorithm selector
		vgt = auto_ptr<RooVoigtian>(new RooVoigtian((name+"_vgt").c_str(), (title+" Voigtian").c_str(), *var, *mean.get(), *width.get(), *sigma.get(), false));
		model = auto_ptr<RooAddPdf>(new RooAddPdf((name+"_model_vgt").c_str(), (title+" model (Voigtian)").c_str(), RooArgList(*vgt.get()), RooArgList()));
	    }
	    void setConstant(bool m, bool s, bool w)
	    {
		mean->setConstant(m);
		sigma->setConstant(s);
		width->setConstant(w);
	    };
	    void plotOnFrame(RooPlot* frame, int nCPU)
	    {	// make sure that the full model is plotted as last function so that chiSquare reports the correct value afterwards
		model->plotOn(frame, RooFit::NumCPU(nCPU));
	    };
	    void writeResultsOnFrame(RooPlot* frame, double scale = 1, std::string unit = "", double sigmaScale = -1, std::string sigmaUnit = "")
	    {
		if (sigmaScale < 0)
		{
		    sigmaScale = scale;
		    sigmaUnit = unit;
		}
		frameWriter fw;
		fw.setFrame(frame);
		fw.addLine("Voigtian fit:");
		fw.addLine("mean: " + roundToString(mean->getVal()*scale,3) + " #pm " + roundToString(mean->getError()*scale,3) + " " + unit);
		fw.addLine("sigma: " + roundToString(sigma->getVal()*sigmaScale,3) + " #pm " + roundToString(sigma->getError()*sigmaScale,3) + " " + unit);
		fw.addLine("width: " + roundToString(width->getVal()*sigmaScale,3) + " #pm " + roundToString(width->getError()*sigmaScale,3) + " " + unit);
		fw.addLine("#chi^{2}/ndof: " + roundToString(frame->chiSquare(),3));
	    };
	};


	struct modelDecayTripleGauss : public modelTripleGauss
	{
	    // TODO: vielleicht braucht es noch einen Normierungsfaktor.
	    std::auto_ptr<RooRealVar> tau;
	    std::auto_ptr<RooDecay> model;
	    std::auto_ptr<RooAddModel> resoModel;
	    void initTau(double val, double lo, double hi)
		{ tau = auto_ptr<RooRealVar>(new RooRealVar((name+"_tau").c_str(), (title+" tau").c_str(), val, lo, hi)); };
	    void initModel()
	    {
		gm1 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm1").c_str(), (title+" Gauss 1").c_str(), *var, *mean.get(), *sigma1.get()));
		gm2 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm2").c_str(), (title+" Gauss 2").c_str(), *var, *mean.get(), *sigma2.get()));
		gm3 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm3").c_str(), (title+" Gauss 3").c_str(), *var, *mean.get(), *sigma3.get()));
		resoModel = auto_ptr<RooAddModel>(new RooAddModel((name+"_resoModel").c_str(), (title+" resolution model").c_str(),
						    RooArgList(*gm1.get(), *gm2.get(), *gm3.get()), RooArgList(*frac1, *frac2)));
		model = auto_ptr<RooDecay>(new RooDecay((name+"_model_dec_trigauss").c_str(), (title+" model (decay plus triple Gauss)").c_str(), *var, *tau.get(),
						    *resoModel.get(), RooDecay::SingleSided));
	    }
	    void setConstant(bool t, bool m, bool s1, bool s2, bool s3, bool f1, bool f2)
	    {
		std::cout << "t m s1 s2 s3 f1 f2: " << t << m << s1 << s2 << s3 << f1 << f2 << std::endl;
		tau->setConstant(t);
		mean->setConstant(m);
		sigma1->setConstant(s1);
		sigma2->setConstant(s2);
		sigma3->setConstant(s3);
		frac1->setConstant(f1);
		frac2->setConstant(f2);
	    };
	    void plotOnFrame(RooPlot* frame, int nCPU)
	    {	// make sure that the full model is plotted as last function so that chiSquare reports the correct value afterwards
		model->plotOn(frame, RooFit::NumCPU(nCPU));
		//model->plotOn(frame, RooFit::Components((name+"_gm1").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
		//model->plotOn(frame, RooFit::Components((name+"_gm2").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kCyan));
		//model->plotOn(frame, RooFit::Components((name+"_gm3").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kMagenta));
	    };
	    void writeResultsOnFrame(RooPlot* frame, double scale = 1, std::string unit = "")
	    {
		frameWriter fw;
		fw.setFrame(frame);
		fw.addLine("#tau: " + roundToString(tau->getVal()*scale,3) + " #pm " + roundToString(tau->getError()*scale,3) + " " + unit);
		fw.textsize *= 0.8; fw.txtLineSpace *= 0.8;
		fw.addLine("");
		fw.addLine("Resolution function:");
		fw.addLine("mean: " + roundToString(mean->getVal()*scale,3) + " #pm " + roundToString(mean->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 1: " + roundToString(sigma1->getVal()*scale,3) + " #pm " + roundToString(sigma1->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 2: " + roundToString(sigma2->getVal()*scale,3) + " #pm " + roundToString(sigma2->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 3: " + roundToString(sigma3->getVal()*scale,3) + " #pm " + roundToString(sigma3->getError()*scale,3) + " " + unit);
		fw.addLine("frac 1: " + roundToString(frac1->getVal(),3) + " #pm " + roundToString(frac1->getError(),3));
		fw.addLine("frac 2: " + roundToString(frac2->getVal(),3) + " #pm " + roundToString(frac2->getError(),3));
		const double restVal = 1 - (frac1->getVal()+frac2->getVal());
		const double restErr = sqrt(frac1->getError()*frac1->getError()+frac2->getError()*frac2->getError());
		fw.addLine("frac 3 (rest): " + roundToString(restVal,3) + " #pm " + roundToString(restErr,3));
		fw.addLine("#chi^{2}/ndof: " + roundToString(frame->chiSquare(),3));
	    };
	};

	struct modelDecayTripleGaussBgrDoubleGauss : public modelTemplate
	{
	    std::auto_ptr<RooRealVar> tau;
	    std::auto_ptr<RooRealVar> reso_mean, reso_sigma1, reso_sigma2, reso_sigma3, reso_frac1, reso_frac2;
	    std::auto_ptr<RooRealVar> bgr_mean, bgr_sigma1, bgr_sigma2, bgr_frac1;
	    std::auto_ptr<RooRealVar> sig_frac;
	    std::auto_ptr<RooGaussModel> reso_gm1, reso_gm2, reso_gm3;
	    std::auto_ptr<RooGaussModel> bgr_gm1, bgr_gm2;
	    std::auto_ptr<RooDecay> sig_model;
	    std::auto_ptr<RooAddModel> reso_model, bgr_model;
	    std::auto_ptr<RooAddPdf> full_model;
	    void initTau(double val, double lo, double hi)
		{ tau = auto_ptr<RooRealVar>(new RooRealVar((name+"_tau").c_str(), (title+" tau").c_str(), val, lo, hi)); };
	    void initResoMean(double val, double lo, double hi)
		{ reso_mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_reso_mean").c_str(), (title+" reso_mean").c_str(), val, lo, hi)); };
	    void initResoSigma1(double val, double lo, double hi)
		{ reso_sigma1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_reso_sigma1").c_str(), (title+" reso_sigma1").c_str(), val, lo, hi)); };
	    void initResoSigma2(double val, double lo, double hi)
		{ reso_sigma2 = auto_ptr<RooRealVar>(new RooRealVar((name+"_reso_sigma2").c_str(), (title+" reso_sigma2").c_str(), val, lo, hi)); };
	    void initResoSigma3(double val, double lo, double hi)
		{ reso_sigma3 = auto_ptr<RooRealVar>(new RooRealVar((name+"_reso_sigma3").c_str(), (title+" reso_sigma3").c_str(), val, lo, hi)); };
	    void initResoFrac1(double val, double lo, double hi)
		{ reso_frac1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_reso_frac1").c_str(), (title+" reso_frac1").c_str(), val, lo, hi)); };
	    void initResoFrac2(double val, double lo, double hi)
		{ reso_frac2 = auto_ptr<RooRealVar>(new RooRealVar((name+"_reso_frac2").c_str(), (title+" reso_frac2").c_str(), val, lo, hi)); };
	    void initBgrMean(double val, double lo, double hi)
		{ bgr_mean = auto_ptr<RooRealVar>(new RooRealVar((name+"_bgr_mean").c_str(), (title+" bgr_mean").c_str(), val, lo, hi)); };
	    void initBgrSigma1(double val, double lo, double hi)
		{ bgr_sigma1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_bgr_sigma1").c_str(), (title+" bgr_sigma1").c_str(), val, lo, hi)); };
	    void initBgrSigma2(double val, double lo, double hi)
		{ bgr_sigma2 = auto_ptr<RooRealVar>(new RooRealVar((name+"_bgr_sigma2").c_str(), (title+" bgr_sigma2").c_str(), val, lo, hi)); };
	    void initBgrGaussFrac1(double val, double lo, double hi)
		{ bgr_frac1 = auto_ptr<RooRealVar>(new RooRealVar((name+"_bgr_frac1").c_str(), (title+" bgr_frac1").c_str(), val, lo, hi)); };
	    void initSigFrac(double val, double lo, double hi)
		{ sig_frac = auto_ptr<RooRealVar>(new RooRealVar((name+"_sig_frac").c_str(), (title+" sig_frac").c_str(), val, lo, hi)); };
	    void initModel()
	    {
		reso_gm1 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_reso_gm1").c_str(), (title+" Reso Gauss 1").c_str(), *var, *reso_mean.get(), *reso_sigma1.get()));
		reso_gm2 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_reso_gm2").c_str(), (title+" Reso Gauss 2").c_str(), *var, *reso_mean.get(), *reso_sigma2.get()));
		reso_gm3 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_reso_gm3").c_str(), (title+" Reso Gauss 3").c_str(), *var, *reso_mean.get(), *reso_sigma3.get()));
		reso_model = auto_ptr<RooAddModel>(new RooAddModel((name+"_reso_model").c_str(), (title+" resolution model").c_str(),
						    RooArgList(*reso_gm1.get(), *reso_gm2.get(), *reso_gm3.get()), RooArgList(*reso_frac1, *reso_frac2)));
		sig_model = auto_ptr<RooDecay>(new RooDecay((name+"_model_dec_trigauss").c_str(), (title+" model (decay plus triple Gauss)").c_str(), *var, *tau.get(),
						    *reso_model.get(), RooDecay::SingleSided));

		bgr_gm1 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_bgr_gm1").c_str(), (title+" Bgr Gauss 1").c_str(), *var, *bgr_mean.get(), *bgr_sigma1.get()));
		bgr_gm2 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_bgr_gm2").c_str(), (title+" Bgr Gauss 2").c_str(), *var, *bgr_mean.get(), *bgr_sigma2.get()));
		bgr_model = auto_ptr<RooAddModel>(new RooAddModel((name+"_bgr_model").c_str(), (title+" Bgr model (double Gauss)").c_str(), RooArgList(*bgr_gm1.get(), *bgr_gm2.get()),
				RooArgList(*bgr_frac1)));

		full_model = auto_ptr<RooAddPdf>(new RooAddPdf((name+"_full_model").c_str(), (title+" model (decay with trigauss reso fct and double Gauss bgr)").c_str(), 
			    RooArgList(*sig_model.get(), *bgr_model.get()), RooArgList(*sig_frac)));
	    };
	    void setConstant(bool t, bool m, bool s1, bool s2, bool s3, bool f1, bool f2)
	    {
		tau->setConstant(t);
		reso_mean->setConstant(m);
		reso_sigma1->setConstant(s1);
		reso_sigma2->setConstant(s2);
		reso_sigma3->setConstant(s3);
		reso_frac1->setConstant(f1);
		reso_frac2->setConstant(f2);
	    };
	    void setConstantBgr(bool m, bool s1, bool s2, bool f1)
	    {
		bgr_mean->setConstant(m);
		bgr_sigma1->setConstant(s1);
		bgr_sigma2->setConstant(s2);
		bgr_frac1->setConstant(f1);
	    };
	    void setConstantSigFrac(bool sf) { sig_frac->setConstant(sf); };
	    void plotOnFrame(RooPlot* frame, int nCPU)
	    {	// make sure that the full model is plotted as last function so that chiSquare reports the correct value afterwards
		full_model->plotOn(frame, RooFit::Components((name+"_model_dec_trigauss").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
		full_model->plotOn(frame, RooFit::Components((name+"_bgr_model").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(kDashed),RooFit::LineColor(kCyan));
		full_model->plotOn(frame, RooFit::NumCPU(nCPU));
	    };
	    void writeResultsOnFrame(RooPlot* frame, double scale = 1, std::string unit = "")
	    {
		frameWriter fw(frame, .65, .80, .05, .04);
		fw.addLine("#tau: " + roundToString(tau->getVal()*scale,3) + " #pm " + roundToString(tau->getError()*scale,3) + " " + unit);
		fw.textsize *= 0.8; fw.txtLineSpace *= 0.8;
		fw.addLine("");
		fw.addLine("Resolution function:");
		fw.addLine("mean: " + roundToString(reso_mean->getVal()*scale,3)      + " #pm " + roundToString(reso_mean->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 1: " + roundToString(reso_sigma1->getVal()*scale,3) + " #pm " + roundToString(reso_sigma1->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 2: " + roundToString(reso_sigma2->getVal()*scale,3) + " #pm " + roundToString(reso_sigma2->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 3: " + roundToString(reso_sigma3->getVal()*scale,3) + " #pm " + roundToString(reso_sigma3->getError()*scale,3) + " " + unit);
		fw.addLine("frac 1: " + roundToString(reso_frac1->getVal(),3)         + " #pm " + roundToString(reso_frac1->getError(),3));
		fw.addLine("frac 2: " + roundToString(reso_frac2->getVal(),3)         + " #pm " + roundToString(reso_frac2->getError(),3));
		const double restVal = 1 - (reso_frac1->getVal()+reso_frac2->getVal());
		const double restErr = sqrt(reso_frac1->getError()*reso_frac1->getError()+reso_frac2->getError()*reso_frac2->getError());
		fw.addLine("frac 3 (rest): " + roundToString(restVal,3) + " #pm " + roundToString(restErr,3));
		fw.addLine("Background function:");
		fw.addLine("mean: " + roundToString(bgr_mean->getVal()*scale,3)      + " #pm " + roundToString(bgr_mean->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 1: " + roundToString(bgr_sigma1->getVal()*scale,3) + " #pm " + roundToString(bgr_sigma1->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 2: " + roundToString(bgr_sigma2->getVal()*scale,3) + " #pm " + roundToString(bgr_sigma2->getError()*scale,3) + " " + unit);
		fw.addLine("frac 1: " + roundToString(bgr_frac1->getVal(),3)         + " #pm " + roundToString(bgr_frac1->getError(),3));
		fw.addLine("sig frac: " + roundToString(sig_frac->getVal(),3)         + " #pm " + roundToString(sig_frac->getError(),3));
		//fw.addLine("#chi^{2}/ndof: " + roundToString(frame->chiSquare(),3));
	    };
	    void markWorkInProgress(RooPlot* frame)
	    {
		frameWriter fw(frame, .40, .85, .05, .06);
		fw.addLine("Work in progress");
	    }
	};

	/*
	struct modelDecayTripleGaussBgrDoubleGaussDUMMY : public virtual modelDecayTripleGauss, modelDoubleGauss
	{
	    std::auto_ptr<RooRealVar> sigFrac;
	    std::auto_ptr<RooAddModel> model;
	    void initSigFrac(double val, double lo, double hi)
		{ sigFrac = auto_ptr<RooRealVar>(new RooRealVar((name+"_sigFrac").c_str(), (title+" signal fraction").c_str(), val, lo, hi)); };
	    void initSigFrac(double val, double lo, double hi, double err)
		{ initBgrFrac(val, lo, hi); sigFrac->setError(err); };
	    void initBgrFrac(double val, double lo, double hi) { initBgrFrac(1.0-val, 1.0-hi, 1.0-lo); };
	    void initModel()
	    {
		modelDecayTripleGauss::initModel();
		modelDoubleGauss::initModel();
		model = auto_ptr<RooAddModel>(new RooAddModel((name+"_model_dec_resotrigauss_bgrdblgauss").c_str(), (title+" model").c_str(), 
			    RooArgList(*modelDecayTripleGauss::model.get(), *modelDoubleGauss::model.get()),
			    RooArgList(*sigFrac.get())));
	    };
	    void setFracConstant(bool f) { sigFrac->setConstant(f); };
	    void plotOnFrame(RooPlot* frame, int nCPU)
	    {	// make sure that the full model is plotted as last function so that chiSquare reports the correct value afterwards
		model->plotOn(frame, RooFit::NumCPU(nCPU));
		model->plotOn(frame, RooFit::Components((name+"_model_dec_trigauss").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(2),RooFit::LineColor(kGreen));
		model->plotOn(frame, RooFit::Components((name+"_model_dblgauss").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(2),RooFit::LineColor(kCyan));
	    };
	    void writeResultsOnFrame(RooPlot* frame, double scale = 1, std::string unit = "")
	    {
		frameWriter fw;
		fw.setFrame(frame);
		fw.addLine("#tau: " + roundToString(tau->getVal()*scale,3) + " #pm " + roundToString(tau->getError()*scale,3) + " " + unit);
		fw.textsize *= 0.8; fw.txtLineSpace *= 0.8;
		fw.addLine("");
		fw.addLine("Resolution function:");
		fw.addLine("mean: " + roundToString(modelDecayTripleGauss::mean->getVal()*scale,3) + " #pm " + roundToString(modelDecayTripleGauss::mean->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 1: " + roundToString(modelDecayTripleGauss::sigma1->getVal()*scale,3) + " #pm " + roundToString(modelDecayTripleGauss::sigma1->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 2: " + roundToString(modelDecayTripleGauss::sigma2->getVal()*scale,3) + " #pm " + roundToString(modelDecayTripleGauss::sigma2->getError()*scale,3) + " " + unit);
		fw.addLine("sigma 3: " + roundToString(modelDecayTripleGauss::sigma3->getVal()*scale,3) + " #pm " + roundToString(modelDecayTripleGauss::sigma3->getError()*scale,3) + " " + unit);
		fw.addLine("frac 1: " + roundToString(modelDecayTripleGauss::frac1->getVal(),3) + " #pm " + roundToString(modelDecayTripleGauss::frac1->getError(),3));
		fw.addLine("frac 2: " + roundToString(modelDecayTripleGauss::frac2->getVal(),3) + " #pm " + roundToString(modelDecayTripleGauss::frac2->getError(),3));
		const double restVal = 1 - (modelDecayTripleGauss::frac1->getVal()+modelDecayTripleGauss::frac2->getVal());
		const double restErr = sqrt(modelDecayTripleGauss::frac1->getError()*modelDecayTripleGauss::frac1->getError()+modelDecayTripleGauss::frac2->getError()*modelDecayTripleGauss::frac2->getError());
		fw.addLine("frac 3 (rest): " + roundToString(restVal,3) + " #pm " + roundToString(restErr,3));
		fw.addLine("#chi^{2}/ndof: " + roundToString(frame->chiSquare(),3));
	    };
	};
	*/

	struct modelDoubleGaussWithFlatBgr : public modelDoubleGauss
	{
	    std::auto_ptr<RooRealVar> c0;
	    std::auto_ptr<RooAddModel> sig;
	    std::auto_ptr<RooChebychev> bgr;
	    std::auto_ptr<RooRealVar> Nsig, Nbgr;
	    Double_t Nsig_sigmas, Nsig_val, Nsig_err, Nbgr_val, Nbgr_err;
	    void initc0(double val, double lo, double hi)
		{ c0 = auto_ptr<RooRealVar>(new RooRealVar((name+"_c0").c_str(), (title+" c0").c_str(), val, lo, hi)); };
	    void initNsig(double val, double lo, double hi)
		{ Nsig = auto_ptr<RooRealVar>(new RooRealVar((name+"_Nsig").c_str(), (title+" Nsig").c_str(), val, lo, hi)); };
	    void initNbgr(double val, double lo, double hi)
		{ Nbgr = auto_ptr<RooRealVar>(new RooRealVar((name+"_Nbgr").c_str(), (title+" Nbgr").c_str(), val, lo, hi)); };
	    void initModel()
	    {
		gm1 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm1").c_str(), (title+" Gauss 1").c_str(), *var, *mean.get(), *sigma1.get()));
		gm2 = auto_ptr<RooGaussModel>(new RooGaussModel((name+"_gm2").c_str(), (title+" Gauss 2").c_str(), *var, *mean.get(), *sigma2.get()));
		sig = auto_ptr<RooAddModel>(new RooAddModel((name+"_sig").c_str(), (title+" Signal shape").c_str(), RooArgList(*gm1.get(), *gm2.get()),
						RooArgList(*frac1)));
		bgr = auto_ptr<RooChebychev>(new RooChebychev((name+"_bgr").c_str(), (title+" Bgr shape").c_str(), *var, RooArgSet(*c0.get())));
		model = auto_ptr<RooAddModel>(new RooAddModel((name+"_model").c_str(), (title+" model").c_str(), RooArgList(*sig.get(),*bgr.get()),
						RooArgList(*Nsig, *Nbgr)));
	    };
	    void setConstant(bool m, bool s1, bool s2, bool f1, bool c, bool ns, bool nb)
	    {
		mean->setConstant(m);
		sigma1->setConstant(s1);
		sigma2->setConstant(s2);
		frac1->setConstant(f1);
		c0->setConstant(c);
		Nsig->setConstant(ns);
		Nbgr->setConstant(nb);
	    };
	    void plotOnFrame(RooPlot* frame, int nCPU)
	    {	// make sure that the full model is plotted as last function so that chiSquare reports the correct value afterwards
		model->plotOn(frame, RooFit::Components((name+"_gm1").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(2),RooFit::LineColor(kGreen));
		model->plotOn(frame, RooFit::Components((name+"_gm2").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(2),RooFit::LineColor(kCyan));
		model->plotOn(frame, RooFit::Components((name+"_sig").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(7),RooFit::LineColor(kBlue));
		model->plotOn(frame, RooFit::Components((name+"_bgr").c_str()),RooFit::NumCPU(nCPU),RooFit::LineStyle(7),RooFit::LineColor(kBlue));
		model->plotOn(frame, RooFit::NumCPU(nCPU),RooFit::LineColor(kBlue));
	    };
	    void calculateRange(Double_t sigmas) // to force that this will be calculated with sig only
	    {
		sigSigmas = sigmas;
		const Double_t area = .5*(TMath::Erf(-sigSigmas/TMath::Sqrt2())+1);
		std::auto_ptr<RooAbsReal> cdf (sig->createCdf(*var));
		sigLo = cdf->findRoot(*var, var->getMin(), var->getMax(), area);
		sigHi = cdf->findRoot(*var, var->getMin(), var->getMax(), 1-area);
		cout << "sigLo: " << sigLo << " sigHi: " << sigHi << " \"sigma\": " << (sigHi-sigLo)/(2*sigSigmas) << " center: " << .5*(sigHi+sigLo) << endl;
	    };
	    void calculateSigInRange(Double_t sigmas)
	    {
		calculateRange(sigmas);
		std::string rangename = name + "_range_" + toString(sigmas) + "sigmas";
		var->setRange(rangename.c_str(), sigLo, sigHi);
		auto_ptr<RooAbsReal> nsigWindow(sig->createIntegral(*var, *var, rangename.c_str()));
		auto_ptr<RooAbsReal> nbgrWindow(bgr->createIntegral(*var, *var, rangename.c_str()));
		Nsig_val = nsigWindow->getVal()*Nsig->getVal();
		Nsig_err = nsigWindow->getVal()*Nsig->getError();
		Nbgr_val = nbgrWindow->getVal()*Nbgr->getVal();
		Nbgr_err = nbgrWindow->getVal()*Nbgr->getError();
		Nsig_sigmas = sigmas;
	    };
	    void writeResultsOnFrame(RooPlot* frame, double scale = 1, std::string unit = "", double sigmaScale = -1, std::string sigmaUnit = "")
	    {
		if (sigmaScale < 0)
		{
		    sigmaScale = scale;
		    sigmaUnit = unit;
		}
		frameWriter fw;
		fw.setFrame(frame);
		fw.addLine("Double Gaussian fit:");
		fw.addLine("mean: " + roundToString(mean->getVal()*scale,3) + " #pm " + roundToString(mean->getError()*scale,3) + " " + unit);
		fw.textsize *= 0.8; fw.txtLineSpace *= 0.8;
		fw.addLine("");
		fw.addLine("sigma 1: " + roundToString(sigma1->getVal()*sigmaScale,3) + " #pm " + roundToString(sigma1->getError()*sigmaScale,3) + " " + sigmaUnit);
		fw.addLine("sigma 2: " + roundToString(sigma2->getVal()*sigmaScale,3) + " #pm " + roundToString(sigma2->getError()*sigmaScale,3) + " " + sigmaUnit);
		fw.addLine("frac 1: " + roundToString(frac1->getVal(),3) + " #pm " + roundToString(frac1->getError(),3));
		const double restVal = 1 - frac1->getVal();
		const double restErr = frac1->getError();
		fw.addLine("frac 2 (rest): " + roundToString(restVal,3) + " #pm " + roundToString(restErr,3));
		fw.addLine("#chi^{2}/ndof: " + roundToString(frame->chiSquare(),3));
		fw.addLine("N_{sig}: " + roundToString(Nsig_val, 1) + " #pm " + roundToString(Nsig_err, 1));
		fw.addLine("N_{bgr}: " + roundToString(Nbgr_val, 1) + " #pm " + roundToString(Nbgr_err, 1));
		fw.addLine("(#pm" + roundToString(Nsig_sigmas, 2) + " #sigma  [" + roundToString(sigLo,3) + "," + roundToString(sigHi,3) + "])");
	    };
	};

	// ======================================================================================================= helper struct
	struct frameWriter
	{
	    frameWriter() : txtPosLeft(.60), txtPosTop(.80), txtLineSpace(.05), textsize(.04), line(0) {};
	    frameWriter(RooPlot* f, double l, double t, double sp, double sz)
		: frame(f), txtPosLeft(l), txtPosTop(t), txtLineSpace(sp), textsize(sz), line(0) {};

	    RooPlot* frame;
	    double txtPosLeft, txtPosTop, txtLineSpace, textsize;
	    int line;

	    void setFrame(RooPlot* f) { frame = f; };
	    void addLine(std::string text)
	    {
		frame->addObject(writeTLatex(text, txtPosLeft, txtPosTop-(line++)*txtLineSpace, textsize));
	    };
	};

	// ======================================================================================================= general methods
	TFile* openFile(std::string filename);
	TTree* getTreeFromFile(std::string treename, TFile* f);
	TTree* getTreeFromFile(std::string treename, std::string rename, TFile* f);
	TTree* openFileAndGetTree(std::string filename, std::string treename, TFile* &file);
	void openFileAndGetTree(DoFitAll01_dataset &ds);
	void openFileAndGetTree(DoFitAll01_dataset &ds, std::string rename);

	// ======================================================================================================= dataset initialisers
	void initMCB0datasets();
	void initDataB0datasets();
	void initDataLbdatasets();

	// ======================================================================================================= specific fitter methods
	void doFitMCB0reso();
	void doFitMCB0matchTau(bool floatingTau = true, bool floatingReso = true);
	void doFitMCB0matchTau(const modelTripleGauss &mtg, bool floatingTau = true, bool floatingReso = true);
	void doFitMCB0mass();
	void doFitDataB0mass();
	void doFitDataLbmass();

	// ======================================================================================================= fit models including variables
	modelTripleGauss mtgMCB0reso;
	modelDecayTripleGauss mdtgMCB0match;
	modelDoubleGauss mdgMCB0mass, mdgDataB0BgrT, mdgDataLbBgrT;
	modelDoubleGaussWithFlatBgr msgDataB0mass, msgDataLbmass;
	modelDecayTripleGaussBgrDoubleGauss mdtgbdgDataB0lifetime, mdtgbdgDataLblifetime;
	//modelVoigtian msgDataB0mass;

    private:
	// ======================================================================================================= specific fitter workers
	//void doFitMCB0resoWorker();
	void doFitMCB0matchTauWorker(const bool &floatingTau, const bool &floatingReso);

	// ------------------------------ datasets
	RooDataSet *dataMCB0_, *dataMCB0match_, *dataMCB0nomatch_;
	RooDataSet *dataMCLb_, *dataDatB0_, *dataDatLb_;

	// ------------------------------ dataset containers
	DoFitAll01_dataset ds_MCB0_, ds_MCB0match_, ds_MCB0nomatch_;
	DoFitAll01_dataset ds_dataB0barrel, ds_dataLbbarrel;

	// ------------------------------ observables from data: mass, time and time diff (MC only, resolution studies)
	RooRealVar *rrvMCB0_m_, *rrvMCB0_t_, *rrvMCB0_dt_; // MC B0
	RooRealVar *rrvMCLb_m_, *rrvMCLb_t_, *rrvMCLb_dt_; // MC Lb
	RooRealVar *rrvDatB0_m_, *rrvDatB0_t_;             // data B0
	RooRealVar *rrvDatLb_m_, *rrvDatLb_t_;             // data Lb

	// ------------------------------ fit variables

	// ------------------------------ frames
	std::auto_ptr<RooPlot> curFrame;
	std::auto_ptr<TCanvas> curCanvas;

	// ------------------------------ steering variables
	int fnCPU_;
	int fVerbose_;
	int fnBins_t_, fnBins_m_;
	bool isPrelim, isWorkInProgress, noTitle;
};

#endif
