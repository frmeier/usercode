#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include "RooGaussModel.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "TMath.h"
#include <string>
#include "TLatex.h"
#include <sstream>
#include <iomanip>
//#include "setTDRStyle_modified.C"


using std::cout;
using std::endl;
using RooFit::Extended;
using RooFit::Name;
using RooFit::Title;
using RooFit::Binning;
using RooFit::Components;
using RooFit::LineStyle;
using namespace RooFit;

struct doCtFitB0_fitresults
{
    double lifetime;
    double lifetimeE;
    string title;
};

doCtFitB0_fitresults doCtFitB0_01(TTree * tree, const double lumi, const string addTitle = "", const bool noTitle = false, const bool prelim = true, const bool noStat = false)
{
    doCtFitB0_fitresults fitresults;

    const int nBins(80);

    RooRealVar t("ct3dB0","measured lifetime B^{0}",1.5e-12, 15e-12,"s");
    RooTruthModel idealres("idealres","Ideal resolution model",t);

    //RooRealVar tau("tau","tau",1.525e-12,1e-17,5e-12,"s");
    RooRealVar tau("tau","tau",-1./1.525e-12,-1./1.4e-12,-1./5e-12,"s-1");
    RooExponential decay("decay","decay",t,tau);
    //RooDecay decay("decay","decay",t,tau,idealres,RooDecay::SingleSided);

    // Getting data from tree
    RooDataSet data("data", "B^{0} lifetime dataset", t, RooFit::Import(*tree));

    // Now do the fit
    decay.fitTo(data);

    fitresults.lifetime = tau.getVal();
    fitresults.lifetimeE = tau.getError();
    cout << "tau: " << -1./tau.getVal() << endl;

    // Plotting
    RooPlot* frame = t.frame(Title("Lifetime of B^{0}"));
    data.plotOn(frame,RooFit::Binning(nBins));
    decay.plotOn(frame,ProjWData(data,true));
    
    frame->Draw();
    gPad->SetLogy();
    return fitresults;
}

