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
#include "TMath.h"
#include <string>
#include <sstream>
#include "TLatex.h"
#include <sstream>
#include <iomanip>
//#include "setTDRStyle_modified.C"
#include "utils.h"

using std::cout;
using std::endl;
using RooFit::Extended;
using RooFit::Name;
using RooFit::Title;
using RooFit::Binning;
using RooFit::Components;
using RooFit::LineStyle;

struct doMassFitJp01_fitresults
{
    double sig;
    double sig_err;
    double bgr;
    double soversqrtsb;
    double soverb;
    double mass;
    double width;
};

// Unbinned likelihood fit for Jp mass
doMassFitJp01_fitresults doMassFitJp01(TTree * tree, const double lumi, const string addTitle = "", const bool noTitle = false, const bool prelim = true, const bool noStat = false, const string toFit = "mjp")
{
    cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-" << endl;
    cout << "This is doMassFitJp01 -- fitting the mass of Jp to the passed tree" << endl;

    const double curEntries = tree->GetEntries();
    cout << "curEntries: " << curEntries << endl;

    // main configs are here
    const double valLo(2.9), valHi(3.3);   // range to cover for fit and histo
    const double binwidth(curEntries > 500 ? 0.005 : 0.01);             // bin width for histo (not fit)
    const int nBins = (valHi-valLo)/binwidth; // width of the bins
    const double sigmas(2.5);                 // region around peak to use for signal estimation

    //const double valLo(5.35), valHi(5.90);
    //const int nBins(33);

    //RooRealVar mass("mass", "#mu#mu mass [GeV/c^{2}]", valLo, valHi);
    RooRealVar mass(toFit.c_str(),"#mu#mu mass [GeV/c^{2}]", valLo, valHi);
    RooDataSet data("data", "dataset", tree, mass);

    //   mass fit
    RooRealVar mean("mean","mean",3.097,3.0,3.2) ;
    RooRealVar sigma("sigma","sigma",0.03,0.01,0.09) ;
    RooGaussian sig("sig","signal p.d.f.",mass,mean,sigma) ;

    RooRealVar c0("c0","coefficient #0", 1.0,-1,2) ;
    RooRealVar c1("c1","coefficient #1", 0.1,-1,1) ;
    RooChebychev bkg("bkg","background p.d.f.",mass,RooArgSet(c0)) ;

    RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries) ;
    RooRealVar nbkg("nbkg","Background fraction", .01* curEntries,0.,curEntries*1.4) ;

    RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg)) ;
    model.fitTo(data,Extended(kTRUE));

    RooPlot* xframe = mass.frame(Name("xframe"),Title(("Mass of J/#psi"+ (addTitle.size()>0 ? " - " + addTitle : "")).c_str())) ;
    if(noTitle) xframe->SetTitle("");
    data.plotOn(xframe,Binning(nBins)) ;
    model.plotOn(xframe);
    model.plotOn(xframe,Components("bkg"),LineStyle(kDashed));
    
    xframe->SetYTitle(("Events / (" + toString(binwidth) + " GeV/c^{2})").c_str());
    xframe->Draw();

    Double_t m=mean.getVal();
    Double_t m_err=mean.getError();
    Double_t s1=sigma.getVal();
    Double_t s1_err=sigma.getError();

    mass.setRange("window",m-sigmas*s1,m+sigmas*s1);
    RooAbsReal* fracSigRange = sig.createIntegral(mass,mass,"window");
    Double_t nsigWindow = fracSigRange->getVal()*nsig.getVal();
    Double_t nsig_err = nsig.getError()*fracSigRange->getVal();
    RooAbsReal* fracBGRange = bkg.createIntegral(mass,mass,"window");
    Double_t nbkgWindow = nbkg.getVal()*fracBGRange->getVal();
    Double_t nbkg_err = nbkg.getError()*fracBGRange->getVal();

    cout << "n_Signal     =  " << nsigWindow << endl;
    cout << "n_Background =  " << nbkgWindow << endl;
    cout << "S/Sqrt(S+B)  =  " << nsigWindow/sqrt(nsigWindow+nbkgWindow) << endl;

    double sb_err = sqrt(nsig_err*nsig_err/(nsigWindow*nsigWindow)+nbkg_err*nbkg_err/(nbkgWindow*nbkgWindow))*nsigWindow/nbkgWindow;
    double s_over_sqrt_err = sqrt(nsig_err*nsig_err/(nsigWindow*nsigWindow)+(nbkg_err*nbkg_err+nsig_err*nsig_err)/(4.0*(nbkgWindow+nsigWindow)*(nbkgWindow+nsigWindow)))*nsigWindow/sqrt(nsigWindow+nbkgWindow);

    // write the results to the canvas
    const double txtPosLeft = .6;
    const double txtPosTop = (lumi>0 ? .70 : .80);
    const double txtLineSpace = .05;
    const double textsize = 0.06;

    if (!noStat)
    {
	xframe->addObject(writeTLatex("n Signal: " + roundToString(nsigWindow,1)+" #pm "+roundToString(nsig_err,1),txtPosLeft,txtPosTop-0*txtLineSpace));
	xframe->addObject(writeTLatex("n Bgr: " + roundToString(nbkgWindow,1)+" #pm "+roundToString(nbkg_err,1), txtPosLeft,txtPosTop-1*txtLineSpace));
	xframe->addObject(writeTLatex("S/#sqrt{S+B}: " + roundToString(nsigWindow/sqrt(nsigWindow+nbkgWindow),1) + " #pm " + roundToString(s_over_sqrt_err,1) ,txtPosLeft,txtPosTop-2*txtLineSpace));
	xframe->addObject(writeTLatex("S/B: " + roundToString(nsigWindow/nbkgWindow,2) + " #pm " + roundToString(sb_err,2),txtPosLeft,txtPosTop-3*txtLineSpace));
	xframe->addObject(writeTLatex("mass: " + roundToString(m,m_err>0.003?3:4) + " #pm " + roundToString(m_err,m_err>0.003?3:4)+" GeV/c^{2}", txtPosLeft,txtPosTop-4*txtLineSpace));
	xframe->addObject(writeTLatex("width: " + roundToString(s1,s1_err>0.003?3:4) + " #pm " + roundToString(s1_err,s1_err>0.003?3:4)+ " GeV/c^{2}",txtPosLeft,txtPosTop-5*txtLineSpace));
    }
    if (lumi>0) xframe->addObject(writeTLatex("#int Ldt = " + toString(lumi) + " pb^{-1}",txtPosLeft,.80,textsize));
    if (prelim) xframe->addObject(writeTLatex("CMS preliminary",.18,.85,textsize));
    xframe->addObject(writeTLatex("#sqrt{s} = 7 TeV",.18,.77,textsize));
    xframe->Draw();

    // preparing results for return to caller
    doMassFitJp01_fitresults res;
    res.sig = nsigWindow;
    res.sig_err = nsig_err;
    res.bgr = nbkgWindow;
    res.soversqrtsb = nsigWindow/sqrt(nsigWindow+nbkgWindow);
    res.soverb = nsigWindow/nbkgWindow;
    res.mass = m;
    res.width = s1;
    return res;
}

