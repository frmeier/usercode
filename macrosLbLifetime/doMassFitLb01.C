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
#include "TLine.h"
#include "TArrow.h"
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
using RooFit::LineColor;

int fnCPU_Lb_(4);

struct doMassFitLb01_fitresults
{
    double sig;
    double bgr;
    double soversqrtsb;
    double soverb;
    double mass;
    double width;
    double twoSigmas, threeSigmas;
};

// Unbinned likelihood fit for Lb mass
doMassFitLb01_fitresults doMassFitLb01(TTree * tree, const double lumi, const string addTitle = "", const bool noTitle = false, const bool prelim = true, const bool noStat = false)
{
    cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-" << endl;
    cout << "This is doMassFitLb01 -- fitting the mass of Lb to the passed tree" << endl;

    const double curEntries = tree->GetEntries();
    cout << "curEntries: " << curEntries << endl;

    // main configs are here
    const double valLo(5.35), valHi(5.895);   // range to cover for fit and histo
    const double binwidth(0.005);             // bin width for histo (not fit)
    const int nBins = (valHi-valLo)/binwidth; // width of the bins
    const double sigmas(2.0);                 // region around peak to use for signal estimation

    //const double valLo(5.35), valHi(5.90);
    //const int nBins(33);

    RooRealVar mass("mbc", "J/#psi #Lambda mass [GeV/c^{2}]", valLo, valHi);
    RooDataSet data("data", "Lambda_b dataset", tree, mass);

    //   mass fit
    RooRealVar mean("mean","mean",5.62,5.60,5.64);
    RooRealVar sigma("sigma","sigma",0.014,0.01,0.03);
    RooGaussian sig("sig","signal p.d.f.",mass,mean,sigma);

    RooRealVar c0("c0","coefficient #0", 1.0,-1,2);
    RooRealVar c1("c1","coefficient #1", 0.1,-1,1);
    RooChebychev bkg("bkg","background p.d.f.",mass,RooArgSet(c0));

    RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries);
    RooRealVar nbkg("nbkg","Background fraction", .01* curEntries,0.,curEntries*1.4);

    RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));
    model.fitTo(data,Extended(kTRUE), RooFit::NumCPU(fnCPU_Lb_));

    RooPlot* xframe = mass.frame(Name("xframe"),Title(("Mass of #Lambda_{b}"+ (addTitle.size()>0 ? " - " + addTitle : "")).c_str()));
    if(noTitle) xframe->SetTitle("");
    data.plotOn(xframe,Binning(nBins));
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
    if (prelim) xframe->addObject(writeTLatex("Work in progress",.18,.85,textsize));
    xframe->addObject(writeTLatex("#sqrt{s} = 7 TeV",.18,.77,textsize));
    xframe->Draw();
    //delete xframe;
    //delete fracSigRange;
    //delete fracBGRange;

    // preparing results for return to caller
    doMassFitLb01_fitresults res;
    res.sig = nsigWindow;
    res.bgr = nbkgWindow;
    res.soversqrtsb = nsigWindow/sqrt(nsigWindow+nbkgWindow);
    res.soverb = nsigWindow/nbkgWindow;
    res.mass = m;
    res.width = s1;
    return res;
}

// Unbinned likelihood fit for B0 mass using a double Gaussian
doMassFitLb01_fitresults doMassFitLb01_DG(TTree * tree, const double lumi, const string addTitle = "", const bool noTitle = false, const bool prelim = true, const bool noStat = false, const bool drawSigmaRanges = false, const string leafname = "mbc")
{
    doMassFitLb01_fitresults res; // container for results

    cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-" << endl;
    cout << "This is doMassFitLb01 -- fitting the mass of Lb to the passed tree with a double Gaussian" << endl;

    const double curEntries = tree->GetEntries();
    cout << "curEntries: " << curEntries << endl;

    // main configs are here
    const double valLo(5.35), valHi(5.895);   // range to cover for fit and histo
    const double binwidth(0.005);             // bin width for histo (not fit)
    const int nBins = (valHi-valLo)/binwidth; // width of the bins
    //const double sigmas(2.5);                 // region around peak to use for signal estimation

    //const double valLo(5.35), valHi(5.90);
    //const int nBins(33);

    RooRealVar mass(leafname.c_str(), "J/#psi #Lambda mass [GeV/c^{2}]", valLo, valHi);
    RooDataSet data("data", "B0 dataset", tree, mass);

    //   mass fit
    RooRealVar mean("mean","mean",5.62,5.60,5.64);
    RooRealVar sigma1("sigma1","sigma1",0.010,0.005,0.015);
    RooRealVar sigma2("sigma2","sigma2",0.030,0.010,0.050);
    RooGaussian sig1("sig1","signal Gaussian 1",mass,mean,sigma1);
    RooGaussian sig2("sig2","signal Gaussian 2",mass,mean,sigma2);
    RooRealVar sig1frac("sig1frac","gauss1 fraction", .4, .1, .9);
    RooAddPdf sig("sig","double gaussian",RooArgList(sig1,sig2),RooArgList(sig1frac));

    RooRealVar c0("c0","coefficient #0", 1.0,-1,2);
    RooRealVar c1("c1","coefficient #1", 0.1,-1,1);
    RooChebychev bkg("bkg","background p.d.f.",mass,RooArgSet(c0));

    RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries);
    RooRealVar nbkg("nbkg","Background fraction", .01* curEntries,0.,curEntries*1.4);

    RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));
    model.fitTo(data,Extended(kTRUE),RooFit::NumCPU(fnCPU_Lb_));

    RooPlot* xframe = mass.frame(Name("xframe"),Title(("Mass of #Lambda_{b}"+ (addTitle.size()>0 ? " - " + addTitle : "")).c_str()));
    if(noTitle) xframe->SetTitle("");
    data.plotOn(xframe,Binning(nBins));
    model.plotOn(xframe);
    model.plotOn(xframe,Components("bkg"),LineStyle(kDashed));
    model.plotOn(xframe,Components("sig1"),LineStyle(kDashed),LineColor(kGreen));
    model.plotOn(xframe,Components("sig2"),LineStyle(kDashed),LineColor(kGreen));
    
    xframe->SetYTitle(("Events / (" + toString(binwidth) + " GeV/c^{2})").c_str());
    xframe->Draw();

    // get range of 2 and 3 sigmas
    const Double_t area2s = .5*(TMath::Erf(-2.0/TMath::Sqrt2())+1); // one sided
    const Double_t area3s = .5*(TMath::Erf(-3.0/TMath::Sqrt2())+1);
    const Double_t areafrac2s = 1.0 - 2*area2s; // two sided

    RooAbsReal* cdf = sig.createCdf(mass);
    const Double_t mlo = cdf->findRoot(mass, mass.getMin(), mass.getMax(), area2s);
    const Double_t mhi = cdf->findRoot(mass, mass.getMin(), mass.getMax(), 1.0-area2s);
    res.twoSigmas = (mhi-mlo)*0.5;
    const Double_t mlo3 = cdf->findRoot(mass, mass.getMin(), mass.getMax(), area3s);
    const Double_t mhi3 = cdf->findRoot(mass, mass.getMin(), mass.getMax(), 1.0-area3s);
    res.threeSigmas = (mhi3-mlo3)*0.5;

    // extract nsig and ngr
    mass.setRange("sigwindow", mlo, mhi);
    RooAbsReal* fracSigRange = sig.createIntegral(mass,mass,"sigwindow");
    const Double_t nsig_val = fracSigRange->getVal()*nsig.getVal();
    const Double_t nsig_err = nsig.getError()*fracSigRange->getVal();
    RooAbsReal* fracBGRange = bkg.createIntegral(mass,mass,"sigwindow");
    const Double_t nbkg_val = nbkg.getVal()*fracBGRange->getVal();
    const Double_t nbkg_err = nbkg.getError()*fracBGRange->getVal();
    cout << "nbkg.getVal(): " << nbkg.getVal() << " nbkg_val: " << nbkg_val << endl;
    cout << "nbkg.getError(): " << nbkg.getError() << " nbkg_err: " << nbkg_err << endl;
    cout << "fracBGRange->getVal(): " << fracBGRange->getVal() << endl;

    cout << "S/Sqrt(S+B)  =  " << nsig_val/sqrt(nsig_val+nbkg_val) << endl;

    double sb_err = sqrt(nsig_err*nsig_err/(nsig_val*nsig_val)+nbkg_err*nbkg_err/(nbkg_val*nbkg_val))*nsig_val/nbkg_val;
    double s_over_sqrt_err = sqrt(nsig_err*nsig_err/(nsig_val*nsig_val)+(nbkg_err*nbkg_err+nsig_err*nsig_err)/(4.0*(nbkg_val+nsig_val)*(nbkg_val+nsig_val)))*nsig_val/sqrt(nsig_val+nbkg_val);

    // write the results to the canvas
    const double txtPosLeft = .6;
    const double txtPosTop = (lumi>0 ? .70 : .80);
    const double txtLineSpace = .05;
    const double textsize = 0.06;

    if (!noStat)
    {
	int line(0);
	xframe->addObject(writeTLatex("n Signal: " + roundToString(nsig_val,1)+" #pm "+roundToString(nsig_err,1)+" (2#sigma window)",txtPosLeft,txtPosTop-(line++)*txtLineSpace));
	xframe->addObject(writeTLatex("n Bgr: " + roundToString(nbkg_val,1)+" #pm "+roundToString(nbkg_err,1)+" (2#sigma window)", txtPosLeft,txtPosTop-(line++)*txtLineSpace));
	xframe->addObject(writeTLatex("S/#sqrt{S+B}: " + roundToString(nsig_val/sqrt(nsig_val+nbkg_val),1) + " #pm " + roundToString(s_over_sqrt_err,1) ,txtPosLeft,txtPosTop-(line++)*txtLineSpace));
	xframe->addObject(writeTLatex("S/B: " + roundToString(nsig_val/nbkg_val,2) + " #pm " + roundToString(sb_err,2),txtPosLeft,txtPosTop-(line++)*txtLineSpace));
	const Double_t m_err=mean.getError();
	xframe->addObject(writeTLatex("mass: " + roundToString(mean.getVal(),m_err>0.003?3:4) + " #pm " + roundToString(m_err,m_err>0.003?3:4)+" GeV/c^{2}", txtPosLeft,txtPosTop-(line++)*txtLineSpace));
	const Double_t s1_err=sigma1.getError();
	xframe->addObject(writeTLatex("width 1: " + roundToString(sigma1.getVal(),s1_err>0.003?3:4) + " #pm " + roundToString(s1_err,s1_err>0.003?3:4)+ " GeV/c^{2}",txtPosLeft,txtPosTop-(line++)*txtLineSpace));
	const Double_t s2_err=sigma2.getError();
	xframe->addObject(writeTLatex("width 2: " + roundToString(sigma2.getVal(),s2_err>0.003?3:4) + " #pm " + roundToString(s2_err,s2_err>0.003?3:4)+ " GeV/c^{2}",txtPosLeft,txtPosTop-(line++)*txtLineSpace));
    }
    if (lumi>0) xframe->addObject(writeTLatex("#int Ldt = " + toString(lumi) + " pb^{-1}",txtPosLeft,.80,textsize));
    if (prelim) xframe->addObject(writeTLatex("Work in progress",.18,.85,textsize));
    xframe->addObject(writeTLatex("#sqrt{s} = 7 TeV",.18,.77,textsize));
    xframe->Draw();

    // draw sigma ranges, if required
    if (drawSigmaRanges)
    {
	const Double_t hfrac = 1.0/(2.0*nsig_val/nbkg_val+1.0);
	const Double_t min = xframe->GetMinimum();
	const Double_t max = xframe->GetMaximum();
	TLine *l_left2 = new TLine(mlo, min, mlo, max*hfrac*.6);
	l_left2->Draw();
	TLine *l_right2 = new TLine(mhi, min, mhi, max*hfrac*.6);
	l_right2->Draw();
	TArrow *arr2 = new TArrow(mhi, max*hfrac*.55, mlo, max*hfrac*.55, .015, "<>");
	arr2->Draw();
	TLatex* txt2 = new TLatex(mass.getVal(), max*hfrac*.55, "2#sigma equiv.");
	txt2->SetTextAlign(21);
	txt2->SetTextSize(0.035);
	txt2->Draw();
	TLine *l_left3 = new TLine(mlo3, min, mlo3, max*hfrac*.45);
	l_left3->Draw();
	TLine *l_right3 = new TLine(mhi3, min, mhi3, max*hfrac*.45);
	l_right3->Draw();
	TArrow *arr3 = new TArrow(mhi3, max*hfrac*.4, mlo3, max*hfrac*.4, .015, "<>");
	arr3->Draw();
	TLatex* txt3 = new TLatex(mass.getVal(), max*hfrac*.38, "3#sigma equiv.");
	txt3->SetTextAlign(23);
	txt3->SetTextSize(0.035);
	txt3->Draw();
    }

    // preparing results for return to caller
    res.sig = nsig_val;
    //res.sig_err = nsig_err;
    res.bgr = nbkg_val;
    res.soversqrtsb = nsig_val/sqrt(nsig_val+nbkg_val);
    res.soverb = nsig_val/nbkg_val;
    res.mass = mean.getVal();
    res.width = sigma1.getVal();
    return res;
}

