/////////////////////////////////////////////////////////////////////////
//
// 'VALIDATION AND MC STUDIES' RooFit tutorial macro #801
// 
// A Toy Monte Carlo study that perform cycles of
// event generation and fittting
//
// 
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH2.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "setTDRStyle_modified.C"

using namespace RooFit ;

// based on rf801_mcstudy
void rooFitValid()
{
  // C r e a t e   m o d e l
  // -----------------------

    // some config
    
	const double curEntries = 600;
	const double valLo(5.22), valHi(6.045);

  // Declare observable x
        RooRealVar mass("mlb", "J/#psi #Lambda mass [GeV]", valLo, valHi);
	mass.setBins(33);

  // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their paramaters
        RooRealVar mean("mean","mean",5.62,5.60,5.64) ;
        RooRealVar sigma("sigma","sigma",0.02,0.01,0.03) ;
        RooGaussian sig("sig","signal p.d.f.",mass,mean,sigma) ;
  
        RooRealVar c0("c0","coefficient #0", 1.0,-1,2) ;
        RooRealVar c1("c1","coefficient #1", 0.1,-1,1) ;
        RooPolynomial bkg("bkg","background p.d.f.",mass,RooArgList(c0,c1),0) ;

        //RooRealVar nsig("nsig","signal fraction", .833*curEntries,0.,curEntries) ;
        //RooRealVar nbkg("nbkg","Background fraction", .167* curEntries) ;
        RooRealVar nsig("nsig","signal fraction", .833 * curEntries, 0, curEntries) ;
        RooRealVar nbkg("nbkg","Background fraction", .167* curEntries, 0, curEntries) ;

        RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg)) ;



  // C r e a t e   m a n a g e r
  // ---------------------------

  // Instantiate RooMCStudy manager on model with x as observable and given choice of fit options
  //
  // The Silence() option kills all messages below the PROGRESS level, leaving only a single message
  // per sample executed, and any error message that occur during fitting
  //
  // The Extended() option has two effects: 
  //    1) The extended ML term is included in the likelihood and 
  //    2) A poisson fluctuation is introduced on the number of generated events 
  //
  // The FitOptions() given here are passed to the fitting stage of each toy experiment.
  // If Save() is specified, the fit result of each experiment is saved by the manager  
  //
  // A Binned() option is added in this example to bin the data between generation and fitting
  // to speed up the study at the expemse of some precision

  RooMCStudy* mcstudy = new RooMCStudy(model,mass,Binned(kFALSE),Silence(),Extended(),
				       FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;
  

  // G e n e r a t e   a n d   f i t   e v e n t s
  // ---------------------------------------------

  // Generate and fit 1000 samples of Poisson(nExpected) events
  mcstudy->generateAndFit(1000) ;



  // E x p l o r e   r e s u l t s   o f   s t u d y 
  // ------------------------------------------------

  //setTDRStyle();
  // Make plots of the distributions of mean, the error on mean and the pull of mean
  RooPlot* frame1 = mcstudy->plotParam(mean,Bins(40)) ;
  RooPlot* frame2 = mcstudy->plotError(mean,Bins(40)) ;
  RooPlot* frame3 = mcstudy->plotPull(mean,Bins(40),FitGauss(kTRUE)) ;

  // Plot distribution of minimized likelihood
  RooPlot* frame4 = mcstudy->plotNLL(Bins(40)) ;

  // Make some histograms from the parameter dataset
  TH1* hh_cor_a0_s1f = mcstudy->fitParDataSet().createHistogram("hh",c1,YVar(mean)) ;
  TH1* hh_cor_a0_a1  = mcstudy->fitParDataSet().createHistogram("hh",c0,YVar(c1)) ;

  // Access some of the saved fit results from individual toys
  TH2* corrHist000 = mcstudy->fitResult(0)->correlationHist("c000") ;
  TH2* corrHist127 = mcstudy->fitResult(127)->correlationHist("c127") ;
  TH2* corrHist953 = mcstudy->fitResult(953)->correlationHist("c953") ;



  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* c = new TCanvas("rooFitValid","rooFitValid",900,900) ;
  c->Divide(3,3) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
  c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
  c->cd(5) ; gPad->SetLeftMargin(0.15) ; hh_cor_a0_s1f->GetYaxis()->SetTitleOffset(1.4) ; hh_cor_a0_s1f->Draw("box") ;
  c->cd(6) ; gPad->SetLeftMargin(0.15) ; hh_cor_a0_a1->GetYaxis()->SetTitleOffset(1.4) ; hh_cor_a0_a1->Draw("box") ;
  c->cd(7) ; gPad->SetLeftMargin(0.15) ; corrHist000->GetYaxis()->SetTitleOffset(1.4) ; corrHist000->Draw("colz") ;
  c->cd(8) ; gPad->SetLeftMargin(0.15) ; corrHist127->GetYaxis()->SetTitleOffset(1.4) ; corrHist127->Draw("colz") ;
  c->cd(9) ; gPad->SetLeftMargin(0.15) ; corrHist953->GetYaxis()->SetTitleOffset(1.4) ; corrHist953->Draw("colz") ;

  // Make RooMCStudy object available on command line after
  // macro finishes
  gDirectory->Add(mcstudy) ;
}
