
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
    
	const double valLo(5.22), valHi(6.045);


  // Declare observable x
        RooRealVar mass("mlb", "J/#psi #Lambda mass [GeV/c^{2}]", valLo, valHi);
	//	mass.setBins(33);

        //   mass fit
        double m=5.62354;
	double m_sig=0.00360244;
        double sg=0.0215277;
	double sg_sig=0.0031824;
	double ns=71.2533;
	double ns_sig=10.8354;
	double nb=255.756;
        double nb_sig=17.3733;
        double cc=-0.335304;
        double cc_sig=0.109179;
	double n_std = 3.0;
        RooRealVar mean("mean","mean",m,5.60,5.64);
        RooRealVar sigma("sigma","sigma",sg,0.01,0.03);
        RooGaussian sig("sig","signal p.d.f.",mass,mean,sigma) ;

        RooRealVar c0("c0","coefficient #0", cc,cc-n_std*cc_sig,cc+n_std*cc_sig) ;
        RooChebychev bkg("bkg","background p.d.f.",mass,RooArgSet(c0)) ;

        RooRealVar nsig("nsig","n_{sig}", ns, ns-n_std*ns_sig, ns+n_std*ns_sig);
        RooRealVar nbkg("nbkg","n_{bg}", nb, nb-n_std*nb_sig, nb+n_std*nb_sig);

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

  RooMCStudy* mcstudy = new RooMCStudy(model,mass,Binned(kFALSE),Silence(), Extended(),
				       FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;
  

  // G e n e r a t e   a n d   f i t   e v e n t s
  // ---------------------------------------------

  // Generate and fit 1000 samples of Poisson(nExpected) events
  mcstudy->generateAndFit(1000) ;


  // E x p l o r e   r e s u l t s   o f   s t u d y 
  // ------------------------------------------------

  setTDRStyle();
  // Make plots of the distributions of mean, the error on mean and the pull of mean
  RooPlot* frame1 = mcstudy->plotPull(nsig, FrameRange(-3.0,3.0),Bins(40),FitGauss(kTRUE)) ;
  RooPlot* frame2 = mcstudy->plotPull(nbkg, FrameRange(-3.0,3.0),Bins(40),FitGauss(kTRUE)) ;
  RooPlot* frame3 = mcstudy->plotPull(mean, FrameRange(-3.0,3.0),Bins(40),FitGauss(kTRUE)) ;
  RooPlot* frame4 = mcstudy->plotPull(sigma, FrameRange(-3.0,3.0),Bins(40),FitGauss(kTRUE)) ;

  RooPlot* frame5 = mcstudy->plotParam(nsig,Bins(40)) ;
  RooPlot* frame6 = mcstudy->plotParam(nbkg,Bins(40)) ;
  RooPlot* frame7 = mcstudy->plotParam(mean,Bins(40));
  RooPlot* frame8 = mcstudy->plotParam(sigma,Bins(40));
  RooPlot* frame9 = mcstudy->plotParam(c0,Bins(40));

  // Plot distribution of minimized likelihood
  RooPlot* frame10 = mcstudy->plotNLL(Bins(40)) ;



  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* c = new TCanvas("rooFitPulls","rooFitPulls",900,900) ;
  c->Divide(2,2) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->SetTitle(""); frame1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->SetTitle(""); frame2->Draw() ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->SetTitle(""); frame3->Draw() ;
  c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->SetTitle(""); frame4->Draw() ;

  TCanvas* c2 = new TCanvas("rooFitValid","rooFitValid",900,900) ;
  c2->Divide(3,2) ;
  c2->cd(1) ; gPad->SetLeftMargin(0.15) ; frame5->GetYaxis()->SetTitleOffset(1.4) ;frame5->SetTitle(""); frame5->Draw() ;
  c2->cd(2) ; gPad->SetLeftMargin(0.15) ; frame6->GetYaxis()->SetTitleOffset(1.4) ;frame6->SetTitle(""); frame6->Draw() ;
  c2->cd(3) ; gPad->SetLeftMargin(0.15) ; frame7->GetYaxis()->SetTitleOffset(1.4) ;frame7->SetTitle(""); frame7->Draw() ;
  c2->cd(4) ; gPad->SetLeftMargin(0.15) ; frame8->GetYaxis()->SetTitleOffset(1.4) ;frame8->SetTitle(""); frame8->Draw() ;
  c2->cd(5) ; gPad->SetLeftMargin(0.15) ; frame9->GetYaxis()->SetTitleOffset(1.4) ;frame9->SetTitle(""); frame9->Draw() ;
  c2->cd(6) ; gPad->SetLeftMargin(0.15) ; frame10->GetYaxis()->SetTitleOffset(1.4) ;frame10->SetTitle(""); frame10->Draw() ;



  // Make RooMCStudy object available on command line after
  // macro finishes
  gDirectory->Add(mcstudy) ;
}
