#include <iomanip>
#include <sstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
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
#include "RooPolynomial.h"
#include "RooTruthModel.h"
#include "RooEffProd.h"
#include "RooAddModel.h"
#include "TH2F.h"
#include "TLatex.h"

#include "utils.h"

using namespace RooFit;
using namespace std;

enum Types {Lambda_b, B0, Bs, Bplus};

const int set_nCPU(6);
const double timeFactor(1.0); // make it 1e-12 if your data is in s, or take 1 if it is in ps
const double invTimeFactor = 1./timeFactor;

TTree *tree;
ofstream result;
bool PerEventError;
bool one_scale;
int type;
bool noplot;
bool validation;
double mass_low;
double mass_high;
double masspeak;
TString title1,title2,title3;
TFile *f;

Double_t yield_nonprompt, yield_sig, yield_prompt;
Double_t two_sigma_upper, two_sigma_lower, three_sigma_upper, three_sigma_lower;
TCanvas *c1, *c2, *c3;
TString name_c1("M_vs_lt"), name_c2("projections"), name_c3("slices");
int Configure(TString arguments);
TLatex* writeTLatex(string text, double x, double y, double size = 0.04);
string roundToString(double v, streamsize precision);


void Fit_2D(TString type, TString saveAs = "") // saveAs: without extension
{
   if(!Configure(type)) return;
   if (!validation) tree=(TTree*) f->Get("fittree");

   RooRealVar mass("mass", "mass [GeV/c^{2}]", mass_low, mass_high);
   RooRealVar t("t", "lifetime [ps]", -1.0*timeFactor, 15*timeFactor);
   RooRealVar tE("tE", "lifetime error", 1e-2*timeFactor, 1*timeFactor);
   if (!validation) RooDataSet data("data", "data", tree, RooArgSet(mass, t, tE));
 
   /*********************************************************************************************
    *    Prompt background = 1st order polynomial (mass) x double Gauss (lifetime)
    *********************************************************************************************/

   // prompt mass polynomial
   RooRealVar prompt_p1("prompt_p1","Linear coefficient of prompt background mass polynomial",0.1, -0.2, 5.0);
   RooPolynomial prompt_mass("prompt_mass", "Linear function for prompt background mass", mass, RooArgList(prompt_p1));
   // prompt lifetime resolution
   if(PerEventError){
         RooRealVar reso_bias("reso_bias", "Gauss mean of resolution function", 0,-10,10);
         RooRealVar reso_sigma("reso_sigma","Gauss sigma core of prompt background", 1.0, 0.1, 10);
         RooGaussModel prompt_lt("prompt_core","Gauss core for prompt background", t, reso_bias, reso_sigma,tE);
   } else {
         RooRealVar reso_mean("reso_mean", "Gauss mean of resolution function", 0);
         RooRealVar prompt_sigma_core("prompt_sigma_core","Gauss sigma core of prompt background", 0.1*timeFactor, 0.05*timeFactor, 0.12*timeFactor);
         RooRealVar prompt_sigma_tail("prompt_sigma_tail","Gauss sigma tail of prompt background", 0.2*timeFactor, 0.12*timeFactor, 0.9*timeFactor);
         RooGaussian prompt_core("prompt_core","Gauss core for prompt background", t, reso_mean, prompt_sigma_core);
         RooGaussian prompt_tail("prompt_tail","Gauss tail for prompt background", t, reso_mean, prompt_sigma_tail);
         RooRealVar core_frac("core_frac","Fraction of core in prompt background",0.8, 0.1, 0.99);
         RooAddModel prompt_lt("prompt_lt","Double Gauss for prompt background lifetime",RooArgList(prompt_core, prompt_tail), core_frac);
   }
   RooProdPdf prompt("prompt","Prompt 2D background model", RooArgSet(prompt_mass, prompt_lt));   


   /*********************************************************************************************
    *    Non-prompt background = 1st order polynomial (mass) x (Decay (x) Gauss)(lifetime)
    *********************************************************************************************/

   // non-prompt mass polynomial
   RooPolynomial nonprompt_mass("nonprompt_mass", "Linear function for nonprompt background mass", mass, RooArgList(prompt_p1));
      
   // non-prompt lifetime model
   RooRealVar tau_bk("tau_bk","Background lifetime",1.2*timeFactor, .005*timeFactor, 2*timeFactor);
   if(PerEventError) {
      RooRealVar npr_reso_bias("npr_reso_bias","Gauss mean of signal resolution", 0,-10,10); 
      RooRealVar npr_reso_sigma("npr_reso_sigma","Gauss sigma of signal resolution", 1.0, 0.05, 5);
      if(one_scale)
         RooGaussModel npr_reso("npr_reso", "Gauss of resolution function", t, reso_bias, reso_sigma,tE);
      else
         RooGaussModel npr_reso("npr_reso", "Gauss of resolution function", t, npr_reso_bias, npr_reso_sigma,tE);
   } else {
      RooRealVar npr_reso_mean("npr_reso_mean","Gauss mean of signal resolution", 0); 
      RooRealVar npr_reso_sigma("npr_reso_sigma","Gauss sigma of signal resolution", 0.1*timeFactor, 0.001*timeFactor, 10*timeFactor);
      RooGaussModel npr_reso("npr_reso", "Gauss of resolution function", t, npr_reso_mean, npr_reso_sigma);
   }
   RooDecay nonprompt_lt("nonprompt_lt","decay (x) double Gauss for background lifetime", t, tau_bk, npr_reso, RooDecay::SingleSided) ;

   RooProdPdf nonprompt("nonprompt", "Non-prompt background PDF", RooArgSet(nonprompt_mass, nonprompt_lt));   
   
   /***************************************************************************************
    *    Signal = Double Gauss (mass) x (Decay (x) double Gauss)(lifetime)
    ***************************************************************************************/

   //  mass peak 
   RooRealVar mass_peak("mass_peak","Gauss mean of signal mass peak", masspeak, masspeak-0.04, masspeak+0.04);
   RooRealVar m_sigma1("m_sigma1","Gauss core sigma for signal mass",0.01,0.001,0.015); 
   RooRealVar m_sigma2("m_sigma2","Gauss tail sigma for signal mass",0.025, 0.02, 0.03);
   //RooRealVar m_sigma1("m_sigma1","Gauss core sigma for signal mass",0.01,0.001,0.015); 
   //RooRealVar m_sigma2("m_sigma2","Gauss tail sigma for signal mass",0.012, 0.008, 0.025);
   RooGaussian m_gauss1("m_gauss1","Core Gauss for signal mass", mass, mass_peak, m_sigma1);
   RooGaussian m_gauss2("m_gauss2","Tail Gauss for signal mass", mass, mass_peak, m_sigma2);
   RooRealVar frac_m_gauss("frac_m_gauss","Fraction of tail Gauss for signal mass",0.2, 0.1, 0.99);
   RooAddPdf mgauss("mgauss","Double gauss for signal mass",RooArgList(m_gauss1, m_gauss2),RooArgList(frac_m_gauss));

   // signal lifetime
   RooRealVar tau("tau","Lambda_b lifetime",1.5*timeFactor, 0.1*timeFactor, 2*timeFactor);
   if(PerEventError) {
      RooRealVar sig_reso_bias("sig_reso_bias","Gauss sigma of signal resolution", 0,-10,10); 
      RooRealVar sig_reso_sigma("sig_reso_sigma","Gauss sigma of signal resolution", 1.0, 0.05, 5.);
      if(one_scale)
         RooGaussModel sig_reso("sig_reso", "Gauss of resolution function", t, reso_bias, reso_sigma, tE);
      else
         RooGaussModel sig_reso("sig_reso", "Gauss of resolution function", t, sig_reso_bias, sig_reso_sigma, tE);
   } else {
      RooRealVar sig_reso_mean("sig_reso_mean","Gauss mean of signal resolution", 0); 
      RooRealVar sig_reso_sigma("sig_reso_sigma","Gauss sigma of signal resolution", 0.1*timeFactor, 0.001*timeFactor, 10.*timeFactor);
      RooGaussModel sig_reso("sig_reso", "Gauss of resolution function", t, sig_reso_mean, sig_reso_sigma);
   }
   RooDecay tau_decay_model("tau_decay_model","decay (x) double Gauss for signal lifetime", t, tau, sig_reso, RooDecay::SingleSided) ;

   RooProdPdf signal("signal", "Signal PDF", RooArgSet(mgauss, tau_decay_model));


   /*********************************************************************************************
    *    Complete 2D model
    *********************************************************************************************/
 
   RooRealVar n_prompt("n_prompt","Number of prompt background events",1000, 0, 50000);
   RooRealVar n_nonprompt("n_nonprompt","Number of non-prompt background events",300, 0, 20000);
   RooRealVar n_signal("n_signal","Number of signal events",4000, 10, 20000);
   RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal, prompt, nonprompt), RooArgList(n_signal, n_prompt, n_nonprompt));

   //   RooFormulaVar eff("eff","1.0 - t*1.09e9",t) ;
   RooFormulaVar eff("eff","0.001574*(1.0 - 0.04739*t+1.038/(1+TMath::Exp(-t/0.1865)))",t) ;
   //   RooEffProd full_model("full_model","model with efficiency", model, eff) ;

   // extract number of parameters in use
   RooArgSet dummy;
   const int Npars = full_model.getParameters(dummy)->getSize();
   cout << "Npars: " << Npars << endl;


   t.setRange("all", -1.0*timeFactor, 15*timeFactor);

   /*********************************************************************************************
    *    Perform fit validation and stop execution
    *********************************************************************************************/ 

   if (validation)
   {
       RooMCStudy* mcstudy = new RooMCStudy(full_model, RooArgSet(mass,t), Binned(kTRUE), Silence(), Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(0), NumCPU(6)));
       mcstudy->generateAndFit(2);

       // E x p l o r e   r e s u l t s   o f   s t u d y 
       // ------------------------------------------------

       // Make plots of the distributions of tau, the error on tau and the pull of tau
       vector<RooPlot*> plotvec;
       plotvec.push_back(mcstudy->plotParam(tau,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(tau,Bins(40)));
       plotvec.push_back(mcstudy->plotParam(tau,Bins(40),FitGauss(kTRUE)));

       plotvec.push_back(mcstudy->plotError(mass_peak,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(mass_peak,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotParam(mass_peak,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(m_sigma1,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(m_sigma1,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotParam(m_sigma1,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(m_sigma2,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(m_sigma2,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotPull(m_sigma2,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(frac_m_gauss,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(frac_m_gauss,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotPull(frac_m_gauss,Bins(40),FitGauss(kTRUE)));

       plotvec.push_back(mcstudy->plotParam(n_prompt,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(n_prompt,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(n_prompt,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotParam(n_nonprompt,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(n_nonprompt,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(n_nonprompt,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotParam(n_signal,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(n_signal,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(n_signal,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotParam(prompt_p1,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(prompt_p1,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(prompt_p1,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotParam(tau_bk,Bins(40),FitGauss(kTRUE)));
       plotvec.push_back(mcstudy->plotError(tau_bk,Bins(40)));
       plotvec.push_back(mcstudy->plotPull(tau_bk,Bins(40),FitGauss(kTRUE)));
       if(PerEventError)
       {
	   plotvec.push_back(mcstudy->plotParam(reso_bias,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(reso_bias,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(reso_bias,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotParam(reso_sigma,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(reso_sigma,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(reso_sigma,Bins(40),FitGauss(kTRUE)));
       }
       if(PerEventError && !one_scale)
       {
	   plotvec.push_back(mcstudy->plotParam(npr_reso_bias,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(npr_reso_bias,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(npr_reso_bias,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotParam(npr_reso_sigma,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(npr_reso_sigma,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(npr_reso_sigma,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotParam(sig_reso_bias,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(sig_reso_bias,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(sig_reso_bias,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotParam(sig_reso_sigma,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(sig_reso_sigma,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(sig_reso_sigma,Bins(40),FitGauss(kTRUE)));
       }
       if(!PerEventError)
       {
	   plotvec.push_back(mcstudy->plotParam(prompt_sigma_core,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(prompt_sigma_core,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(prompt_sigma_core,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotParam(prompt_sigma_tail,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(prompt_sigma_tail,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(prompt_sigma_tail,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotParam(core_frac,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(core_frac,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(core_frac,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotParam(npr_reso_sigma,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(npr_reso_sigma,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(npr_reso_sigma,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotParam(sig_reso_sigma,Bins(40),FitGauss(kTRUE)));
	   plotvec.push_back(mcstudy->plotError(sig_reso_sigma,Bins(40)));
	   plotvec.push_back(mcstudy->plotPull(sig_reso_sigma,Bins(40),FitGauss(kTRUE)));
       }
       
       
       // Plot distribution of minimized likelihood
       plotvec.push_back(mcstudy->plotNLL(Bins(40)));
       
      // Draw all plots on a canvas
      gStyle->SetPalette(1);
      gStyle->SetOptStat(0);
      TCanvas* c = new TCanvas("rf801_mcstudy","rf801_mcstudy",600,600);
      int i=0;
      for (vector<RooPlot*>::const_iterator it=plotvec.begin(); it!=plotvec.end(); it++, i++)
      {
	  gPad->SetLeftMargin(0.15);
	  (*it)->GetYaxis()->SetTitleOffset(1.4);
	  (*it)->Draw();
	  gPad->SaveAs(("doFit2DNPlots_validation_"+toString(i)+".pdf").c_str());
      }

       return;
       // Make some histograms from the parameter dataset
       TH1* hh_cor_a0_s1f = mcstudy->fitParDataSet().createHistogram("hh",a1,YVar(sig1frac)) ;
       TH1* hh_cor_a0_a1  = mcstudy->fitParDataSet().createHistogram("hh",a0,YVar(a1)) ;
       
       // Access some of the saved fit results from individual toys
       TH2* corrHist000 = mcstudy->fitResult(0)->correlationHist("c000") ;
       TH2* corrHist127 = mcstudy->fitResult(127)->correlationHist("c127") ;
       TH2* corrHist953 = mcstudy->fitResult(953)->correlationHist("c953") ;

   }

   /*********************************************************************************************
    *    Do the fitting
    *********************************************************************************************/ 

   if (PerEventError)
     full_model.fitTo(data, ConditionalObservables(tE), Save(kTRUE));
   else
     full_model.fitTo(data, Save(kTRUE));

   RooAbsReal *cdf=mgauss.createCdf(mass);
   two_sigma_upper=cdf->findRoot(mass, mass_low, mass_high, 0.97725);
   two_sigma_lower=cdf->findRoot(mass, mass_low, mass_high, 0.02275);
   three_sigma_upper=cdf->findRoot(mass, mass_low, mass_high, 0.99865);
   three_sigma_lower=cdf->findRoot(mass, mass_low, mass_high, 1.0-0.99865);

   /*********************************************************************************************
    *    Report results in one line for automated result collection
    *********************************************************************************************/ 

   cout << "FITRESULTS: ";
   cout << tau.getVal()*invTimeFactor << " " << tau.getError()*invTimeFactor << " "; // ok
   cout << mass_peak.getVal() << " " << mass_peak.getError() << " "; // ok
   cout << m_sigma1.getVal() << " " << m_sigma1.getError() << " "; // ok
   cout << m_sigma2.getVal() << " " << m_sigma2.getError() << " "; // ok
   cout << frac_m_gauss.getVal() << " " << frac_m_gauss.getError() << " "; // ok
   if(!PerEventError) cout << core_frac.getVal() << " " << core_frac.getError() << " "; // ok
   cout << tau_bk.getVal()*invTimeFactor << " " << tau_bk.getError()*invTimeFactor << " "; // ok
   cout << n_nonprompt.getVal() << " " << n_nonprompt.getError() << " "; // ok
   cout << n_prompt.getVal() << " " << n_prompt.getError() << " "; // ok
   cout << n_signal.getVal() << " " << n_signal.getError() << " "; // ok
   cout << prompt_p1.getVal() << " " << prompt_p1.getError() << " "; // ok
   if(!PerEventError) cout << prompt_sigma_core.getVal() << " " << prompt_sigma_core.getError() << " "; // ok
   if(!PerEventError) cout << prompt_sigma_tail.getVal() << " " << prompt_sigma_tail.getError() << " "; // ok
   if(PerEventError) cout << reso_bias.getVal() << " " << reso_bias.getError() << " "; // ok
   if(PerEventError) cout << reso_sigma.getVal() << " " << reso_sigma.getError() << " "; // ok
   if(PerEventError && !one_scale) cout << npr_reso_bias.getVal() << " " << npr_reso_bias.getError() << " "; // ok
   if(PerEventError && !one_scale) cout << npr_reso_sigma.getVal() << " " << npr_reso_sigma.getError() << " "; // ok
   if(PerEventError && !one_scale) cout << sig_reso_bias.getVal() << " " << sig_reso_bias.getError() << " "; // ok
   if(PerEventError && !one_scale) cout << sig_reso_sigma.getVal() << " " << sig_reso_sigma.getError() << " "; // ok
   if(!PerEventError) cout << npr_reso_sigma.getVal() << " " << npr_reso_sigma.getError() << " "; // ok
   if(!PerEventError) cout << sig_reso_sigma.getVal() << " " << sig_reso_sigma.getError() << " "; // ok
   //cout << chi2_m << " " << chi2_t << " ";
   cout << "Settings: ";
   if(PerEventError) cout << "PEE ";
   if(one_scale) cout << "1sc "; else cout << "3sc ";
   cout << endl;
   if (noplot) return;

   /*********************************************************************************************
    *     Plot results
    *********************************************************************************************/  

   TH2* hd = data.createHistogram("hd",t,Binning(60),YVar(mass,Binning(20)));
   TH2* hf = full_model.createHistogram("hf",t,Binning(60),YVar(mass,Binning(20))) ;

   c1=new TCanvas(name_c1, title1);
   c1->Divide(2) ;
   c1->cd(1); 
   gPad->SetLogz();
   hd->Draw("surf") ;
   c1->cd(2); 
   gPad->SetLogz();
   hf->Draw("surf") ;
   if (saveAs.Length() !=0 ) c1->SaveAs(saveAs+"_2d.pdf");

   c2=new TCanvas(name_c2, title2, 1100, 850);
   c2->Divide(2) ;

   c2->cd(1);
   TPad *m_pad1 = new TPad("m_pad1","m_pad1",0,0.3,1,1);
   m_pad1->SetBottomMargin(0.002);
   m_pad1->SetLeftMargin(0.12);
   m_pad1->SetRightMargin(0.004);
   m_pad1->Draw();
   m_pad1->cd();

   RooPlot* framex = mass.frame(Title("Mass projection")) ;
   data.plotOn(framex, Name("data_m")) ;
   if(PerEventError) {
      full_model.plotOn(framex, ProjWData(data), Range("all"), LineColor(kBlue), Name("model_m"));
      full_model.plotOn(framex, ProjWData(data), Components(signal), LineColor(kRed));
      full_model.plotOn(framex, ProjWData(data), Components(prompt), LineStyle(kDashed), LineColor(kBlack));
      full_model.plotOn(framex, ProjWData(data), Components(nonprompt), LineStyle(kDotted), LineColor(kBlack));
   } else {
      full_model.plotOn(framex, Range("all"), LineColor(kBlue), Name("model_m"));
      full_model.plotOn(framex, Components(signal), LineColor(kRed));
      full_model.plotOn(framex, Components(prompt), LineStyle(kDashed), LineColor(kBlack));
      full_model.plotOn(framex, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack));
   }
   // calculate chi2
   const double chi2_m = framex->chiSquare("model_m","data_m", Npars) ;
   framex->addObject(writeTLatex("mass: "+roundToString(mass_peak.getVal(),3)+" #pm "+roundToString(mass_peak.getError()+0.0004, 3)+" GeV/c^{2}", 0.32, 0.85, 0.05));
   framex->addObject(writeTLatex("#chi^{2}/ndof: "+roundToString(chi2_m, 3), 0.5, 0.80, 0.04));
   framex->GetYaxis()->SetTitleOffset(1.6);
   framex->SetMinimum(0.01);
   framex->Draw();

   c2->cd(1);
   TPad *m_pad2 = new TPad("m_pad2","m_pad2",0,0,1,0.3);
   m_pad2->SetTopMargin(0.002);
   m_pad2->SetBottomMargin(0.2);
   m_pad2->SetLeftMargin(0.12);
   m_pad2->SetRightMargin(0.004);
   m_pad2->Draw();
   m_pad2->cd();

   RooHist *ratioplot_mass = framex->residHist("data_m", "model_m", true);
   RooPlot* framex_ratio = mass.frame(Title("  "));
   framex_ratio->addPlotable(ratioplot_mass, "P");
   framex_ratio->GetYaxis()->SetLabelSize(0.07);
   framex_ratio->GetYaxis()->SetTitleSize(0.07);
   framex_ratio->GetYaxis()->SetTitleOffset(0.70);
   framex_ratio->GetYaxis()->SetTitle("pull");
   framex_ratio->GetXaxis()->SetLabelSize(0.07);
   framex_ratio->GetXaxis()->SetTitleSize(0.07);
   framex_ratio->Draw();


   c2->cd(2);
   TPad *t_pad1 = new TPad("t_pad1","t_pad1",0,0.3,1,1);
   t_pad1->SetBottomMargin(0.002);
   t_pad1->SetLeftMargin(0.12);
   t_pad1->SetRightMargin(0.004);
   t_pad1->Draw();
   t_pad1->cd();

   gPad->SetLogy(); 
   RooPlot* framey = t.frame(Title("Lifetime projection")) ;
   data.plotOn(framey, Name("data_t"), Binning(100)) ;
   if(PerEventError) {
      full_model.plotOn(framey, ProjWData(data), Range("all"), LineColor(kBlue), Name("model_t")) ;
      full_model.plotOn(framey, ProjWData(data), Components(signal), LineColor(kRed));
      full_model.plotOn(framey, ProjWData(data), Components(prompt), LineStyle(kDashed), LineColor(kBlack));
      full_model.plotOn(framey, ProjWData(data), Components(nonprompt), LineStyle(kDotted), LineColor(kBlack));
   } else {
      full_model.plotOn(framey, Range("all"), LineColor(kBlue), Name("model_t")) ;
      full_model.plotOn(framey, Components(signal), LineColor(kRed));
      full_model.plotOn(framey, Components(prompt), LineStyle(kDashed), LineColor(kBlack));
      full_model.plotOn(framey, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack));
   }
   // calculate chi2
   const double chi2_t = framey->chiSquare("model_t","data_t", Npars) ;

   framey->addObject(writeTLatex("#tau: "+roundToString(tau.getVal()*invTimeFactor, 3)+" #pm "+roundToString(tau.getError()*invTimeFactor, 3)+" ps", 0.5, 0.85, 0.05));
   framey->addObject(writeTLatex("#chi^{2}/ndof: "+roundToString(chi2_t, 3), 0.5, 0.80, 0.04));
   gPad->SetLeftMargin(0.10);
   gPad->SetRightMargin(0.003);
   framey->GetYaxis()->SetTitleOffset(1.5);
   framey->Draw();

   c2->cd(2);
   TPad *t_pad2 = new TPad("t_pad2","t_pad2",0,0,1,0.3);
   t_pad2->SetTopMargin(0.002);
   t_pad2->SetBottomMargin(0.2);
   t_pad2->SetLeftMargin(0.10);
   t_pad2->SetRightMargin(0.004);
   t_pad2->Draw();
   t_pad2->cd();

   RooHist *ratioplot_t = framey->residHist("data_t", "model_t", true);
   RooPlot* framey_ratio = t.frame(Title("  "));
   framey_ratio->addPlotable(ratioplot_t , "P");
   framey_ratio->GetYaxis()->SetLabelSize(0.07);
   framey_ratio->GetYaxis()->SetTitleSize(0.07);
   framey_ratio->GetYaxis()->SetTitleOffset(0.70);
   framey_ratio->GetYaxis()->SetTitle("pull");
   framey_ratio->GetXaxis()->SetLabelSize(0.07);
   framey_ratio->GetXaxis()->SetTitleSize(0.07);
   framey_ratio->Draw();

   if (saveAs.Length() !=0 ) c2->SaveAs(saveAs+"_m_t.pdf");

   //--------------- Sideband projection plots
   c3=new TCanvas(name_c3, title3, 1100, 550);
   c3->cd();
   c3->Divide(3);
   c3->cd(1);

   TPad *sbl_pad1 = new TPad("sbl_pad1","sbl_pad1",0,0.3,1,1);
   sbl_pad1->SetBottomMargin(0.002);
   sbl_pad1->SetLeftMargin(0.12);
   sbl_pad1->SetRightMargin(0.004);
   sbl_pad1->Draw();
   sbl_pad1->cd();

   gPad->SetLogy();
   mass.setRange("sbl", mass_low, three_sigma_lower);
   RooPlot* framesbl = t.frame(Title("Lower sideband projection"), Range("sbl")) ;
   data.plotOn(framesbl, CutRange("sbl")) ;
   if(PerEventError) {
      full_model.plotOn(framesbl, ProjWData(tE, data), ProjectionRange("sbl")) ;
      full_model.plotOn(framesbl, ProjWData(tE, data), Components(signal), LineColor(kRed), ProjectionRange("sbl"));
      full_model.plotOn(framesbl, ProjWData(tE, data), Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sbl"));
      full_model.plotOn(framesbl, ProjWData(tE, data), Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sbl"));
   } else {
      full_model.plotOn(framesbl, ProjectionRange("sbl")) ;
      full_model.plotOn(framesbl, Components(signal), LineColor(kRed), ProjectionRange("sbl"));
      full_model.plotOn(framesbl, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sbl"));
      full_model.plotOn(framesbl, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sbl"));
   }
   gPad->SetLeftMargin(0.12);
   gPad->SetRightMargin(0.004);
   framesbl->GetYaxis()->SetTitleOffset(1.2);
   framesbl->GetYaxis()->SetTitleSize(0.045);
   framesbl->GetXaxis()->SetTitleSize(0.045);
   framesbl->Draw();

   c3->cd(1);
   TPad *sbl_pad2 = new TPad("sbl_pad2","sbl_pad2",0,0,1,0.3);
   sbl_pad2->SetTopMargin(0.002);
   sbl_pad2->SetBottomMargin(0.2);
   sbl_pad2->SetLeftMargin(0.12);
   sbl_pad2->SetRightMargin(0.004);
   sbl_pad2->Draw();
   sbl_pad2->cd();

   if(PerEventError)
      RooHist *ratioplot_t = framesbl->residHist("h_data_CutRange[sbl]", "full_model_Int[mass|sbl]_Norm[mass,t]_DataAvg[tE]",true);
   else
      RooHist *ratioplot_t = framesbl->residHist("h_data_CutRange[sbl]", "full_model_Int[mass|sbl]_Norm[mass,t]", true);
   RooPlot* framesbl_ratio = t.frame(Title("  "));
   framesbl_ratio->addPlotable(ratioplot_t , "P");
   framesbl_ratio->GetYaxis()->SetLabelSize(0.10);
   framesbl_ratio->GetYaxis()->SetTitleSize(0.10);
   framesbl_ratio->GetYaxis()->SetTitleOffset(0.45);
   framesbl_ratio->GetYaxis()->SetTitle("pull");
   framesbl_ratio->GetXaxis()->SetLabelSize(0.10);
   framesbl_ratio->GetXaxis()->SetTitleSize(0.10);
   framesbl_ratio->Draw();

   c3->cd(2);
   TPad *sig_pad1 = new TPad("sig_pad1","sig_pad1",0,0.3,1,1);
   sig_pad1->SetBottomMargin(0.002);
   sig_pad1->SetLeftMargin(0.12);
   sig_pad1->SetRightMargin(0.004);
   sig_pad1->Draw();
   sig_pad1->cd();

   gPad->SetLogy();
   mass.setRange("sig", two_sigma_lower, two_sigma_upper);
   RooPlot* framesig = t.frame(Title("Signal region projection"),Range("sig")) ;
   data.plotOn(framesig, CutRange("sig")) ;
   if(PerEventError) {
      full_model.plotOn(framesig, ProjWData(tE, data), ProjectionRange("sig")) ;
      full_model.plotOn(framesig, ProjWData(tE, data), Components(signal), LineColor(kRed), ProjectionRange("sig"));
      full_model.plotOn(framesig, ProjWData(tE, data), Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sig"));
      full_model.plotOn(framesig, ProjWData(tE, data), Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sig"));
   } else {
      full_model.plotOn(framesig, ProjectionRange("sig")) ;
      full_model.plotOn(framesig, Components(signal), LineColor(kRed), ProjectionRange("sig"));
      full_model.plotOn(framesig, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sig"));
      full_model.plotOn(framesig, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sig"));
   }
   framesig->GetYaxis()->SetTitleOffset(1.2);
   framesig->GetYaxis()->SetTitleSize(0.045);
   framesig->GetXaxis()->SetTitleSize(0.045);
   framesig->Draw();

   c3->cd(2);
   TPad *sig_pad2 = new TPad("sig_pad2","sig_pad2",0,0,1,0.3);
   sig_pad2->SetTopMargin(0.002);
   sig_pad2->SetBottomMargin(0.2);
   sig_pad2->SetLeftMargin(0.12);
   sig_pad2->SetRightMargin(0.004);
   sig_pad2->Draw();
   sig_pad2->cd();

   if(PerEventError)
      RooHist *ratioplot_t = framesig->residHist("h_data_CutRange[sig]", "full_model_Int[mass|sig]_Norm[mass,t]_DataAvg[tE]", true);
   else
      RooHist *ratioplot_t = framesig->residHist("h_data_CutRange[sig]", "full_model_Int[mass|sig]_Norm[mass,t]", true);
   RooPlot* framesig_ratio = t.frame(Title("  "));
   framesig_ratio->addPlotable(ratioplot_t , "P");
   framesig_ratio->GetYaxis()->SetLabelSize(0.10);
   framesig_ratio->GetYaxis()->SetTitleSize(0.10);
   framesig_ratio->GetYaxis()->SetTitleOffset(0.45);
   framesig_ratio->GetYaxis()->SetTitle("pull");
   framesig_ratio->GetXaxis()->SetLabelSize(0.10);
   framesig_ratio->GetXaxis()->SetTitleSize(0.10);
   framesig_ratio->Draw();


   c3->cd(3);
   TPad *sbh_pad1 = new TPad("sbh_pad1","sbh_pad1",0,0.3,1,1);
   sbh_pad1->SetBottomMargin(0.002);
   sbh_pad1->SetLeftMargin(0.12);
   sbh_pad1->SetRightMargin(0.004);
   sbh_pad1->Draw();
   sbh_pad1->cd();

   gPad->SetLogy();
   mass.setRange("sb", three_sigma_upper, mass_high);
   RooPlot* framesbh = t.frame(Title("Upper sideband projection"),Range("sb")) ;
   data.plotOn(framesbh, CutRange("sb")) ;
   if(PerEventError) {
      full_model.plotOn(framesbh, ProjWData(tE, data), ProjectionRange("sb")) ;
      full_model.plotOn(framesbh, ProjWData(tE, data), Components(signal), LineColor(kRed), ProjectionRange("sb"));
      full_model.plotOn(framesbh, ProjWData(tE, data), Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sb"));
      full_model.plotOn(framesbh, ProjWData(tE, data), Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sb"));
   } else {
      full_model.plotOn(framesbh, ProjectionRange("sb")) ;
      full_model.plotOn(framesbh, Components(signal), LineColor(kRed), ProjectionRange("sb"));
      full_model.plotOn(framesbh, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sb"));
      full_model.plotOn(framesbh, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sb"));
   }
   gPad->SetLeftMargin(0.12);
   gPad->SetRightMargin(0.004);
   framesbh->GetYaxis()->SetTitleOffset(1.2);
   framesbh->GetYaxis()->SetTitleSize(0.045);
   framesbh->GetXaxis()->SetTitleSize(0.045);
   framesbh->Draw();

   c3->cd(3);
   TPad *sbh_pad2 = new TPad("sbh_pad2","sbh_pad2",0,0,1,0.3);
   sbh_pad2->SetTopMargin(0.002);
   sbh_pad2->SetBottomMargin(0.2);
   sbh_pad2->SetLeftMargin(0.12);
   sbh_pad2->SetRightMargin(0.004);
   sbh_pad2->Draw();
   sbh_pad2->cd();

   if(PerEventError)
      RooHist *ratioplot_t = framesbh->residHist("h_data_CutRange[sb]", "full_model_Int[mass|sb]_Norm[mass,t]_DataAvg[tE]", true);
   else
      RooHist *ratioplot_t = framesbh->residHist("h_data_CutRange[sb]", "full_model_Int[mass|sb]_Norm[mass,t]", true);
   RooPlot* framesbh_ratio = t.frame(Title("  "));
   framesbh_ratio->addPlotable(ratioplot_t , "P");
   framesbh_ratio->GetYaxis()->SetLabelSize(0.10);
   framesbh_ratio->GetYaxis()->SetTitleSize(0.10);
   framesbh_ratio->GetYaxis()->SetTitleOffset(0.45);
   framesbh_ratio->GetYaxis()->SetTitle("pull");
   framesbh_ratio->GetXaxis()->SetLabelSize(0.10);
   framesbh_ratio->GetXaxis()->SetTitleSize(0.10);
   framesbh_ratio->Draw();

   if (saveAs.Length() !=0 ) c3->SaveAs(saveAs+"_sidebands.pdf");

   cout << "Fitresults: t: "<< tau.getVal()*invTimeFactor << " +- "<<tau.getError()*invTimeFactor<<" ps ";
   cout << " chisq(m): " << framex->chiSquare() << " chisq(t): " << framey->chiSquare();
   cout << " mass: " << mass_peak.getVal() << " m_sigma1: " << m_sigma1.getVal() << " m_sigma2: " << m_sigma2.getVal() << " tau_bk: " << tau_bk.getVal();
   cout <<  " nsig: " << n_signal.getVal()<<" npr: "<<n_prompt.getVal()<<" nnpr: "<<n_nonprompt.getVal()<<endl;
   //f->Close();
   //f->Delete();
   //delete f;

}


int Configure(TString arguments)
{
  PerEventError=false;
  one_scale=true;
  noplot=false;
  validation=false;
  type = Lambda_b;
  if(arguments.Contains("pee1",TString::kIgnoreCase) ||
     arguments.Contains("pee",TString::kIgnoreCase)) {
    PerEventError = true;
    one_scale=true;
  }
  if(arguments.Contains("pee3",TString::kIgnoreCase)) {
    PerEventError = true;
    one_scale=false;
  }
  if(arguments.Contains("B0",TString::kIgnoreCase)) type=B0;
  if(arguments.Contains("noplot",TString::kIgnoreCase)) noplot=true;
  if(arguments.Contains("validation",TString::kIgnoreCase)) validation=true;

  TObjArray* tokens=arguments.Tokenize(" ");
  bool FileNameEntered=false;
  for(Int_t i=0; i<tokens->GetEntries(); i++){
    TObjString* str =dynamic_cast<TObjString*> (tokens->At(i));
    if(!str) continue;
    TString token=str->GetString();
    if(token.Contains(".root")) {
      FileNameEntered=true;
      f=new TFile(token);
    }
  }
  tokens->Clear();
  delete tokens;

  switch(type) {
  case B0:                    
     mass_low=5.16;
     mass_high=5.75;
     masspeak=5.28;
     if(!FileNameEntered && !validation) f=new TFile("vrt_r479_B0_data_B008.root");
     title1=TString("B^{0} mass, data, barrel trigger");  
     title2=TString("B^{0} lifetime sideband, data, barrel trigger");  
     title3=TString("B^{0} lifetime, data, barrel trigger"); 
     break;
  case Lambda_b:
  default:
     mass_low=5.4;
     mass_high=6.0;
     masspeak=5.62;
     if(!FileNameEntered && !validation) f=new TFile("vrt_r548_lb_data_lb14_ps.root");
     title1=TString("#Lambda_{b} mass, data, barrel trigger");  
     title2=TString("#Lambda_{b} lifetime sideband, data, barrel trigger");  
     title3=TString("#Lambda_{b} lifetime, data, barrel trigger"); 
  }
  if(!validation && f->IsZombie()) return 0;
  else return 1;
}


void doFit2DNPlots(TString filestem, int N, TString arguments)
{
    for(int i=0; i!=N; i++)
    {
	cout << "Working on file no. " << i << endl;
	Fit_2D(filestem+TString(toString(i))+".root " + arguments);
    }

    return;
    
    //Fit_2D("../data/vrt_r548_lb_data_lb14.root");
}

