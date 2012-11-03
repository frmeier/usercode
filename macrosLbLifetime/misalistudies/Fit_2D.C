#include "Configure.C"
#include <iomanip>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "utils.h"

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

using namespace RooFit;
using namespace std;

TTree *tree;
ofstream result;
ConfigData cfg;

Double_t yield_nonprompt, yield_sig, yield_prompt;
Double_t two_sigma_upper, two_sigma_lower, three_sigma_upper, three_sigma_lower;
TCanvas *c1, *c2, *c3;
TString name_c1("M_vs_lt"), name_c2("projections"), name_c3("slices");

void Fit_2D(TString type)
{
   if(!Configure(type, cfg)) return;
   tree=(TTree*) cfg.f->Get("fittree");

   RooRealVar mass("mass", "mass [GeV/c^{2}]", cfg.mass_low, cfg.mass_high);
   RooRealVar t("t", "lifetime [ps]", -1.0e-12, 15e-12);
   RooRealVar tE("tE", "lifetime error", 0, 1e-12);
   RooRealVar weight("weight","Event weight", 1, 10);
   RooDataSet data("data", "data", tree, RooArgSet(mass, t, tE));
   //   RooDataSet data("data", "data", tree, RooArgSet(mass, t, weight), "", "weight");
 
   /*********************************************************************************************
    *    Non-prompt background = 1st order polynomial (mass) x (Decay (x) double Gauss)(lifetime)
    *    Prompt background = 1st order polynomial (mass) x double Gauss (lifetime)
    *********************************************************************************************/

   if(cfg.prompt_bk){
      RooRealVar reso_mean("reso_mean", "mean of resolution function", 0);
   // prompt mass polynomial
      RooRealVar prompt_p1("prompt_p1","Linear coefficient of prompt background mass polynomial",0.1, -0.2, 5.0);
      RooPolynomial prompt_mass("prompt_mass", "Linear function for prompt background mass", mass, RooArgList(prompt_p1));
   // prompt lifetime resolution
      RooRealVar prompt_sigma_core("prompt_sigma_core","Gauss sigma core of prompt background", 1e-13, 5e-14, 1.2e-13);
      RooRealVar prompt_sigma_tail("prompt_sigma_tail","Gauss sigma tail of prompt background", 2e-13, 1.2e-13, 9e-13);
      RooGaussian prompt_core("prompt_core","Gauss core for prompt background", t, reso_mean, prompt_sigma_core);
      RooGaussian prompt_tail("prompt_tail","Gauss tail for prompt background", t, reso_mean, prompt_sigma_tail);
      RooRealVar core_frac("core_frac","Fraction of core in prompt background",0.8, 0.5, 0.99);
      RooAddPdf prompt_lt("prompt_lt","Double Gauss for prompt background lifetime",RooArgList(prompt_core, prompt_tail), core_frac);
      RooProdPdf prompt("prompt","Prompt 2D background model", RooArgSet(prompt_mass, prompt_lt));   

      RooGaussModel reso_core("reso_core", "First Gauss of resolution function", t, reso_mean, prompt_sigma_core);
      RooGaussModel reso_tail("reso_tail", "Second Gauss of resolution function", t, reso_mean, prompt_sigma_tail);
      RooAddModel tau_reso("tau_reso", "Double Gauss of resolution function", RooArgList(reso_core, reso_tail), core_frac);
   } else {
      RooRealVar tau_reso_sigma("tau_reso_sigma","Sigma of lifetime resolution function",0.1e-13, 1e-15,9.0e-13) ;
      //      tau_reso_sigma=6.26e-14;
      //      tau_reso_sigma.setConstant(kTRUE);
      RooGaussModel tau_reso("tau_reso","Gauss of lifetime resolution function", t, RooConst(0), tau_reso_sigma) ;
   }

   if(cfg.nonprompt_bk){
      // non-prompt mass polynomial
      if(cfg.constant_prompt_fraction){
	 cout<<"Constant prompt-to-nonprompt fraction in mass."<<endl;
         RooPolynomial nonprompt_mass("nonprompt_mass", "Linear function for nonprompt background mass", mass, RooArgList(prompt_p1));
         RooPolynomial nonprompt_mass2("nonprompt_mass2", "Linear function for nonprompt background mass", mass, RooArgList(prompt_p1));
      } else {
	 cout<<"Different mass polynomials for prompt and nonprompt backgrounds."<<endl;
         RooRealVar nonprompt_p1("nonprompt_p1","Linear coefficient of nonprompt background mass polynomial", 0.1, -0.2, 0.3);
         RooRealVar nonprompt_p2("nonprompt_p2","Linear coefficient of nonprompt background mass polynomial", 0.1, -0.2, 0.3);
         RooPolynomial nonprompt_mass("nonprompt_mass", "Linear function for nonprompt background mass", mass, RooArgList(nonprompt_p1));
         RooPolynomial nonprompt_mass2("nonprompt_mass", "Linear function for nonprompt background mass", mass, RooArgList(nonprompt_p2));
      }
      // non-prompt lifetime model
      RooRealVar tau_bk("tau_bk","Background lifetime",1.2e-12, 5e-15, 2e-12);
      RooDecay nonprompt_lt("nonprompt_lt","decay (x) double Gauss for background lifetime", t, tau_bk, tau_reso, RooDecay::SingleSided) ;
      if(cfg.single_lt){
	 cout<<"Nonprompt background has one lifetime component."<<endl;
         RooProdPdf nonprompt("nonprompt", "Non-prompt background PDF", RooArgSet(nonprompt_mass, nonprompt_lt));   
      } else {
	 cout<<"Nonprompt background has two lifetime components."<<endl;
         RooRealVar tau_bk2("tau_bk2","Second background lifetime",0.6e-12, 5e-15, 2e-12);
         RooDecay nonprompt_lt2("nonprompt_lt2","decay (x) double Gauss for second background lifetime", t, tau_bk2, tau_reso, RooDecay::SingleSided) ;
         RooProdPdf nonprompt1("nonprompt1", "Non-prompt background PDF", RooArgSet(nonprompt_mass, nonprompt_lt));   
         RooProdPdf nonprompt2("nonprompt2", "Non-prompt background PDF", RooArgSet(nonprompt_mass, nonprompt_lt2));
         RooRealVar lt_frac("lt_frac", "fraction of longer background lifetime", 0.2, 0.01, 0.99);
         RooAddPdf nonprompt("nonprompt", "Non-prompt background PDF", RooArgList(nonprompt1, nonprompt2), RooArgList(lt_frac));
      }
   }
   
   /***************************************************************************************
    *    Signal = Double Gauss (mass) x (Decay (x) double Gauss)(lifetime)
    ***************************************************************************************/

   //  mass peak 
   RooRealVar mass_peak("mass_peak","Gauss mean of signal mass peak", cfg.masspeak, cfg.masspeak-0.04, cfg.masspeak+0.04);
   RooRealVar m_sigma1("m_sigma1","Gauss core sigma for signal mass",0.007,0.001,0.012);
   RooRealVar m_sigma2("m_sigma2","Gauss tail sigma for signal mass",0.02, 0.010, 0.035);
   RooGaussian m_gauss1("m_gauss1","Core Gauss for signal mass", mass, mass_peak, m_sigma1);
   RooGaussian m_gauss2("m_gauss2","Tail Gauss for signal mass", mass, mass_peak, m_sigma2);
   RooRealVar frac_m_gauss("frac_m_gauss","Fraction of tail Gauss for signal mass",0.2, 0.1, 0.99);
   if (cfg.single_sig)
       RooGaussian mgauss("mgauss","Single gauss for signal mass", mass, mass_peak, m_sigma1);
   else
       RooAddPdf mgauss("mgauss","Double gauss for signal mass",RooArgList(m_gauss1, m_gauss2),RooArgList(frac_m_gauss));

   // signal lifetime
   RooRealVar tau("tau","Lambda_b lifetime",1.5e-12, 0.1e-12, 2e-12);
   RooDecay tau_decay_model("tau_decay_model","decay (x) double Gauss for signal lifetime", t, tau, tau_reso, RooDecay::SingleSided) ;
   RooProdPdf signal("signal", "Signal PDF", RooArgSet(mgauss, tau_decay_model));

   /*********************************************************************************************
    *    Function for efficiency
    *********************************************************************************************/

   double effSlope(0);
   if (cfg.type==B0) effSlope = -7.019e9;
   if (cfg.type==Lambda_b) effSlope = +1.316e9;
   effSlope*=.5;
   RooFormulaVar eff("eff",("1.0 + TMath::Abs(t)*"+toString(effSlope)).c_str(),t) ;

   //RooRealVar eff0("eff0","eff0",1);
   //RooRealVar eff1("eff1","eff1",cfg.type==B0?-7.019e9:+1.316e9);
   //RooPolynomial eff("eff","eff",t,RooArgSet(eff0, eff1));

   //RooEffProd modelEff("modelEff","model with efficiency", full_model, eff) ;


   /*********************************************************************************************
    *    Complete 2D model
    *********************************************************************************************/
 
   RooRealVar n_prompt("n_prompt","Number of prompt background events",1000, 0, 50000);
   RooRealVar n_nonprompt("n_nonprompt","Number of non-prompt background events",300, 0, 20000);
   RooRealVar n_signal("n_signal","Number of signal events",4000, 10, 20000);
   if(cfg.no_bkgd && !cfg.eff) RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal), RooArgList(n_signal));
   else if(cfg.no_bkgd && cfg.eff)
   {
       RooAddPdf intermed_model("intermed_model", "Full 2D model", RooArgList(signal), RooArgList(n_signal));
       RooEffProd full_model("full_model", "Full 2D model", intermed_model, eff);
   }
   else if(!cfg.nonprompt_bk) RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal, prompt), RooArgList(n_signal, n_prompt));
   else if(!cfg.prompt_bk) RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal, nonprompt), RooArgList(n_signal, n_nonprompt));
   else if(cfg.eff)
   {
       RooAddPdf intermed_model("intermed_model", "Full 2D model", RooArgList(signal, prompt, nonprompt), RooArgList(n_signal, n_prompt, n_nonprompt));
       RooEffProd full_model("full_model", "Full 2D model", intermed_model, eff);
   }
   else RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal, prompt, nonprompt), RooArgList(n_signal, n_prompt, n_nonprompt));


   /*********************************************************************************************
    *    Do the fitting
    *********************************************************************************************/ 

   t.setRange("all", -1.0e-12, 15e-12);

   full_model.fitTo(data, Save(kTRUE), NumCPU(6));
   
   RooAbsReal *cdf=mgauss.createCdf(mass);
   two_sigma_upper=cdf->findRoot(mass, cfg.mass_low, cfg.mass_high, 0.97725);
   two_sigma_lower=cdf->findRoot(mass, cfg.mass_low, cfg.mass_high, 0.02275);
   three_sigma_upper=cdf->findRoot(mass, cfg.mass_low, cfg.mass_high, 0.99865);
   three_sigma_lower=cdf->findRoot(mass, cfg.mass_low, cfg.mass_high, 1.0-0.99865);


   /*********************************************************************************************
    *     Plot results
    *********************************************************************************************/  

   TH2* hd = data.createHistogram("hd",t,Binning(60),YVar(mass,Binning(20)));
   TH2* hf = full_model.createHistogram("hf",t,Binning(60),YVar(mass,Binning(20))) ;

   c1=new TCanvas(name_c1, cfg.title1);
   c1->Divide(2) ;
   c1->cd(1); 
   gPad->SetLogz();
   hd->Draw("surf") ;
   c1->cd(2); 
   gPad->SetLogz();
   hf->Draw("surf") ;

   c2=new TCanvas(name_c2, cfg.title2, 1100, 550);
   c2->Divide(2) ;
   c2->cd(1);
   RooPlot* framex = mass.frame(Title("Mass projection")) ;
   data.plotOn(framex) ;
   full_model.plotOn(framex, Range("all"), LineColor(kBlue));
   full_model.plotOn(framex, Components(signal), LineColor(kRed));
   if(cfg.prompt_bk) full_model.plotOn(framex, Components(prompt), LineStyle(kDashed), LineColor(kBlack));
   if(cfg.nonprompt_bk) full_model.plotOn(framex, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack));
   framex->addObject(writeTLatex("mass: " +
   roundToString(mass_peak.getVal(), 3) + " #pm " +
   roundToString(mass_peak.getError()+0.0004, 3) + " GeV/c^{2}", 0.32, 0.87, 0.05));
   framex->Draw();
   c2->cd(2);
   gPad->SetLogy(); 
   RooPlot* framey = t.frame(Title("Lifetime projection")) ;
   data.plotOn(framey) ;
   full_model.plotOn(framey, Range("all"), LineColor(kBlue)) ;
   full_model.plotOn(framey, Components(signal), LineColor(kRed));
   if(cfg.prompt_bk) full_model.plotOn(framey, Components(prompt), LineStyle(kDashed), LineColor(kBlack));
   if(cfg.nonprompt_bk) full_model.plotOn(framey, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack));

   framey->addObject(writeTLatex("#tau: " +
   roundToString(tau.getVal()*1e12, 3) + " #pm " +
   roundToString(tau.getError()*1e12, 3) + " ps", 0.5,
   0.87, 0.05));
   framey->Draw();
    
   c3=new TCanvas(name_c3, cfg.title3, 1100, 550);
   c3->cd();
   c3->Divide(3);
   c3->cd(1);
   gPad->SetLogy();
   mass.setRange("sbl", cfg.mass_low, three_sigma_lower);
   RooPlot* framesbl = t.frame(Title("Lower sideband"), Range("sbl")) ;
   data.plotOn(framesbl, CutRange("sbl")) ;
   full_model.plotOn(framesbl, ProjectionRange("sbl")) ;
   full_model.plotOn(framesbl, Components(signal), LineColor(kRed), ProjectionRange("sbl"));
   if(cfg.prompt_bk) full_model.plotOn(framesbl, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sbl"));
   if(cfg.nonprompt_bk) full_model.plotOn(framesbl, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sbl"));
   framesbl->Draw();

   c3->cd(2);
   gPad->SetLogy();
   mass.setRange("sig", two_sigma_lower, two_sigma_upper);
   RooPlot* framesig = t.frame(Title("Signal region"),Range("sig")) ;
   data.plotOn(framesig, CutRange("sig")) ;
   full_model.plotOn(framesig, ProjectionRange("sig")) ;
   full_model.plotOn(framesig, Components(signal), LineColor(kRed), ProjectionRange("sig"));
   if(cfg.prompt_bk) full_model.plotOn(framesig, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sig"));
   if(cfg.nonprompt_bk) full_model.plotOn(framesig, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sig"));
   framesig->Draw();

   c3->cd(3);
   gPad->SetLogy();
   mass.setRange("sb", three_sigma_upper, cfg.mass_high);
   RooPlot* framesbh = t.frame(Title("Upper sideband"),Range("sb")) ;
   data.plotOn(framesbh, CutRange("sb")) ;
   full_model.plotOn(framesbh, ProjectionRange("sb")) ;
   full_model.plotOn(framesbh, Components(signal), LineColor(kRed), ProjectionRange("sb"));
   if(cfg.prompt_bk) full_model.plotOn(framesbh, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sb"));
   if(cfg.nonprompt_bk) full_model.plotOn(framesbh, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sb"));
   framesbh->Draw();

   //cout << "Fitresults: t: "<< tau.getVal()*1e12 << " +- "<<tau.getError()*1e12<<" ps ";
   //cout << " chisq(m): " << framex->chiSquare() << " chisq(t): " << framey->chiSquare();
   //cout << " mass: " << mass_peak.getVal() << " m_sigma1: " << m_sigma1.getVal() << " m_sigma2: " << m_sigma2.getVal() << " tau_bk: " << tau_bk.getVal();
   //cout <<  " nsig: " << n_signal.getVal()<<" npr: "<<n_prompt.getVal()<<" nnpr: "<<n_nonprompt.getVal()<<endl;
   
   cfg.f->Close();
   cfg.f->Delete();
   delete cfg.f;

   cout << "FITRESULTS: ";
   cout << tau.getVal()*1e12 << " " << tau.getError()*1e12 << " ";
   cout << mass_peak.getVal() << " " << mass_peak.getError() << " ";
   cout << m_sigma1.getVal() << " " << m_sigma1.getError() << " ";
   cout << m_sigma2.getVal() << " " << m_sigma2.getError() << " ";
   cout << frac_m_gauss.getVal() << " " << frac_m_gauss.getError() << " ";
   cout << core_frac.getVal() << " " << core_frac.getError() << " ";
   cout << tau_bk.getVal()*1e12 << " " << tau_bk.getError()*1e12 << " ";
   cout << n_nonprompt.getVal() << " " << n_nonprompt.getError() << " ";
   cout << n_prompt.getVal() << " " << n_prompt.getError() << " ";
   cout << n_signal.getVal() << " " << n_signal.getError() << " ";
   cout << prompt_p1.getVal() << " " << prompt_p1.getError() << " ";
   cout << prompt_sigma_core.getVal() << " " << prompt_sigma_core.getError() << " ";
   cout << prompt_sigma_tail.getVal() << " " << prompt_sigma_tail.getError() << " ";
   cout << endl;
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

