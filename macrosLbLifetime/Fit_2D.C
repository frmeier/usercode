#include "Configure.C"
#include <iomanip>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
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

const int set_nCPU(6);
const double timeFactor(1.0e-12); // make it 1e-12 if your data is in s, or take 1 if it is in ps
const double invTimeFactor = 1./timeFactor;

TTree *tree;
ofstream result;
ConfigData cfg;

Double_t yield_nonprompt, yield_sig, yield_prompt;
Double_t two_sigma_upper, two_sigma_lower, three_sigma_upper, three_sigma_lower;
TCanvas *c1, *c2, *c3;
TString name_c1("M_vs_lt"), name_c2("projections"), name_c3("slices");

void Fit_2D(TString type, TString saveAs = "") // saveAs: without extension
{
   if(!Configure(type, cfg)) return;
   tree=(TTree*) cfg.f->Get("fittree");

   RooRealVar mass("mass", "mass [GeV/c^{2}]", cfg.mass_low, cfg.mass_high);
   RooRealVar t("t", "lifetime [ps]", -1.0*timeFactor, 15*timeFactor);
   RooRealVar tE("tE", "lifetime error", 1e-5*timeFactor, 1*timeFactor);
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
      RooRealVar prompt_sigma_core("prompt_sigma_core","Gauss sigma core of prompt background", 0.1*timeFactor, 0.05*timeFactor, 0.12*timeFactor);
      RooRealVar prompt_sigma_tail("prompt_sigma_tail","Gauss sigma tail of prompt background", 0.2*timeFactor, 0.12*timeFactor, 0.9*timeFactor);
      RooGaussian prompt_core("prompt_core","Gauss core for prompt background", t, reso_mean, prompt_sigma_core);
      RooGaussian prompt_tail("prompt_tail","Gauss tail for prompt background", t, reso_mean, prompt_sigma_tail);
      RooRealVar core_frac("core_frac","Fraction of core in prompt background",0.8, 0.5, 0.99);
      RooAddPdf prompt_lt("prompt_lt","Double Gauss for prompt background lifetime",RooArgList(prompt_core, prompt_tail), core_frac);
      RooProdPdf prompt("prompt","Prompt 2D background model", RooArgSet(prompt_mass, prompt_lt));   

      if(!cfg.PerEventError)
      {
	  RooGaussModel reso_core("reso_core", "First Gauss of resolution function", t, reso_mean, prompt_sigma_core);
	  RooGaussModel reso_tail("reso_tail", "Second Gauss of resolution function", t, reso_mean, prompt_sigma_tail);
	  RooAddModel tau_reso("tau_reso", "Double Gauss of resolution function", RooArgList(reso_core, reso_tail), core_frac);
	  //RooAddModel tau_reso("tau_reso", "Double Gauss of resolution function", reso_core, core_frac);
      }
   } else {
      RooRealVar tau_reso_sigma("tau_reso_sigma","Sigma of lifetime resolution function",0.01*timeFactor, 0.001*timeFactor,0.9*timeFactor) ;
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
      RooRealVar tau_bk("tau_bk","Background lifetime",1.2*timeFactor, .005*timeFactor, 2*timeFactor);
      if(!cfg.PerEventError)
      {
	  RooDecay nonprompt_lt("nonprompt_lt","decay (x) double Gauss for background lifetime", t, tau_bk, tau_reso, RooDecay::SingleSided) ;
      }
      else
      {
	  RooRealVar npr_reso_bias("npr_reso_bias","Gauss sigma of signal resolution", 0, -10, 10); // has meaning of a pull
	  RooRealVar npr_reso_sigma("npr_reso_sigma","Gauss sigma of signal resolution", 1, 0.1, 10);
	  RooGaussModel npr_reso("npr_reso", "Gauss of resolution function", t, npr_reso_bias, npr_reso_sigma, tE);
	  RooDecay nonprompt_lt("nonprompt_lt","decay (x) double Gauss for background lifetime", t, tau_bk, npr_reso, RooDecay::SingleSided) ;
      }
      if(cfg.single_lt){
	 cout<<"Nonprompt background has one lifetime component."<<endl;
         RooProdPdf nonprompt("nonprompt", "Non-prompt background PDF", RooArgSet(nonprompt_mass, nonprompt_lt));   
      } else {
	 cout<<"Nonprompt background has two lifetime components."<<endl;
         RooRealVar tau_bk2("tau_bk2","Second background lifetime",0.2*timeFactor, .005*timeFactor, 2*timeFactor);
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
   //RooRealVar m_sigma1("m_sigma1","Gauss core sigma for signal mass",0.007,0.001,0.012); // bessere Werte
   //RooRealVar m_sigma2("m_sigma2","Gauss tail sigma for signal mass",0.02, 0.010, 0.035); //  dito
   RooRealVar m_sigma1("m_sigma1","Gauss core sigma for signal mass",0.007,0.001,0.012); // Hadi, alter Zustand
   RooRealVar m_sigma2("m_sigma2","Gauss tail sigma for signal mass",0.02);//, 0.010, 0.035); // Hadi, alter Zustand
   RooGaussian m_gauss1("m_gauss1","Core Gauss for signal mass", mass, mass_peak, m_sigma1);
   RooGaussian m_gauss2("m_gauss2","Tail Gauss for signal mass", mass, mass_peak, m_sigma2);
   RooRealVar frac_m_gauss("frac_m_gauss","Fraction of tail Gauss for signal mass",0.2, 0.1, 0.99);
   if (cfg.single_sig)
       RooGaussian mgauss("mgauss","Single gauss for signal mass", mass, mass_peak, m_sigma1);
   else
       RooAddPdf mgauss("mgauss","Double gauss for signal mass",RooArgList(m_gauss1, m_gauss2),RooArgList(frac_m_gauss));

   // signal lifetime
   RooRealVar tau("tau","Lambda_b lifetime",1.5*timeFactor, 0.1*timeFactor, 2*timeFactor);
   if(!cfg.PerEventError)
   {
       RooDecay tau_decay_model("tau_decay_model","decay (x) double Gauss for signal lifetime", t, tau, tau_reso, RooDecay::SingleSided) ;
   }
   else
   {
       RooRealVar sig_reso_bias("sig_reso_bias","Gauss sigma of signal resolution", 0, -10, 10); // has meaning of a pull
       RooRealVar sig_reso_sigma("sig_reso_sigma","Gauss sigma of signal resolution", 1, 0.1, 10);
       RooGaussModel sig_reso("sig_reso", "Gauss of resolution function", t, sig_reso_bias, sig_reso_sigma, tE);
       RooDecay tau_decay_model("tau_decay_model","decay (x) double Gauss for signal lifetime", t, tau, sig_reso, RooDecay::SingleSided) ;
   }
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
   // extract number of parameters in use
   RooArgSet dummy;
   const int Npars = full_model->getParameters(dummy).getSize();
   cout << "Npars: " << Npars << endl;


   /*********************************************************************************************
    *    Do the fitting
    *********************************************************************************************/ 

   t.setRange("all", -1.0*timeFactor, 15*timeFactor);

   if (!cfg.PerEventError)
       full_model.fitTo(data, Save(kTRUE), NumCPU(set_nCPU));
   else
       full_model.fitTo(data, ConditionalObservables(tE), Save(kTRUE), NumCPU(set_nCPU));
   
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

   // ---------------------------------------------- 2d plots
   c1=new TCanvas(name_c1, cfg.title1);
   c1->Divide(2) ;
   c1->cd(1); 
   gPad->SetLogz();
   hd->Draw("surf") ;
   c1->cd(2); 
   gPad->SetLogz();
   hf->Draw("surf") ;
   if (saveAs.Length() !=0 ) c1->SaveAs(saveAs+"_2d.pdf");

   c2=new TCanvas(name_c2, cfg.title2, 1100, cfg.m_t_plot_small?550:850);
   c2->Divide(2) ;

   // ---------------------------------------------- Main plots: mass and time projection
   c2->cd(1);
   if (cfg.ratioplots) 
   {
       TPad *m_pad1 = new TPad("m_pad1","m_pad1",0,0.3,1,1);
       m_pad1->SetBottomMargin(0.002);
       m_pad1->SetLeftMargin(0.12);
       m_pad1->SetRightMargin(0.004);
       m_pad1->Draw();
       m_pad1->cd();
   }

   RooPlot* framex = mass.frame(Title("Mass projection")) ;
   data.plotOn(framex, Name("data_m"), NumCPU(set_nCPU)) ;
   full_model.plotOn(framex, Range("all"), LineColor(kBlue), Name("model_m"), NumCPU(set_nCPU));
   full_model.plotOn(framex, Components(signal), LineColor(kRed), NumCPU(set_nCPU));
   if(cfg.prompt_bk) full_model.plotOn(framex, Components(prompt), LineStyle(kDashed), LineColor(kBlack), NumCPU(set_nCPU));
   if(cfg.nonprompt_bk) full_model.plotOn(framex, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), NumCPU(set_nCPU));
   // calculate chi2
   const double chi2_m = framex->chiSquare("model_m","data_m", Npars) ;
   if (cfg.publication) 
   {
       framex->addObject(writeTLatex("CMS preliminary", 0.5, 0.8, 0.05));
       framex->GetYaxis()->SetTitleOffset(1.4);
   }
   else
   {
       framex->addObject(writeTLatex("mass: "+roundToString(mass_peak.getVal(),3)+" #pm "+roundToString(mass_peak.getError()+0.0004, 3)+" GeV/c^{2}", 0.32, 0.85, 0.05));
       framex->addObject(writeTLatex("#chi^{2}/ndof: "+roundToString(chi2_m, 3), 0.5, 0.80, 0.04));
       framex->GetYaxis()->SetTitleOffset(1.6);
   }
   framex->SetMinimum(0.01);
   framex->Draw();

   if (cfg.ratioplots)
   {
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
   }


   c2->cd(2);
   if (cfg.ratioplots) 
   {
       TPad *t_pad1 = new TPad("t_pad1","t_pad1",0,0.3,1,1);
       t_pad1->SetBottomMargin(0.002);
       t_pad1->SetLeftMargin(0.12);
       t_pad1->SetRightMargin(0.004);
       t_pad1->Draw();
       t_pad1->cd();
   }

   gPad->SetLogy(); 
   RooPlot* framey = t.frame(Title("Lifetime projection")) ;
   data.plotOn(framey, Name("data_t"), NumCPU(set_nCPU), Binning(100)) ;
   full_model.plotOn(framey, Range("all"), LineColor(kBlue), Name("model_t"), NumCPU(set_nCPU)) ;
   full_model.plotOn(framey, Components(signal), LineColor(kRed), NumCPU(set_nCPU));
   if(cfg.prompt_bk) full_model.plotOn(framey, Components(prompt), LineStyle(kDashed), LineColor(kBlack), NumCPU(set_nCPU) );
   if(cfg.nonprompt_bk) full_model.plotOn(framey, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), NumCPU(set_nCPU));

   // calculate chi2
   const double chi2_t = framey->chiSquare("model_t","data_t", Npars) ;

   if (cfg.publication) 
   {
       framey->addObject(writeTLatex("CMS preliminary", 0.5, 0.8, 0.05));
       framey->GetYaxis()->SetTitleOffset(1.5);
   }
   else
   {
       framey->addObject(writeTLatex("#tau: "+roundToString(tau.getVal()*invTimeFactor, 3)+" #pm "+roundToString(tau.getError()*invTimeFactor, 3)+" ps", 0.5, 0.85, 0.05));
       framey->addObject(writeTLatex("#chi^{2}/ndof: "+roundToString(chi2_t, 3), 0.5, 0.80, 0.04));
       gPad->SetLeftMargin(0.10);
       gPad->SetRightMargin(0.003);
       framey->GetYaxis()->SetTitleOffset(1.5);
   }
   framey->Draw();

   if (cfg.ratioplots)
   {
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
   }

   if (saveAs.Length() !=0 ) c2->SaveAs(saveAs+"_m_t.pdf");
    
   //--------------- Sideband projection plots
   c3=new TCanvas(name_c3, cfg.title3, 1100, 550);
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
   mass.setRange("sbl", cfg.mass_low, three_sigma_lower);
   RooPlot* framesbl = t.frame(Title("Lower sideband projection"), Range("sbl")) ;
   data.plotOn(framesbl, CutRange("sbl"), NumCPU(set_nCPU)) ;
   full_model.plotOn(framesbl, ProjectionRange("sbl"), NumCPU(set_nCPU)) ;
   full_model.plotOn(framesbl, Components(signal), LineColor(kRed), ProjectionRange("sbl"), NumCPU(set_nCPU));
   if(cfg.prompt_bk) full_model.plotOn(framesbl, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sbl"), NumCPU(set_nCPU));
   if(cfg.nonprompt_bk) full_model.plotOn(framesbl, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sbl"), NumCPU(set_nCPU));
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

   RooHist *ratioplot_t = framesbl->residHist("h_data_CutRange[sbl]", cfg.PerEventError ? "full_model_Int[mass,tE|sbl]_Norm[mass,t,tE]" : "full_model_Int[mass|sbl]_Norm[mass,t]", true);
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
   data.plotOn(framesig, CutRange("sig"), NumCPU(set_nCPU)) ;
   full_model.plotOn(framesig, ProjectionRange("sig"), NumCPU(set_nCPU)) ;
   full_model.plotOn(framesig, Components(signal), LineColor(kRed), ProjectionRange("sig"), NumCPU(set_nCPU));
   if(cfg.prompt_bk) full_model.plotOn(framesig, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sig"), NumCPU(set_nCPU));
   if(cfg.nonprompt_bk) full_model.plotOn(framesig, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sig"), NumCPU(set_nCPU));
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

   RooHist *ratioplot_t = framesig->residHist("h_data_CutRange[sig]", cfg.PerEventError ? "full_model_Int[mass,tE|sig]_Norm[mass,t,tE]" : "full_model_Int[mass|sig]_Norm[mass,t]", true);
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
   mass.setRange("sb", three_sigma_upper, cfg.mass_high);
   RooPlot* framesbh = t.frame(Title("Upper sideband projection"),Range("sb")) ;
   data.plotOn(framesbh, CutRange("sb"), NumCPU(set_nCPU)) ;
   full_model.plotOn(framesbh, ProjectionRange("sb"), NumCPU(set_nCPU)) ;
   full_model.plotOn(framesbh, Components(signal), LineColor(kRed), ProjectionRange("sb"), NumCPU(set_nCPU));
   if(cfg.prompt_bk) full_model.plotOn(framesbh, Components(prompt), LineStyle(kDashed), LineColor(kBlack), ProjectionRange("sb"), NumCPU(set_nCPU));
   if(cfg.nonprompt_bk) full_model.plotOn(framesbh, Components(nonprompt), LineStyle(kDotted), LineColor(kBlack), ProjectionRange("sb"), NumCPU(set_nCPU));
   gPad->SetLeftMargin(0.12);
   gPad->SetRightMargin(0.004);
   framesbh->GetYaxis()->SetTitleOffset(1.2);
   framesbh->GetYaxis()->SetTitleSize(0.045);
   framesbh->GetXaxis()->SetTitleSize(0.045);
   framesbh->Draw();

   /*
    cout << "Frame objects:\n";
    for (int i=0; i<framesbh->numItems(); i++) {
	const string obj_name = framesbh->nameOf(i);
	if (obj_name=="") continue;
	cout << i << ".: " << obj_name << endl;
    }
   return;
   */

   c3->cd(3);
   TPad *sbh_pad2 = new TPad("sbh_pad2","sbh_pad2",0,0,1,0.3);
   sbh_pad2->SetTopMargin(0.002);
   sbh_pad2->SetBottomMargin(0.2);
   sbh_pad2->SetLeftMargin(0.12);
   sbh_pad2->SetRightMargin(0.004);
   sbh_pad2->Draw();
   sbh_pad2->cd();

   RooHist *ratioplot_t = framesbh->residHist("h_data_CutRange[sb]", cfg.PerEventError ? "full_model_Int[mass,tE|sb]_Norm[mass,t,tE]" : "full_model_Int[mass|sb]_Norm[mass,t]", true);
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

   //cfg.f->Close();
   //cfg.f->Delete();
   //delete cfg.f;

   cout << "FITRESULTS: ";
   cout << tau.getVal()*invTimeFactor << " " << tau.getError()*invTimeFactor << " ";
   cout << mass_peak.getVal() << " " << mass_peak.getError() << " ";
   cout << m_sigma1.getVal() << " " << m_sigma1.getError() << " ";
   cout << m_sigma2.getVal() << " " << m_sigma2.getError() << " ";
   cout << frac_m_gauss.getVal() << " " << frac_m_gauss.getError() << " ";
   cout << core_frac.getVal() << " " << core_frac.getError() << " ";
   cout << tau_bk.getVal()*invTimeFactor << " " << tau_bk.getError()*invTimeFactor << " ";
   cout << n_nonprompt.getVal() << " " << n_nonprompt.getError() << " ";
   cout << n_prompt.getVal() << " " << n_prompt.getError() << " ";
   cout << n_signal.getVal() << " " << n_signal.getError() << " ";
   cout << prompt_p1.getVal() << " " << prompt_p1.getError() << " ";
   cout << prompt_sigma_core.getVal() << " " << prompt_sigma_core.getError() << " ";
   cout << prompt_sigma_tail.getVal() << " " << prompt_sigma_tail.getError() << " ";
   cout << chi2_m << " " << chi2_t;
   cout << endl;
}

void doSomePlots(TString filestem, TString arguments)
{
    vector<TString> todolist;
    todolist.push_back(TString(""));
    todolist.push_back(TString("_phiplus"));
    todolist.push_back(TString("_phiminus"));
    todolist.push_back(TString("_phiLTpihalve"));
    todolist.push_back(TString("_phiGTpihalve"));
    todolist.push_back(TString("_etaplus"));
    todolist.push_back(TString("_etaminus"));
    todolist.push_back(TString("_PVlo"));
    todolist.push_back(TString("_PVhi"));
    todolist.push_back(TString("_ptlo"));
    todolist.push_back(TString("_pthi"));
    todolist.push_back(TString("_runA"));
    todolist.push_back(TString("_runB"));
    todolist.push_back(TString("_cow"));
    todolist.push_back(TString("_sea"));
    todolist.push_back(TString("_cowA"));
    todolist.push_back(TString("_cowB"));
    todolist.push_back(TString("_seaA"));
    todolist.push_back(TString("_seaB"));
    todolist.push_back(TString("_lb"));
    todolist.push_back(TString("_lbbar"));

    for(vector<TString>::const_iterator it = todolist.begin(); it!=todolist.end(); it++)
    {
	cout << "Working on: " << *it << endl;
	Fit_2D(filestem+(*it)+".root plotsmall " + arguments, filestem+(*it));
    }

    return;
    
    //Fit_2D("../data/vrt_r548_lb_data_lb14.root");
}

void doNPlots(TString filestem, int N, TString arguments)
{
    for(int i=0; i!=N; i++)
    {
	cout << "Working on file no. " << i << endl;
	Fit_2D(filestem+TString(toString(i))+".root " + arguments, filestem+TString(toString(i)));
    }

    return;

    //Fit_2D("../data/vrt_r548_lb_data_lb14.root");
}

