using namespace RooFit;
using namespace TString;
using namespace std;

#include <iomanip>

Double_t mass_low, mass_high;
Double_t masspeak;
Double_t yield_nonprompt, yield_sig, yield_prompt;
Double_t two_sigma_upper, two_sigma_lower, three_sigma_upper, three_sigma_lower;

TFile *f;
TTree *tree;
TCanvas *c1, *c2, *c3;
TString name1("M_vs_lt"), name2("projections"), name3("slices");
TString title1, title2, title3;

TLatex* writeTLatex(std::string text, double x, double y, double size = 0.04)
{
    TLatex* txt = new TLatex(x,y,text.c_str());
    txt->SetTextSize(size);
    txt->SetNDC(true);
    return txt;
}

string roundToString(double v, streamsize precision)
{
    ostringstream oss;
    oss << setprecision(precision) << fixed << v;
    return oss.str();
}

Fit_2D(TString type)
{
  bool no_prompt = false;
  bool no_nonprompt = false;
  bool no_bkgd = false;
  bool reduced_lt = true;
  bool displaced = false;
  if(type.Contains("noprompt",kIgnoreCase)) no_prompt = true;
  if(type.Contains("nonon",kIgnoreCase)) no_nonprompt = true;
  if(type.Contains("nobk",kIgnoreCase)) no_bkgd = true;
  if(type.Contains("full",kIgnoreCase)) reduced_lt = false;
  if(type.Contains("disp",kIgnoreCase)) displaced = true;

  if(type.Contains("B0",kIgnoreCase)) {          // B0
     mass_low=5.16;
     mass_high=6.0;
     masspeak=5.3;
     if(type.Contains("mc", kIgnoreCase)){       // B0 MC
        if(reduced_lt){
          if(displaced) {
             f=new TFile("data/vrt_r337_B0_MC_B004_displ_red.root");
             title1=TString("B^{0} mass vs reduced lifetime, Monte Carlo, displaced vertex trigger");  
             title2=TString("B^{0} mass and reduced lifetime projections, Monte Carlo, displaced vertex trigger");  
             title3=TString("B^{0} reduced lifetime projections, Monte Carlo, displaced vertex trigger");  
          } else {
             f=new TFile("data/vrt_r337_B0_MC_B004_red.root");
             title1=TString("B^{0} mass vs reduced lifetime, Monte Carlo, barrel trigger");  
             title2=TString("B^{0} mass and reduced lifetime projections, Monte Carlo, barrel trigger");  
             title3=TString("B^{0} reduced lifetime projections, Monte Carlo, barrel trigger");  
          }
        } else {
           if(displaced) {
              printf("No input file defined\n");
              return 0;
           }
           //f=new TFile("data/vrt_r337_372_376_377_378_MCmix.root");
           f=new TFile("data/vrt_r337_B0_MC_B004.root");
           title1=TString("B^{0} mass vs lifetime, Monte Carlo, barrel trigger");  
           title2=TString("B^{0} mass and lifetime projections, Monte Carlo, barrel trigger");  
           title3=TString("B^{0} lifetime projections, Monte Carlo, barrel trigger");  
        }
     } else {                                    // B0 data
        if(reduced_lt) {
          if(displaced) {
             f=new TFile("data/vrt_r355_359_data_B004_displ_red.root");
             title1=TString("B^{0} mass vs reduced lifetime, data, displaced vertex trigger");  
             title2=TString("B^{0} mass and reduced lifetime projections, data, displaced vertex trigger");  
             title3=TString("B^{0} reduced lifetime projections, data, displaced vertex trigger");  
          } else {
             f=new TFile("data/vrt_r355_359_data_B004_red.root");
             title1=TString("B^{0} mass vs reduced lifetime, data, barrel trigger");  
             title2=TString("B^{0} mass and reduced lifetime projections, data, barrel trigger");  
             title3=TString("B^{0} reduced lifetime projections, data, barrel trigger");  
          }
        } else {
           if(displaced) {
              printf("No input file defined\n");
              return 0;
           }
           f=new TFile("data/vrt_r355_359_data_B004.root");
           title1=TString("B^{0} mass vs lifetime, data, barrel trigger");  
           title2=TString("B^{0} mass and lifetime projections, data, barrel trigger");  
           title3=TString("B^{0} lifetime projections, data, barrel trigger");  
        }
     }
  } else {                                       // Lambda_b
     mass_low=5.4;
     mass_high=6.0;
     masspeak=5.62;
     if(type.Contains("mc", kIgnoreCase)){       // Lambda_b MC
        if(reduced_lt){
          if(displaced) {
              printf("No input file defined\n");
              return 0;
          } else {
              printf("No input file defined\n");
              return 0;  
          }
        } else {
           if(displaced) {
              printf("No input file defined\n");
              return 0;
           }
           f=new TFile("data/vrt_MC5mix_332_346_347_438_349_cuts_notrigsel.root");
           //f=new TFile("data/vrt_r332_360_361.root");
           title1=TString("#Lambda_{b} mass vs lifetime, Monte Carlo, barrel trigger");  
           title2=TString("#Lambda_{b} mass and lifetime projections, Monte Carlo, barrel trigger");  
           title3=TString("#Lambda_{b} lifetime projections, Monte Carlo, barrel trigger");  
        }
     } else {                                    // Lambda_b data
        if(reduced_lt) {
           if(displaced) {
              printf("No input file defined\n");
              return 0;
           } else {
              printf("No input file defined\n");
              return 0; 
           }
        } else {
           if(displaced) {
              printf("No input file defined\n");
              return 0;
           }
           //f=new TFile("data/vrt_r350_354_plusMCbgr.root"); 
           f=new TFile("data/vrt_r350_354.root"); 
           title1=TString("#Lambda_{b} mass vs lifetime, data, barrel trigger");  
           title2=TString("#Lambda_{b} mass and lifetime projections, data, barrel trigger");  
           title3=TString("#Lambda_{b} lifetime projections, data, barrel trigger");  
        }
     }
   }
   tree=(TTree*)f->Get("fittree");

   RooRealVar mass("mass","mass",mass_low,mass_high);
   RooRealVar t("t","lifetime",-1.0e-12,15e-12);
   RooRealVar tRed("tRed","Reduced lifetime",-2.0e-12,15e-12);
   RooRealVar tE("tE","lifetime error",0,1e-12);
   RooDataSet data("data","data",tree,RooArgSet(mass,t,tRed,tE));


   /***************************************************************************************
   /*    Signal = Double Gauss (mass) x (Decay (x) double Gauss)(lifetime)
   /***************************************************************************************/

   //  mass peak 
   RooRealVar mass_peak("mass_peak","Gauss mean of signal mass peak", masspeak, masspeak-0.04, masspeak+0.04);
   RooRealVar m_sigma1("m_sigma1","Gauss core sigma for signal mass",0.009,0.005,0.015);
   RooRealVar m_sigma2("m_sigma2","Gauss tail sigma for signal mass",0.02, 0.01, 0.05);
   RooGaussian m_gauss1("m_gauss1","Core Gauss for signal mass",mass,mass_peak,m_sigma1);
   RooGaussian m_gauss2("m_gauss2","Tail Gauss for signal mass",mass,mass_peak,m_sigma2);
   RooRealVar frac_m_gauss("frac_m_gauss","Fraction of tail Gauss for signal mass",0.4, 0.3, 0.99);
   RooAddPdf mgauss("mgauss","Double gauss for signal mass",RooArgList(m_gauss1, m_gauss2),RooArgList(frac_m_gauss));

   // signal lifetime
   RooRealVar tau("tau","Lambda_b lifetime",1.2e-12, 0.1e-12, 2e-12);
   RooRealVar tau_reso_mean("reso_mean","Mean of resolution function",0) ;
   RooRealVar tau_reso_sigma("tau_reso_sigma","Sigma of lifetime resolution function",0.1e-13, 1e-15,9.0e-13) ;
   if(reduced_lt){
      RooGaussModel tau_reso_gauss("tau_reso_gauss","Gauss of lifetime resolution function",tRed,tau_reso_mean,tau_reso_sigma) ;
      RooDecay tau_decay_model("tau_decay_model","decay (x) double Gauss for signal lifetime", tRed, tau, tau_reso_gauss, RooDecay::SingleSided) ;
   } else {
      RooGaussModel tau_reso_gauss("tau_reso_gauss","Gauss of lifetime resolution function",t,tau_reso_mean,tau_reso_sigma) ;
      RooDecay tau_decay_model("tau_decay_model","decay (x) double Gauss for signal lifetime", t, tau, tau_reso_gauss, RooDecay::SingleSided) ;
   }
   RooProdPdf signal("signal", "Signal PDF", RooArgSet(mgauss, tau_decay_model));


   /*********************************************************************************************
   /*    Non-prompt background = 1st order polynomial (mass) x (Decay (x) double Gauss)(lifetime)
   /*    Prompt background = 1st order polynomial (mass) x double Gauss (lifetime)
   /*********************************************************************************************/

   // non-prompt mass polynomial
   RooRealVar nonprompt_p1("nonprompt_p1","Linear coefficient of nonprompt background mass polynomial", 0.1, -0.2, 0.3);
   RooPolynomial nonprompt_mass("nonprompt_mass", "Linear function for nonprompt background mass", mass, RooArgList(nonprompt_p1));
   // non-prompt lifetime model
   RooRealVar tau_bk("tau_bk","Background lifetime",1.2e-12, 5e-15, 2e-12);
   RooRealVar tau_reso_sigma_bk("tau_reso_sigma_bk","Sigma of lifetime resolution function for background",1e-14, 1e-15,1.0e-12) ;
   if(reduced_lt){
      RooGaussModel tau_reso_bkgd("tau_reso_bkgd","Gauss of lifetime resolution function for background",tRed,tau_reso_mean,tau_reso_sigma_bk) ;
      RooDecay tau_bkg_model("tau_bkg_model","decay (x) double Gauss for background lifetime", tRed, tau_bk, tau_reso_bkgd, RooDecay::SingleSided) ;
   } else {
      RooGaussModel tau_reso_bkgd("tau_reso_bkgd","Gauss of lifetime resolution function for background",t,tau_reso_mean,tau_reso_sigma_bk) ;
      RooDecay tau_bkg_model("tau_bkg_model","decay (x) double Gauss for background lifetime", t, tau_bk, tau_reso_bkgd, RooDecay::SingleSided) ;}
   RooProdPdf nonprompt("nonprompt", "Non-prompt background PDF", RooArgSet(nonprompt_mass, tau_bkg_model));   

   // prompt mass polynomial
   RooRealVar prompt_p1("prompt_p1","Linear coefficient of prompt background mass polynomial",0.1, -0.2, 5.0);
   RooPolynomial prompt_mass("prompt_mass", "Linear function for prompt background mass", mass, RooArgList(prompt_p1));
   // prompt lifetime resolution
   RooRealVar prompt_peak("prompt_peak","Gauss mean for prompt background",0.);
   RooRealVar prompt_sigma("prompt_sigma","Sigma for prompt background lifetime",1e-13,5e-14,1e-12);
   if(reduced_lt) RooGaussian prompt_gauss("prompt_gauss","Gauss for prompt background lifetime",tRed, prompt_peak, prompt_sigma);
   else RooGaussian prompt_gauss("prompt_gauss","Gauss for prompt background lifetime",t, prompt_peak, prompt_sigma);
   RooProdPdf prompt("prompt", "Prompt background model", RooArgSet(prompt_mass, prompt_gauss));


   /*********************************************************************************************
   /*    Complete 2D model
   /*********************************************************************************************/
 
   RooRealVar n_prompt("n_prompt","Number of prompt background events",1000, 0, 50000);
   RooRealVar n_nonprompt("n_nonprompt","Number of non-prompt background events",300, 0, 20000);
   RooRealVar n_signal("n_signal","Number of signal events",1000, 200, 20000);
   if(no_bkgd) RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal), RooArgList(n_signal));
   else if(no_prompt) RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal, nonprompt), RooArgList(n_signal, n_nonprompt));
   else if(no_nonprompt) RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal, prompt), RooArgList(n_signal, n_prompt));
   else RooAddPdf full_model("full_model", "Full 2D model", RooArgList(signal, prompt, nonprompt), RooArgList(n_signal, n_prompt, n_nonprompt));


   /*********************************************************************************************
   /*    Do the fitting
   /*********************************************************************************************/ 

   if(reduced_lt) {
      tRed.setRange("reduced", 5.0e-13, 15e-12);
      t.setRange("reduced",-1e-12, 15e-12);
   } else {
      t.setRange("reduced", -1.0e-12, 15e-12);
      tRed.setRange("reduced", -2.0e-12, 15e-12);
   }
   full_model.fitTo(data, Range("reduced"), Save(kTRUE));
   
   RooAbsReal *cdf=mgauss.createCdf(mass);
   two_sigma_upper=cdf->findRoot(mass,mass_low, mass_high, 0.97725);
   two_sigma_lower=cdf->findRoot(mass,mass_low, mass_high, 0.02275);
   three_sigma_upper=cdf->findRoot(mass, mass_low, mass_high, 0.99865);
   three_sigma_lower=cdf->findRoot(mass, mass_low, mass_high, 1.0-0.99865);


   /*********************************************************************************************
   /*     Plot results
   /*********************************************************************************************/  

   if(reduced_lt){
      TH2* hd = data.createHistogram("hd",tRed,Binning(60),YVar(mass,Binning(20)));
      TH2* hf = full_model.createHistogram("hf",tRed,Binning(60),YVar(mass,Binning(20))) ;
   } else {
      TH2* hd = data.createHistogram("hd",t,Binning(60),YVar(mass,Binning(20)));
      TH2* hf = full_model.createHistogram("hf",t,Binning(60),YVar(mass,Binning(20))) ;
   }
   c1=new TCanvas(name1, title1);
   c1->Divide(2) ;
   c1->cd(1); 
   gPad->SetLogz();
   hd->Draw("surf") ;
   c1->cd(2); 
   gPad->SetLogz();
   hf->Draw("surf") ;

   c2=new TCanvas(name2, title2);
   c2->Divide(2) ;
   c2->cd(1);
   RooPlot* framex = mass.frame(Title("Mass projection")) ;
   data.plotOn(framex) ;
   full_model.plotOn(framex) ;
   framex->addObject(writeTLatex("mass: " +
   roundToString(mass_peak.getVal(), 3) + " #pm " +
   roundToString(mass_peak.getError()+0.0004, 3) + " GeV/c^{2}", 0.15, 0.84, 0.06));
   framex->Draw();
   c2->cd(2);
   gPad->SetLogy();
   if(reduced_lt) RooPlot* framey = tRed.frame(Title("Lifetime projection")) ;
   else RooPlot* framey = t.frame(Title("Lifetime projection")) ;
   data.plotOn(framey) ;
   full_model.plotOn(framey) ;
   framey->addObject(writeTLatex("#tau: " +
   roundToString(tau.getVal()*1e12, 3) + " #pm " +
   roundToString(tau.getError()*1e12, 3) + " ps", 0.3,
   0.8, 0.06));
   framey->Draw();
    
   return 0;
   c3=new TCanvas(name3, title3);
   c3->cd();
   c3->Divide(3);
   c3->cd(1);
   gPad->SetLogy();
   mass.setRange("sbl", mass_low, three_sigma_lower);
   if (reduced_lt) RooPlot* framesbl = tRed.frame(Title("Lower sideband"),Range("sbl")) ;
   else RooPlot* framesbl = t.frame(Title("Lower sideband"),Range("sbl")) ;
   data.plotOn(framesbl, CutRange("sbl")) ;
   full_model.plotOn(framesbl, ProjectionRange("sbl")) ;
   framesbl->Draw();

   c3->cd(2);
   gPad->SetLogy();
   mass.setRange("sig", two_sigma_lower, two_sigma_upper);
   if(reduced_lt) RooPlot* framesig = tRed.frame(Title("Signal region"),Range("sig")) ;
   else RooPlot* framesig = t.frame(Title("Signal region"),Range("sig")) ;
   data.plotOn(framesig, CutRange("sig")) ;
   full_model.plotOn(framesig, ProjectionRange("sig")) ;
   framesig->Draw();

   c3->cd(3);
   gPad->SetLogy();
   mass.setRange("sb", three_sigma_upper, mass_high);
   if(reduced_lt) RooPlot* framesb = tRed.frame(Title("Upper sideband"),Range("sb")) ;
   else RooPlot* framesb = t.frame(Title("Upper sideband"),Range("sb")) ;
   data.plotOn(framesb, CutRange("sb")) ;
   full_model.plotOn(framesb, ProjectionRange("sb")) ;
   framesb->Draw();

   cout << "Lifetime="<< tau.getVal()*1e12 << " +- "<<tau.getError()*1e12<<" ps"<<endl;
}

