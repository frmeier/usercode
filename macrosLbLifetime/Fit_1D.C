using namespace RooFit;
using namespace TString;
using namespace std;

TFile *f;
TTree *tree;
ofstream result;

#include <iomanip>

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

Double_t nsig, nbkg;
Double_t masspeak;
Double_t two_sigma_upper, two_sigma_lower, three_sigma, three_sigma_lower;
Double_t mass_low, mass_high;

TString title1, title2, title3;
TString name1("mass"), name2("sideband"), name3("lifetime");
TCanvas *c1, *c2, *c3;


Fit_1D(TString type)
{
  bool reduced_lt = true;
  bool displaced = false;
  bool MC = false;
  if(type.Contains("full",kIgnoreCase)) reduced_lt = false;
  if(type.Contains("disp",kIgnoreCase)) displaced = true;
  if(type.Contains("mc",kIgnoreCase)) MC = true;
  
  if(type.Contains("B0",kIgnoreCase)) {          // B0
     mass_low=5.16;
     mass_high=6.0;
     masspeak=5.28;
     if(type.Contains("mc", kIgnoreCase)){       // B0 MC
        if(reduced_lt){
	  if(displaced) {
	     result.open("Results_b0_mc_displaced.tex",ios::out);
	     f=new TFile("../data/vrt_r398__399_B0_MCmix_B0_B005_displ.root");
	     title1=TString("B^{0} mass, Monte Carlo, reduced lifetime, displaced vertex trigger");  
	     title2=TString("B^{0} reduced lifetime sideband, Monte Carlo, displaced vertex trigger");  
	     title3=TString("B^{0} reduced lifetime, Monte Carlo, displaced vertex trigger");  
	  } else {
	     result.open("Results_b0_mc_barrel.tex",ios::out);
	     f=new TFile("../data/vrt_r398__399_B0_MCmix_B0_B005_barrel.root");
	     title1=TString("B^{0} mass, Monte Carlo, reduced lifetime, barrel trigger");  
	     title2=TString("B^{0} reduced lifetime sideband, Monte Carlo, barrel trigger");  
	     title3=TString("B^{0} reduced lifetime, Monte Carlo, barrel trigger");  
	  }
        } else {
	   if(displaced) {
	      f=new TFile("../data/vrt_r398__399_B0_MCmix_B0_B005_displ.root");
	      title1=TString("B^{0} mass, Monte Carlo, displaced vertex trigger");  
	      title2=TString("B^{0} lifetime sideband, Monte Carlo, displaced vertex trigger");  
	      title3=TString("B^{0} lifetime, Monte Carlo, displaced vertex trigger"); 
	   } else {
       	      f=new TFile("../data/vrt_r398__399_B0_MCmix_B0_B005_barrel.root");
	      title1=TString("B^{0} mass, Monte Carlo, full lifetime, barrel trigger");  
	      title2=TString("B^{0} lifetime sideband, Monte Carlo, barrel trigger");  
	      title3=TString("B^{0} lifetime, Monte Carlo, barrel trigger"); 
	   } 
	}
     } else {                                    // B0 data
        if(reduced_lt) {
	  if(displaced) {
	     result.open("Results_b0_data_displaced.tex",ios::out);
	     f=new TFile("../data/vrt_r393__397_B0_data_B005_displ.root");
	     title1=TString("B^{0} mass, data, reduced lifetime, displaced vertex trigger");  
	     title2=TString("B^{0} reduced lifetime sideband, data, displaced vertex trigger");  
	     title3=TString("B^{0} reduced lifetime, data, displaced vertex trigger");  
	  } else {
	     result.open("Results_b0_data_barrel.tex",ios::out);
	     f=new TFile("../data/vrt_r393__397_B0_data_B005_barrel.root");
	     title1=TString("B^{0} mass, data, reduced lifetime, barrel trigger");  
	     title2=TString("B^{0} reduced lifetime sideband, data, barrel trigger");  
	     title3=TString("B^{0} reduced lifetime, data, barrel trigger");  
	  }
        } else {
	   if(displaced) {
	     f=new TFile("../data/vrt_r393__397_B0_data_B005_displ.root");
	     title1=TString("B^{0} mass, data, displaced vertex trigger");  
	     title2=TString("B^{0} lifetime sideband, data, displaced vertex trigger");  
	     title3=TString("B^{0} lifetime, data, displaced vertex trigger"); 
	   } else {
 	      f=new TFile("../data/vrt_r393__397_B0_data_barrel.root");
 	      title1=TString("B^{0} mass, data, barrel trigger");  
	      title2=TString("B^{0} lifetime sideband, data, barrel trigger");  
	      title3=TString("B^{0} lifetime, data, barrel trigger"); 
	   } 
	}
     }
  } else {                                       // Lambda_b
     mass_low=5.4;
     mass_high=6.0;
     masspeak=5.62;
     if(type.Contains("mc", kIgnoreCase)){       // Lambda_b MC
        if(reduced_lt){
	  if(displaced) {
	     f=new TFile("../data/vrt_r385__387_Lb_MCmix_Lb_lb11_displ.root");
	     title1=TString("#Lambda_{b} mass, Monte Carlo, reduced lifetime, displaced vertex trigger");  
	     title2=TString("#Lambda_{b} reduced lifetime sideband, Monte Carlo, displaced vertex trigger");  
	     title3=TString("#Lambda_{b} reduced lifetime, Monte Carlo, displaced vertex trigger");  
	  } else {
	     f=new TFile("../data/vrt_r385__387_Lb_MCmix_Lb_lb11_barrel.root");
	     title1=TString("#Lambda_{b} mass, Monte Carlo, reduced lifetime, barrel trigger");  
	     title2=TString("#Lambda_{b} reduced lifetime sideband, Monte Carlo, barrel trigger");  
	     title3=TString("#Lambda_{b} reduced lifetime, Monte Carlo, barrel trigger"); 
	  }
        } else {
	   if(displaced) {
	      f=new TFile("../data/vrt_r385__387_Lb_MCmix_Lb_lb11_displ.root");
	      title1=TString("#Lambda_{b} mass, Monte Carlo, full lifetime, displaced vertex trigger");  
	      title2=TString("#Lambda_{b} lifetime sideband, Monte Carlo, displaced vertex trigger");  
	      title3=TString("#Lambda_{b} lifetime, Monte Carlo, displaced verrtex trigger");  
	   } else {
              f=new TFile("../data/vrt_r385__387_Lb_MCmix_Lb_lb11_barrel.root");
	      title1=TString("#Lambda_{b} mass, Monte Carlo, full lifetime, barrel trigger");  
	      title2=TString("#Lambda_{b} lifetime sideband, Monte Carlo, barrel trigger");  
	      title3=TString("#Lambda_{b} lifetime, Monte Carlo, barrel trigger"); 
	   } 
	}
     } else {                                    // Lambda_b data
        if(reduced_lt) {
	   if(displaced) {
	      f=new TFile("../data/vrt_r380__384_Lb_data_lb11_displ.root");
	      title1=TString("#Lambda_{b} mass, data, reduced lifetime, displaced vertex trigger");  
	      title2=TString("#Lambda_{b} reduced lifetime sideband, data, displaced vertex trigger");  
	      title3=TString("#Lambda_{b} reduced lifetime, data, displaced vertex trigger");  
	   } else {
	     f=new TFile("../data/vrt_r380__384_Lb_data_lb11_barrel.root");
	      title1=TString("#Lambda_{b} mass, data, reduced lifetime, barrel trigger");  
	      title2=TString("#Lambda_{b} reduced lifetime sideband, data, barrel trigger");  
	      title3=TString("#Lambda_{b} reduced lifetime, data, barrel trigger"); 
	   }
        } else {
	   if(displaced) {
	      f=new TFile("../data/vrt_r380__384_Lb_data_lb11_displ.root");
	      title1=TString("#Lambda_{b} mass, data, full lifetime, displaced vertex trigger");  
	      title2=TString("#Lambda_{b} lifetime sideband, data, displaced vertex trigger");  
	      title3=TString("#Lambda_{b} lifetime, data, displaced vertex trigger");  
	   } else {
	      f=new TFile("../data/vrt_r380__384_Lb_data_lb11_barrel.root");
	      title1=TString("#Lambda_{b} mass, data, barrel trigger");  
	      title2=TString("#Lambda_{b} lifetime sideband, data, barrel trigger");  
	      title3=TString("#Lambda_{b} lifetime, data, barrel trigger"); 
	   } 
	}
     }
   }
   tree=(TTree*)f->Get("fittree");

   c1=new TCanvas(name1,title1);
   c1->cd();

   RooRealVar mass("mass","mass",mass_low,mass_high);
   RooRealVar t("t","lifetime",-1e-12,15e-12);
   RooRealVar tRed("tRed","Reduced lifetime",-1.0e-12,15e-12);
   RooDataSet data("data","data",tree,RooArgSet(mass,t,tRed));

   RooPlot *frame1=mass.frame(mass_low, mass_high);

   RooRealVar peak("peak","Gauss mean",masspeak, masspeak-0.1,masspeak+0.1);
   RooRealVar sigma1("sigma1","sigma1",0.005,0.003,0.01);
   RooRealVar sigma2("sigma2","sigma2",0.03, 0.02, 0.05);
   RooGaussian gauss1("gauss1","Gauss 1",mass,peak,sigma1);
   RooGaussian gauss2("gauss2","Gauss 2",mass,peak,sigma2);
   RooRealVar frac_gauss("frac_gauss","Fraction of 2nd Gauss",0.5);
   RooAddPdf dgauss("dgauss","Double gauss",RooArgList(gauss1, gauss2),RooArgList(frac_gauss));
   RooRealVar yield_dgauss("yield_dgauss","Yield of double Gauss",5000, 500,50000);
   RooRealVar Bspeak("Bspeak","Gauss mean of Bs background", 5.365, 5.355, 5.37);
   RooRealVar Bs_sigma("Bs_sigma","Bs sigma",0.005,0.001,0.02);
   RooGaussian Bs_gauss("Bs_gauss","Gauss of Bs background", mass, Bspeak, Bs_sigma);
   RooRealVar Bs_yield("Bs_yield","Yield of Bs background",100,0,1000);
   RooRealVar bkgd_yield("bkgd_yield","Yield of linear Background",3000,0,50000);
   RooRealVar p1("p1","p1",0, -1,1);
   RooPolynomial bkgd("bkgd", "Linear function", mass, RooArgList(p1));
   if(type.Contains("B0",kIgnoreCase)) 
         RooAddPdf sum("sum","Double gauss + linear",RooArgList(dgauss, bkgd, Bs_gauss),RooArgList(yield_dgauss, bkgd_yield, Bs_yield));
   else RooAddPdf sum("sum","Double gauss + linear",RooArgList(dgauss, bkgd),RooArgList(yield_dgauss, bkgd_yield));

   Bs_sigma=0.017;
   Bs_sigma.setConstant(kTRUE);
   RooFitResult* fitres =sum.fitTo(data,Save(kTRUE));
   frame1->addObject(writeTLatex("Mass: " +
   roundToString(peak.getVal(), 4) + " #pm " +
   roundToString(sigma1.getError(), 4) + " GeV/c^{2}", 0.3, 0.8, 0.06));

   data.plotOn(frame1);
   sum.plotOn(frame1,LineColor(kRed));
   sum.plotOn(frame1, Components(bkgd),LineStyle(kDashed));
   sum.plotOn(frame1, Components(dgauss));
   sum.paramOn(frame1);
   frame1->SetTitle(title1);
   frame1->Draw();

   RooPlot *frame2=mass.frame(mass_low, mass_high);
   RooAbsReal *cdf=dgauss.createCdf(mass);
   two_sigma_upper=cdf->findRoot(mass,mass_low, mass_high, 0.97725);
   two_sigma_lower=cdf->findRoot(mass,mass_low, mass_high, 0.02275);
   three_sigma=cdf->findRoot(mass, mass_low, mass_high, 0.99865);
   three_sigma_lower=cdf->findRoot(mass, mass_low, mass_high, 1.0-0.99865);

   mass.setRange("sideband",three_sigma, mass_high);
   mass.setRange("sideband_low",mass_low,three_sigma_lower);
   mass.setRange("signal",two_sigma_lower, two_sigma_upper);
   RooAbsReal *int_sig=dgauss.createIntegral(mass, NormSet(mass), Range("signal"));
   nsig=int_sig->getVal()*yield_dgauss.getVal();
   int_sig=bkgd.createIntegral(mass, NormSet(mass), Range("signal"));
   nbkg=int_sig->getVal()*bkgd_yield.getVal();
   
   /***************************************************************************
   /* now starting sideband fit for resolution fuction
   /***************************************************************************/

   char mass_cut[25];
   if(type.Contains("B0",kIgnoreCase)) sprintf(mass_cut, "mass>%f && mass<%f", 5.4, mass_high);  // exclude Bs peak
   else sprintf(mass_cut, "mass>%f && mass<%f", three_sigma, mass_high);
   RooDataSet *reduced=data.reduce(mass_cut);

   c2=new TCanvas(name2,title2);
   c2->cd();
   gPad->SetLogy();
   if(reduced_lt) RooPlot *frame2=tRed.frame(0, 10e-12, 60);
   else RooPlot *frame2=t.frame(-1e-12, 5e-12, 60);

   RooRealVar prompt_peak("prompt_peak","Gauss mean for prompt background", 0.);
   RooRealVar prompt_sigma("prompt_sigma","Gauss sigma of prompt background", 1e-13, 0.5e-13, 4e-13);
   RooGaussian prompt("prompt","Gauss for prompt background", t, prompt_peak, prompt_sigma);

   RooRealVar tau("tau","Lambda_b lifetime",1.2e-12, 0.1e-12, 2e-12);
   RooRealVar mean_reso("mean_reso","mean of resolution function",0) ;
   RooRealVar sigma_reso("sigma_reso","sigma of resolution function",0.2e-13, 1e-15,3.0e-13) ;
   if (reduced_lt) {
      RooGaussModel gauss_reso("gauss_reso","Gauss of resolution function", tRed, mean_reso, sigma_reso) ;
      RooDecay background_model("background_model","decay (x) gauss", tRed, tau, gauss_reso, RooDecay::SingleSided) ;
      tRed.setRange("sideband", 0.0e-13, 10e-12);
   } else { 
      RooGaussModel gauss_reso("gauss_reso","Gauss of resolution function",t,mean_reso,sigma_reso) ;
      RooDecay bk_decay_model("bk_decay_model","decay (x) gauss", t, tau, gauss_reso, RooDecay::SingleSided) ;
      t.setRange("sideband",-1e-12,10e-12);
      RooRealVar prompt_frac("prompt_frac","Fraction of prompt background",0.2, 0.05, 0.99);
      RooAddPdf background_model("background_model","decay (x) gauss + gauss",
			RooArgList(bk_decay_model, prompt), prompt_frac);
   }

   RooFitResult* fitres2 =background_model.fitTo(*reduced, Range("sideband"), Save(kTRUE));
   reduced->plotOn(frame2);
   background_model.plotOn(frame2);
   frame2->SetTitle(title2);
   frame2->Draw();

   /***************************************************************************
   /* and finally the signal region
   /***************************************************************************/

   c3=new TCanvas(name3, title3);
   c3->cd();
   gPad->SetLogy();
   sprintf(mass_cut, "mass>%f && mass<%f", two_sigma_lower, two_sigma_upper);
   RooDataSet *signal=data.reduce(mass_cut);
   if(reduced_lt) RooPlot *frame3=tRed.frame(0, 15e-12, 60);
   else RooPlot *frame3=t.frame(-1e-12, 15e-12, 60); 
   frame3->SetTitle(title3);

   RooRealVar mean_reso_sig("mean_reso_sig","mean of signal resolution function",0) ;
   RooRealVar sigma_reso_sig("sigma_reso_sig","sigma of signal resolution function",0.2e-12, 0.1e-13,3.0e-13) ;

   Double_t f=1.0-nbkg/(nsig+nbkg);
   if(reduced_lt) {
      RooGaussModel gauss_reso_sig("gauss_reso_sig","Gauss of signal resolution function",tRed,mean_reso_sig,sigma_reso_sig) ;
      RooDecay decay_model("decay_model","decay (x) gauss", tRed, tau, gauss_reso_sig, RooDecay::SingleSided) ;
   } else {
      RooGaussModel gauss_reso_sig("gauss_reso_sig","Gauss of signal resolution function",t,mean_reso_sig,sigma_reso_sig) ;
      RooDecay decay_model("decay_model","decay (x) gauss", t, tau, gauss_reso_sig, RooDecay::SingleSided) ;
   }
   RooRealVar bkg_frac("bkg_frac","Background fraction",f, 0.1, 0.999);
   if(MC && displaced){
      RooRealVar n_signal("n_signal","Number of signal events",1000, 200, 20000);
      RooAddPdf full_model("full_model","decay (x) gauss + background",
			   RooArgList(decay_model),RooArgList(n_signal));
   } else {
      RooAddPdf full_model("full_model","decay (x) gauss + background",
			RooArgList(decay_model, background_model), bkg_frac);
   }

   prompt_sigma.setConstant(kTRUE);                            // background function
   sigma_reso.setConstant(kTRUE);
   //prompt_frac.setConstant(kTRUE);
   
   //   bkg_frac.setConstant(kTRUE);                         // background fraction from mass fit
   if(reduced_lt){
     if(displaced) tRed.setRange("sig_all",5.0e-13,15e-12);
     else tRed.setRange("sig_all",5.0e-13,15e-12);
   } else {
      t.setRange("sig_all",-1.0e-12,15e-12);
   }
   RooFitResult* fitres3 = full_model.fitTo(*signal,Range("sig_all"), Save(kTRUE));
   signal->plotOn(frame3, Binning(150));
   full_model.plotOn(frame3);
   full_model.plotOn(frame3, Components(background_model),LineStyle(kDashed), LineColor(kRed));
   full_model.plotOn(frame3, Components(decay_model),LineStyle(kDashed), LineColor(kBlack));
   //   full_model.paramOn(frame3, Format("NEU",AutoPrecision(2)));

   Double_t lt=tau.getVal();
   Double_t lte=tau.getError();
   Double_t fr=bkg_frac.getVal();
   printf("Background fraction = %f\n",fr);
   printf("B0 lifetime = %e +- %e\n", lt, lte);
   frame3->addObject(writeTLatex("#tau: " +
   roundToString(tau.getVal()*1e12, 3) + " #pm " +
   roundToString(tau.getError()*1e12, 3) + " ps", 0.3, 0.8, 0.06));
   frame3->Draw();
   
   result << "\vdef{lumi:int}{\ensuremath{{"<<4.94<<"} } }"<<endl;
   result.close();
}

