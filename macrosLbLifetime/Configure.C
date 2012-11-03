#include "TFile.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

enum Types {Lambda_b, B0, Bs, Bplus};

typedef struct ConfigData {
  bool no_bkgd;
  bool prompt_bk;
  bool nonprompt_bk;
  bool displaced;
  bool MC;
  bool eff;
  bool PerEventError;
  bool constant_prompt_fraction;
  bool single_lt;
  bool single_sig;
  bool m_t_plot_small;
  bool ratioplots;
  bool publication;
  int type;
  double mass_low;
  double mass_high;
  double masspeak;
  TString title1,title2,title3;
  TFile *f;
};


int Configure(TString arguments, ConfigData &cfg)
{
  cfg.prompt_bk = true;
  cfg.nonprompt_bk = true;
  cfg.no_bkgd = false;
  cfg.displaced = false;
  cfg.MC = false;
  cfg.eff = false;
  cfg.PerEventError=false;
  cfg.single_lt=true;
  cfg.single_sig=false;
  cfg.m_t_plot_small=false;
  cfg.ratioplots=true;
  cfg.publication=false;
  cfg.constant_prompt_fraction=true;
  cfg.type = Lambda_b;
  if(arguments.Contains("plotsmall",TString::kIgnoreCase)) cfg.m_t_plot_small= true;
  if(arguments.Contains("noratioplots",TString::kIgnoreCase)) cfg.ratioplots= true;
  if(arguments.Contains("pub",TString::kIgnoreCase)) { cfg.ratioplots= false; cfg.m_t_plot_small= true; cfg.publication=true; }
  if(arguments.Contains("noprompt",TString::kIgnoreCase)) cfg.prompt_bk = false;
  if(arguments.Contains("nonon",TString::kIgnoreCase)) cfg.nonprompt_bk = false;
  if(arguments.Contains("nobk",TString::kIgnoreCase)) {
    cfg.no_bkgd = true;
    cfg.nonprompt_bk = false;
    cfg.prompt_bk = false;
  }
  if(arguments.Contains("nobgr",TString::kIgnoreCase)) {
    cfg.no_bkgd = true;
    cfg.nonprompt_bk = false;
    cfg.prompt_bk = false;
  }
  if(arguments.Contains("eff",TString::kIgnoreCase)) cfg.eff = true;
  if(arguments.Contains("pee",TString::kIgnoreCase)) cfg.PerEventError = true;
  if(arguments.Contains("double_lt",TString::kIgnoreCase)) cfg.single_lt = false;
  if(arguments.Contains("single_sig",TString::kIgnoreCase)) cfg.single_sig = true;
  if(arguments.Contains("double_mass",TString::kIgnoreCase)) cfg.constant_prompt_fraction = false;
  if(arguments.Contains("disp",TString::kIgnoreCase)) cfg.displaced = true;
  if(arguments.Contains("mc",TString::kIgnoreCase)) cfg.MC = true;
  if(arguments.Contains("B0",TString::kIgnoreCase)) cfg.type=B0;
  if(arguments.Contains("Bs",TString::kIgnoreCase)) cfg.type=Bs;
  if(arguments.Contains("B+",TString::kIgnoreCase)) cfg.type=Bplus;

  TObjArray* tokens=arguments.Tokenize(" ");
  bool FileNameEntered=false;
  for(Int_t i=0; i<tokens->GetEntries(); i++){
    TObjString* str =dynamic_cast<TObjString*> (tokens->At(i));
    if(!str) continue;
    TString token=str->GetString();
    if(token.Contains(".root")) {
      FileNameEntered=true;
      cfg.f=new TFile(token);
    }
  }
  tokens->Clear();
  delete tokens;

  switch(cfg.type) {
  case B0:                    
     cfg.mass_low=5.16;
     //cfg.mass_high=6.0;
     cfg.mass_high=5.75;
     cfg.masspeak=5.28;
     if(cfg.MC){                                                // B0 MC
        if(cfg.displaced){
	  //	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r398__399_B0_MCmix_B0_B005_displ.root");
	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r460_B0_sigMC_B0_B006_displ.root");
           cfg.title1=TString("B^{0} mass, Monte Carlo, displaced vertex trigger");  
           cfg.title2=TString("B^{0} lifetime sideband, Monte Carlo, displaced vertex trigger");  
           cfg.title3=TString("B^{0} lifetime, Monte Carlo, displaced vertex trigger"); 
        } else {                                  //barrel
	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r460_B0_MCmix_B0_B006_barrel.root");
           cfg.title1=TString("B^{0} mass, Monte Carlo, barrel trigger");  
           cfg.title2=TString("B^{0} lifetime sideband, Monte Carlo, barrel trigger");  
           cfg.title3=TString("B^{0} lifetime, Monte Carlo, barrel trigger"); 
        }
     } else {                                    // B0 data
        if(cfg.displaced) {
	  //           if(!FileNameEntered) cfg.f=new TFile("data/vrt_r393__397_B0_data_B005_displ.root");
	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r461_B0_data_B0_B006_displ.root"); 
           cfg.title1=TString("B^{0} mass, data, displaced vertex trigger");  
           cfg.title2=TString("B^{0} lifetime sideband, data, displaced vertex trigger");  
           cfg.title3=TString("B^{0} lifetime, data, displaced vertex trigger"); 
        } else {
	   //if(!FileNameEntered) cfg.f=new TFile("data/new/vrt_r393__397_B0_data_B005_barrel.root");
 	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r461_B0_data_B0_B006_barrel.root");
           cfg.title1=TString("B^{0} mass, data, barrel trigger");  
           cfg.title2=TString("B^{0} lifetime sideband, data, barrel trigger");  
           cfg.title3=TString("B^{0} lifetime, data, barrel trigger"); 
        }
     }
     break;

  case Bplus:                    
     cfg.mass_low=5.16;
     cfg.mass_high=6.0;
     // cfg.mass_high=5.7;
     cfg.masspeak=5.28;
     if(cfg.MC){                                            // B+ MC
        if(cfg.displaced){
	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r398__399_Bp_MCmix_B0_B005_displ.root");
           cfg.title1=TString("B^{+} mass, Monte Carlo, displaced vertex trigger");  
           cfg.title2=TString("B^{+} lifetime sideband, Monte Carlo, displaced vertex trigger");  
           cfg.title3=TString("B^{+} lifetime, Monte Carlo, displaced vertex trigger"); 
        } else {                                            //barrel
           if(!FileNameEntered) cfg.f=new TFile("data/vrt_r398__399_Bp_MCmix_B0_B005_barrel.root");
           cfg.title1=TString("B^{0} mass, Monte Carlo, barrel trigger");  
           cfg.title2=TString("B^{0} lifetime sideband, Monte Carlo, barrel trigger");  
           cfg.title3=TString("B^{0} lifetime, Monte Carlo, barrel trigger"); 
        }
     } else {                                               // B+ data
        if(cfg.displaced) {
           if(!FileNameEntered) cfg.f=new TFile("data/vrt_r393__397_Bp_data_B005_displ.root");
           cfg.title1=TString("B^{+} mass, data, displaced vertex trigger");  
           cfg.title2=TString("B^{+} lifetime sideband, data, displaced vertex trigger");  
           cfg.title3=TString("B^{+} lifetime, data, displaced vertex trigger"); 
        } else {
           if(!FileNameEntered) cfg.f=new TFile("data/new/vrt_r393__397_Bp_data_B005_barrel.root");
           cfg.title1=TString("B^{+} mass, data, barrel trigger");  
           cfg.title2=TString("B^{+} lifetime sideband, data, barrel trigger");  
           cfg.title3=TString("B^{+} lifetime, data, barrel trigger"); 
        }
     }
     break; 
 
  case Bs:                  
     cfg.mass_low=5.16;
     cfg.mass_high=6.0;
     // cfg.mass_high=5.7;
     cfg.masspeak=5.28;
     if(cfg.MC){                                            // Bs MC
        if(cfg.displaced){
	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r398__399_B0_MCmix_Bs_B005_displ.root");
           cfg.title1=TString("B_{s} mass, Monte Carlo, displaced vertex trigger");  
           cfg.title2=TString("B_{s} lifetime sideband, Monte Carlo, displaced vertex trigger");  
           cfg.title3=TString("B_{s} lifetime, Monte Carlo, displaced vertex trigger"); 
        } else {                                            //barrel
           if(!FileNameEntered) cfg.f=new TFile("data/vrt_r398__399_B0_MCmix_Bs_B005_barrel.root");
           cfg.title1=TString("B_{s} mass, Monte Carlo, barrel trigger");  
           cfg.title2=TString("B_{s} lifetime sideband, Monte Carlo, barrel trigger");  
           cfg.title3=TString("B_{s} lifetime, Monte Carlo, barrel trigger"); 
        }
     } else {                                               // Bs data
        if(cfg.displaced) {
           if(!FileNameEntered) cfg.f=new TFile("data/vrt_r393__397_Bs_data_B005_displ.root");
           cfg.title1=TString("B_{s} mass, data, displaced vertex trigger");  
           cfg.title2=TString("B_{s} lifetime sideband, data, displaced vertex trigger");  
           cfg.title3=TString("B_{s} lifetime, data, displaced vertex trigger"); 
        } else {
           if(!FileNameEntered) cfg.f=new TFile("data/new/vrt_r393__397_Bs_data_B005_barrel.root");
           cfg.title1=TString("B_{s} mass, data, barrel trigger");  
           cfg.title2=TString("B_{s} lifetime sideband, data, barrel trigger");  
           cfg.title3=TString("B_{s} lifetime, data, barrel trigger"); 
        }
     }
     break;

  case Lambda_b:
  default:
     cfg.mass_low=5.4;
     cfg.mass_high=6.0;
     cfg.masspeak=5.62;
     if(cfg.MC){                                      // Lambda_b MC
        if(cfg.displaced){
           //if(!FileNameEntered) cfg.f=new TFile("data/vrt_r385__387_Lb_MCmix_Lb_lb11_displ.root");
           if(!FileNameEntered) cfg.f=new TFile("data/vrt_r456_Lb_sigMC_Lb_lb12_displ.root");
           cfg.title1=TString("#Lambda_{b} mass, Monte Carlo, displaced vertex trigger");  
           cfg.title2=TString("#Lambda_{b} lifetime sideband, Monte Carlo, displaced vertex trigger");  
           cfg.title3=TString("#Lambda_{b} lifetime, Monte Carlo, displaced verrtex trigger");  
        } else {
           if(!FileNameEntered) cfg.f=new TFile("data/vrt_r456_Lb_MCmix_Lb_lb12_barrel.root");
           cfg.title1=TString("#Lambda_{b} mass, Monte Carlo, barrel trigger");  
           cfg.title2=TString("#Lambda_{b} lifetime sideband, Monte Carlo, barrel trigger");  
           cfg.title3=TString("#Lambda_{b} lifetime, Monte Carlo, barrel trigger"); 
        }
     } else {                                    // Lambda_b data
        if(cfg.displaced) {
	   //if(!FileNameEntered) cfg.f=new TFile("data/vrt_r380__384_Lb_data_lb11_displ.root");
	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r458_Lb_data_Lb_lb12_displ.root");
           cfg.title1=TString("#Lambda_{b} mass, data, displaced vertex trigger");  
           cfg.title2=TString("#Lambda_{b} lifetime sideband, data, displaced vertex trigger");  
           cfg.title3=TString("#Lambda_{b} lifetime, data, displaced vertex trigger");  
        } else {
	   //if(!FileNameEntered) cfg.f=new TFile("data/new/vrt_r380__384_Lb_data_lb11_barrel.root");
	   if(!FileNameEntered) cfg.f=new TFile("data/vrt_r458_Lb_data_Lb_lb12_barrel.root");
           cfg.title1=TString("#Lambda_{b} mass, data, barrel trigger");  
           cfg.title2=TString("#Lambda_{b} lifetime sideband, data, barrel trigger");  
           cfg.title3=TString("#Lambda_{b} lifetime, data, barrel trigger"); 
        }
     }
  }
  if(cfg.f->IsZombie()) return 0;
  else return 1;
}
