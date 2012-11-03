#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooProdPdf.h>
#include <RooEffProd.h>
#include <RooFormulaVar.h>
#include <RooAddModel.h>

using namespace RooFit;
TRandom3 rnd;

//#define PS                  1.0E-12
#define PS                  1

// General fitting options
#define NUMBER_OF_CPU       8
#define DO_MINOS            kFALSE
// 0 - w/o DISPLAY
// 1 - w/  DISPLAY
#define DISPLAY             0

// PDF options
// 0 - static resolution, 1 - PEE resolution
#define PER_EVENT_ERROR     1
// 0/1 - switch for matching to Fit_2D_light.C
#define MATCH_MAIN_ANALYSIS 0
// 0 - no time efficiency
// 1 - pol3 + break at 6 ps
// 2 - p0(1 + p1*t + p2/(1+exp(-t/p3))), now from Frank's fit
// 3 - pol1 + break at 6 ps
#define PROD_T_EFFICIENCY   2

//Parameter sets for Lambda_b
#define MASS_MIN            5.4
#define MASS_MAX            6.0
#define MASS_PEAK           5.62
#define SOURCE              "../data/vrt_r548_lb_data_lb14_ps.root"

/*
//Parameter sets for B0
#define MASS_MIN            5.16
#define MASS_MAX            5.75
#define MASS_PEAK           5.28
#define SOURCE              "vrt_r479_B0_data_B008.root"
*/

// idx = -1, reads data from SOURCE
// idx >= 0, reads data from pseudo-[idx].txt
// do_generation = true, produce toy MC, save to pseudo-[idx].txt
// do_generation = false, fit data/toy MC
void fitter_core(int idx = -1, bool do_generation = false);

void myfitter()
{
    // Fit to data
    fitter_core();
    
    // Generate 200 sets of toy MC
    //for(int i=0;i<200;i++) fitter_core(i,true);
    
    // Fit to 200 sets of toy MC, fitting results are stored in "results.log"
    //for(int i=0;i<200;i++) fitter_core(i,false);
}

void fitter_core(int idx, bool do_generation)
{
    
    RooRealVar mass("mass","mass",MASS_MIN,MASS_MAX);
    RooRealVar t("t","t",-1.0*PS,15.0*PS);
    RooRealVar tE("tE","tE", 0.01*PS,1.0*PS);
    
    TFile *fin;
    TTree *fittree;
    RooDataSet *data;
    
    if (!do_generation) {
        if (idx<0) {
            fin = new TFile(SOURCE);
            fittree = (TTree*)fin->Get("fittree");

            data = new RooDataSet("data","data",fittree,RooArgSet(mass,t,tE));
        }else {
            data =  RooDataSet::read(TString::Format("pseudo-%03d.ascii",idx),RooArgSet(mass,t,tE));
        }
    }
        
    //-----------------------------------------------------------------
    // signal PDF 
    
    RooRealVar m_mean("m_mean","m_mean",MASS_PEAK,MASS_MIN,MASS_MAX);
    RooRealVar m_sigma1("m_sigma1","m_sigma1",0.010,0.001,0.020);
    RooRealVar m_sigma2("m_sigma2","m_sigma2",0.020,0.001,0.040);
    RooRealVar m_fraction("m_fraction","m_fraction",0.8,0.0,1.0);
    
    RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",mass,m_mean,m_sigma1);
    RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",mass,m_mean,m_sigma2);
    
    RooAddPdf pdf_m_signal("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
    
#if PER_EVENT_ERROR
    RooRealVar res_sig_mean("res_sig_mean","res_sig_mean",0.0,-1.0,1.0);
    RooRealVar res_sig_sigma("res_sig_sigma","res_sig_sigma",1.0,0.1,3.0);
    RooGaussModel res_signal("res_signal","res_signal",t,res_sig_mean,res_sig_sigma,tE);
#else
    RooRealVar res_sig_mean("res_sig_mean","res_sig_mean",0.0,-0.1*PS,0.1*PS);
    RooRealVar res_sig_sigma("res_sig_sigma","res_sig_sigma",0.1*PS,0.01*PS,0.3*PS);
    RooGaussModel res_signal("res_signal","res_signal",t,res_sig_mean,res_sig_sigma);    
#endif
  
#if PROD_T_EFFICIENCY
/*
 Model 1: pol3
 FCN=53.0782 FROM MINOS     STATUS=SUCCESSFUL     28 CALLS         225 TOTAL
 EDM=1.60441e-15    STRATEGY= 1      ERROR MATRIX ACCURATE
 EXT PARAMETER                  PARABOLIC         MINOS ERRORS
 NO.   NAME      VALUE            ERROR      NEGATIVE      POSITIVE
 1  p0           9.62934e-01   3.59552e-02  -3.59553e-02   3.59553e-02
 2  p1           7.67341e-02   7.49200e-02  -7.49204e-02   7.49204e-02
 3  p2          -3.26842e-02   3.75762e-02  -3.75764e-02   3.75764e-02
 4  p3           2.73130e-03   4.94998e-03  -4.95000e-03   4.95000e-03

 Model 2: pol1 + 1/(1+exp)
 FCN=49.1975 FROM MINOS     STATUS=SUCCESSFUL    127 CALLS         776 TOTAL
 EDM=1.79413e-06    STRATEGY= 1      ERROR MATRIX ACCURATE
 EXT PARAMETER                  PARABOLIC         MINOS ERRORS
 NO.   NAME      VALUE            ERROR      NEGATIVE      POSITIVE
 1  p0           1.05293e+00   3.29172e-02  -3.16440e-02   3.48577e-02
 2  p1          -3.16776e-02   1.49128e-02  -1.52634e-02   1.46539e-02
 3  p2          -3.86317e-01   1.59835e-01  -1.79220e-01   1.52816e-01
 4  p3          -1.27957e-01   7.87220e-02  -1.07349e-01   6.92176e-02
 
 Model 3: pol1
 FCN=56.0949 FROM MINOS     STATUS=SUCCESSFUL     12 CALLS          85 TOTAL
 EDM=8.75213e-21    STRATEGY= 1      ERROR MATRIX ACCURATE
 EXT PARAMETER                  PARABOLIC         MINOS ERRORS
 NO.   NAME      VALUE            ERROR      NEGATIVE      POSITIVE
 1  p0           1.00570e+00   2.20531e-02  -2.20531e-02   2.20531e-02
 2  p1          -1.37914e-02   1.19529e-02  -1.19529e-02   1.19529e-02
*/
    RooRealVar tau("tau","tau",1.5*PS, 0.1*PS, 3*PS);
    RooDecay pdf_t_signal_raw("pdf_t_signal_raw","pdf_t_signal_raw",t,tau,res_signal,RooDecay::SingleSided);
    
    if (PS == 1)
    {
	RooFormulaVar eff_model1("eff_model1","t<6E-12?(9.62934e-01+7.67341e-02*t-3.26842e-02*t*t*1E24+2.73130e-03*t*t*t*1E36)*0.9:8.366682e-01*0.9",RooArgList(t));
	//RooFormulaVar eff_model2("eff_model2","t<6E-12?(1.05293e+00-3.16776e-02*t-3.86317e-01/(1.0+exp(t/1.27957e-01)))*0.9:8.628644e-01*0.9",RooArgList(t));
	RooFormulaVar eff_model2("eff_model2","0.002143*(1.0-0.01331*t+0.3757/(1.0+exp(t/-0.1411)))",RooArgList(t));
	RooFormulaVar eff_model3("eff_model3","t<6E-12?(1.00570e+00-1.37914e-02*t)*0.9:9.229534e-01*0.9",RooArgList(t));
    }
    else
    {
	RooFormulaVar eff_model1("eff_model1","t<6E-12?(9.62934e-01+7.67341e-02*t*1E12-3.26842e-02*t*t*1E24+2.73130e-03*t*t*t*1E36)*0.9:8.366682e-01*0.9",RooArgList(t));
	//RooFormulaVar eff_model2("eff_model2","t<6E-12?(1.05293e+00-3.16776e-02*t*1E12-3.86317e-01/(1.0+exp(t*1E12/1.27957e-01)))*0.9:8.628644e-01*0.9",RooArgList(t));
	RooFormulaVar eff_model2("eff_model2","0.002143*(1.0-0.01331*t*1E12+0.3757/(1.0+exp(t*1E12/-0.1411)))",RooArgList(t));
	RooFormulaVar eff_model3("eff_model3","t<6E-12?(1.00570e+00-1.37914e-02*t*1E12)*0.9:9.229534e-01*0.9",RooArgList(t));
    }
    
    RooEffProd pdf_t_signal("pdf_t_signal","pdf_t_signal",pdf_t_signal_raw,(PROD_T_EFFICIENCY==1?eff_model1:(PROD_T_EFFICIENCY==2?eff_model2:eff_model3)));
#else
    RooRealVar tau("tau","tau",1.5*PS, 0.1*PS, 3*PS);
    RooDecay pdf_t_signal("pdf_t_signal","pdf_t_signal",t,tau,res_signal,RooDecay::SingleSided);
#endif
        
    RooProdPdf pdf_signal("pdf_signal","pdf_signal", RooArgSet(pdf_m_signal, pdf_t_signal));
    
    //-----------------------------------------------------------------
    // background (prompt,non-prompt) PDF
    
    RooRealVar m_par1("m_par1","m_par1",0.,-10.,+10.);
    RooChebychev pdf_m_background("pdf_m_background","pdf_m_background",mass,RooArgList(m_par1));
    
#if PER_EVENT_ERROR    
    RooRealVar res_prt_mean("res_prt_mean","res_prt_mean",0.0,-1.0,1.0);
    RooRealVar res_prt_sigma("res_prt_sigma","res_prt_sigma",1.0,0.1,3.0);
    RooGaussModel res_prompt("res_prompt","res_prompt",t,res_prt_mean,res_prt_sigma,tE);
    
    RooRealVar res_npt_mean("res_npt_mean","res_npt_mean",0.0,-2.0,2.0);
    RooRealVar res_npt_sigma("res_npt_sigma","res_npt_sigma",1.0,0.1,5.0);
    RooGaussModel res_nonprompt("res_nonprompt","res_nonprompt",t,res_npt_mean,res_npt_sigma,tE);
#else
    RooRealVar res_prt_mean("res_prt_mean","res_prt_mean",0.0,-0.2*PS,0.2*PS);
    RooRealVar res_prt_sigma1("res_prt_sigma1","res_prt_sigma1",0.1*PS,0.01*PS,0.3*PS);
    RooRealVar res_prt_sigma2("res_prt_sigma2","res_prt_sigma2",0.2*PS,0.01*PS,0.6*PS);
    RooRealVar res_prt_fraction("res_prt_fraction","res_prt_fraction",0.8,0.0,1.0);

    RooGaussModel res_prt_gaussian1("res_prt_gaussian1","res_prt_gaussian1",t,res_prt_mean,res_prt_sigma1);
    RooGaussModel res_prt_gaussian2("res_prt_gaussian2","res_prt_gaussian2",t,res_prt_mean,res_prt_sigma2);
    
    RooAddModel res_prompt("res_prompt","res_prompt",RooArgList(res_prt_gaussian1,res_prt_gaussian2),RooArgList(res_prt_fraction));
    
    RooRealVar res_npt_mean("res_npt_mean","res_npt_mean",0.0,-0.2*PS,0.2*PS);
    RooRealVar res_npt_sigma("res_npt_sigma","res_npt_sigma",0.1*PS,0.001*PS,0.5*PS);
    RooGaussModel res_nonprompt("res_nonprompt","res_nonprompt",t,res_npt_mean,res_npt_sigma);
#endif    
    
    RooRealVar tau_npt("tau_npt","tau_npt",1.0*PS, 0.1*PS, 3.0*PS);
    RooDecay pdf_t_nonprompt("pdf_t_nonprompt","pdf_t_nonprompt",t,tau_npt,res_nonprompt,RooDecay::SingleSided);

    RooProdPdf pdf_prompt("pdf_prompt","pdf_prompt",RooArgSet(pdf_m_background,res_prompt));
    
    RooProdPdf pdf_nonprompt("pdf_nonprompt","pdf_nonprompt",RooArgSet(pdf_m_background, pdf_t_nonprompt));
    
    //-----------------------------------------------------------------
    // full model
    
    RooRealVar n_signal("n_signal","n_signal",1000.,0.,20000.);
    RooRealVar n_prompt("n_prompt","n_prompt",1000.,0.,20000.);
    RooRealVar n_nonprompt("n_nonprompt","n_nonprompt",1000.,0.,20000.);
    
    RooAddPdf model("model","model",
                    RooArgList(pdf_signal, pdf_prompt, pdf_nonprompt),
                    RooArgList(n_signal, n_prompt, n_nonprompt));
    
#if PER_EVENT_ERROR
    
    // Set some of the initial values
    m_fraction.setVal(      4.82199e-01 );
    m_mean.setVal(          5.61965e+00 );
    m_par1.setVal(         -9.16941e-02 );
    m_sigma1.setVal(        7.54107e-03 );
    m_sigma2.setVal(        0.020       );    
    n_nonprompt.setVal(     1.09443e+03 );
    n_prompt.setVal(        2.57430e+03 );
    n_signal.setVal(        9.97147e+02 );
    res_npt_mean.setVal(   -1.16289e+00 );
    res_npt_sigma.setVal(   2.39883e+00 );
    res_prt_mean.setVal(    1.35161e-01 );
    res_prt_sigma.setVal(   9.92769e-01 );
    res_sig_mean.setVal(   -2.21912e-01 );
    res_sig_sigma.setVal(   7.92625e-01 );
    tau.setVal(             1.51266*PS );
    tau_npt.setVal(         1.15003*PS );    
    
#if MATCH_MAIN_ANALYSIS
    m_sigma2.setVal(0.020);
    m_sigma2.setConstant(kTRUE);
#endif
    
    if (do_generation) {
        
        int ntotal = rnd.Poisson(n_signal.getVal()+n_prompt.getVal()+n_nonprompt.getVal());
        
        /* Model of tE
         FCN=-138575 FROM MINOS     STATUS=SUCCESSFUL    766 CALLS        1035 TOTAL
         EDM=0.000147931    STRATEGY= 1      ERROR MATRIX ACCURATE
         EXT PARAMETER                  PARABOLIC         MINOS ERRORS
         NO.   NAME      VALUE            ERROR      NEGATIVE      POSITIVE
         1  fraction     8.75478e-01   2.95202e-02  -3.20404e-02   2.71812e-02
         2  mean1        8.27692e-14   6.97100e-16  -7.21794e-16   6.72467e-16
         3  mean2        1.14172e-13   4.00315e-15  -3.74682e-15   4.31621e-15
         4  sigma1       1.98424e-14   5.25742e-16  -5.33754e-16   5.17362e-16
         5  sigma2       3.57622e-14   2.16486e-15  -1.97930e-15   2.38415e-15
         */
        
        RooRealVar tE_mean1("tE_mean1","tE_mean1",8.27692e-14);
        RooRealVar tE_mean2("tE_mean2","tE_mean2",1.14172e-13);
        RooRealVar tE_sigma1("tE_sigma1","tE_sigma1",1.98424e-14);
        RooRealVar tE_sigma2("tE_sigma2","tE_sigma2",3.57622e-14);
        RooRealVar tE_fraction("tE_fraction","tE_fraction",8.75478e-01);
        
        RooGaussian tE_gaussian1("tE_gaussian1","tE_gaussian1",tE,tE_mean1,tE_sigma1);
        RooGaussian tE_gaussian2("tE_gaussian2","tE_gaussian2",tE,tE_mean2,tE_sigma2);
        
        RooAddPdf tE_model("tE_model","tE_model",RooArgList(tE_gaussian1,tE_gaussian2),RooArgList(tE_fraction));
        RooDataSet* tE_data = tE_model.generate(RooArgSet(tE),ntotal);
        
        data = model.generate(RooArgSet(mass,t),ProtoData(*tE_data));
        data->write(TString::Format("pseudo-%03d.ascii",idx));
        
    } else model.fitTo(*data,ConditionalObservables(tE),Minos(DO_MINOS),NumCPU(NUMBER_OF_CPU));
#else
    
    // Set some of the initial values    
    m_fraction.setVal(       4.95217e-01  );
    m_mean.setVal(           5.61964e+00  );
    m_par1.setVal(          -9.47916e-02  );
    m_sigma1.setVal(         7.60537e-03  );
    m_sigma2.setVal(         0.020        );
    n_nonprompt.setVal(      9.77171e+02  );
    n_prompt.setVal(         2.70706e+03  );
    n_signal.setVal(         9.81786e+02  );
    res_npt_mean.setVal(     0.0          );
    res_npt_sigma.setVal(    3.52986e-14  );
    res_prt_fraction.setVal( 8.93634e-01  );
    res_prt_mean.setVal(     0.0          );
    res_prt_sigma1.setVal(   9.26863e-14  );
    res_prt_sigma2.setVal(   2.57134e-13  );
    res_sig_mean.setVal(     0.0          );
    res_sig_sigma.setVal(    6.74784e-14  );
    tau.setVal(              1.51357*PS  );
    tau_npt.setVal(          1.17838*PS  );
    
#if MATCH_MAIN_ANALYSIS
    m_sigma2.setConstant(kTRUE);
    res_sig_mean.setConstant(kTRUE);
    res_prt_mean.setConstant(kTRUE);
    res_npt_mean.setConstant(kTRUE);
#endif
    res_npt_mean.setConstant(kTRUE);
    
    if (do_generation) {
        
        int ntotal = rnd.Poisson(n_signal.getVal()+n_prompt.getVal()+n_nonprompt.getVal());
        data = model.generate(RooArgSet(mass,t,tE),ntotal);
        data->write(TString::Format("pseudo-%03d.ascii",idx));
        
    } else model.fitTo(*data,Minos(DO_MINOS),NumCPU(NUMBER_OF_CPU));
#endif

    // save lifetime results to log for toy MC fit
    if (idx>=0 && !do_generation) {
        FILE *fp = fopen("results.log","a");
        fprintf(fp,"%d %g %g\n",idx,tau.getVal(),tau.getError());
        fclose(fp);
    }

#if DISPLAY
    TCanvas *c1 = new TCanvas("c1","c1",1000,400);
    c1->Divide(2);
        
    // Display
    c1->cd(1);
    RooPlot* frame_m = mass.frame();
    data->plotOn(frame_m,Binning(120));
    model.plotOn(frame_m);
    model.plotOn(frame_m,Components(pdf_prompt),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_m,Components(pdf_nonprompt),LineColor(kGreen+1),LineWidth(2),LineStyle(7));
    frame_m->SetTitle("");
    frame_m->Draw();

    c1->cd(2);
    RooPlot* frame_t = t.frame();
    data->plotOn(frame_t,Binning(128));
#if PER_EVENT_ERROR
    model.plotOn(frame_t,ProjWData(tE,*data));
    model.plotOn(frame_t,ProjWData(tE,*data),Components(pdf_prompt),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_t,ProjWData(tE,*data),Components(pdf_nonprompt),LineColor(kGreen+1),LineWidth(2),LineStyle(7));
#else
    model.plotOn(frame_t,Binning(128));
    model.plotOn(frame_t,Components(pdf_prompt),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_t,Components(pdf_nonprompt),LineColor(kGreen+1),LineWidth(2),LineStyle(7));
#endif
    frame_t->SetTitle("");
    frame_t->Draw();
    c1->GetPad(2)->SetLogy();

    TCanvas *c2 = new TCanvas("c2","c2",1000,400);
    c2->Divide(2);
    
    mass.setRange("sg",MASS_PEAK-0.045,MASS_PEAK+0.045);
    mass.setRange("sb_lo",MASS_MIN,MASS_PEAK-0.045);
    mass.setRange("sb_hi",MASS_PEAK+0.045,MASS_MAX);
    
    c2->cd(1);
    RooPlot* frame_t1 = t.frame();
    data->plotOn(frame_t1,CutRange("sg"),Binning(128));
#if PER_EVENT_ERROR
    model.plotOn(frame_t1,ProjectionRange("sg"),ProjWData(tE,*data));
    model.plotOn(frame_t1,ProjectionRange("sg"),ProjWData(tE,*data),Components(pdf_prompt),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_t1,ProjectionRange("sg"),ProjWData(tE,*data),Components(pdf_nonprompt),LineColor(kGreen+1),LineWidth(2),LineStyle(7));
#else
    model.plotOn(frame_t1,ProjectionRange("sg"),Binning(128));
    model.plotOn(frame_t1,ProjectionRange("sg"),Components(pdf_prompt),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_t1,ProjectionRange("sg"),Components(pdf_nonprompt),LineColor(kGreen+1),LineWidth(2),LineStyle(7));
#endif
    frame_t1->SetTitle("Signal region");
    frame_t1->Draw();
    c2->GetPad(1)->SetLogy();
    
    c2->cd(2);
    RooPlot* frame_t2 = t.frame();
    data->plotOn(frame_t2,CutRange("sb_lo,sb_hi"),Binning(128));
#if PER_EVENT_ERROR
    model.plotOn(frame_t2,ProjectionRange("sb_lo,sb_hi"),ProjWData(tE,*data));
    model.plotOn(frame_t2,ProjectionRange("sb_lo,sb_hi"),ProjWData(tE,*data),Components(pdf_prompt),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_t2,ProjectionRange("sb_lo,sb_hi"),ProjWData(tE,*data),Components(pdf_nonprompt),LineColor(kGreen+1),LineWidth(2),LineStyle(7));
#else
    model.plotOn(frame_t2,ProjectionRange("sb_lo,sb_hi"),Binning(128));
    model.plotOn(frame_t2,ProjectionRange("sb_lo,sb_hi"),Components(pdf_prompt),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    model.plotOn(frame_t2,ProjectionRange("sb_lo,sb_hi"),Components(pdf_nonprompt),LineColor(kGreen+1),LineWidth(2),LineStyle(7));
#endif
    frame_t2->SetTitle("Sideband");
    frame_t2->Draw();
    c2->GetPad(2)->SetLogy();
#endif    
    
}
