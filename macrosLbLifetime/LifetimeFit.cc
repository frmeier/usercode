#include <iostream>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include <string>
#include "TLatex.h"
#include "TText.h"
#include <sstream>
#include <iomanip>

#include "cut.h"
#include "cuts.C"

using namespace std;
using namespace RooFit;
TCanvas *canvas;


std::string roundToString(double v, std::streamsize precision)
{   
    ostringstream oss;
    oss << std::setprecision(precision) << std::fixed << v;
    return oss.str();
}

void LifetimeFit(){
  TLatex* txt;
  canvas= new TCanvas("mass","Lambda_b mass",800,800);
  cout <<"Canvas opened"<<endl;

  std::string strTitle = "MC #Lambda_{b}";
  //TFile f("/scratch/frmeier/run129.root");
  //TFile f("/scratch/frmeier/run127.root");
  // TFile f("/scratch/frmeier/run152.root");
  TFile f("cl0099.chain.root");
  TTree *tree=(TTree*)f.Get("events");

  Cuts cut;
  //  cut.selectCut("acc03","an05","HLT_matched_01");
  cut.selectCut("acc03","an05");

  TTree *subtree, *lttree;
  subtree = tree->CopyTree(cut.getCut().c_str());
  cout << "subtree: " << subtree->GetEntries() << " entries of " << tree->GetEntries()<< endl;
  
  RooRealVar mass("mlb", "J/#psi #Lambda mass [GeV]", 5.22, 6.045);
  RooDataSet data("data", "Lambda_b dataset", subtree, mass);

  //   mass fit

  RooRealVar mean("mean","mean",5.62,5.6,5.7) ;
  RooRealVar sigma("sigma","sigma",0.02,0.01,0.04) ;
  RooGaussian sig("sig","signal p.d.f.",mass,mean,sigma) ;

  RooRealVar c0("c0","coefficient #0", 1.0,-1.0,2.0) ;
  RooChebychev bkg("bkg","background p.d.f.",mass,RooArgSet(c0)) ;

  RooRealVar nsig("nsig","signal fraction", 55.,0.,250.) ;
  RooRealVar nbkg("nbkg","Background fraction", 520.,0.,1050.) ;

  RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg)) ;
  model.fitTo(data,Extended(kTRUE));

  RooPlot* xframe = mass.frame(Name("xframe"),Title("#Lambda_{b} mass - data")) ;
  xframe->SetTitle("");
  data.plotOn(xframe,Binning(33)) ;
  model.plotOn(xframe);
  model.plotOn(xframe,Components("bkg"),LineStyle(kDashed));

  xframe->Draw();

  Double_t m=mean.getVal();
  Double_t s=sigma.getVal();

  mass.setRange("window",m-2.5*s,m+2.5*s);
  RooAbsReal* fracSigRange = sig.createIntegral(mass,mass,"window");
  Double_t nsigWindow = fracSigRange->getVal()*nsig.getVal();
  RooAbsReal* fracBGRange = bkg.createIntegral(mass,mass,"window");
  Double_t nbkgWindow = nbkg.getVal()*fracBGRange->getVal();

  cout << "n_Signal     =  "<<nsigWindow<<endl;
  cout << "n_Background =  "<<nbkgWindow<<endl;
  cout << "S/Sqrt(S+B)  =  "<<nsigWindow/sqrt(nsigWindow+nbkgWindow)<<endl;


    const double txtPosLeft = 5.69;
    //const double txtPosTop = 50;
    //const double txtLineSpace = 5;
    const double txtPosTop = 100/4.;
    const double txtLineSpace = 10/4.;
    /*
    txt = new TLatex(txtPosLeft,txtPosTop-0*txtLineSpace,("n Signal: " + roundToString(nsigWindow,0)).c_str());
    txt->SetNDC(true);
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    txt = new TLatex(txtPosLeft,txtPosTop-1*txtLineSpace,("n Bgr: " + roundToString(nbkgWindow,1)).c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    txt = new TLatex(txtPosLeft,txtPosTop-2*txtLineSpace,("S/#sqrt{S+B}: " + roundToString(nsigWindow/sqrt(nsigWindow+nbkgWindow),0)).c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    */
    txt = new TLatex(txtPosLeft,txtPosTop-3*txtLineSpace,("mass: " + roundToString(m,3) + " GeV/c^{2}").c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    txt = new TLatex(txtPosLeft,txtPosTop-4*txtLineSpace,("width: " + roundToString(s,3) + " GeV/c^{2}").c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    xframe->Draw();
    canvas->SaveAs("Mass_Data.pdf");

  //  lifetime fit
  Cuts lb_window;
  lb_window.cs.addCut(new cutSymWindow("mlb",m));
  lb_window.parvec.push_back(2.5*s);
  lttree = subtree->CopyTree(lb_window.getCut().c_str());
  cout << "lttree: " << lttree->GetEntries() << " entries of " << tree->GetEntries()<< endl;
  lttree->Scan("nRef1G:nRef2G:run:event");
  RooRealVar lb_lifetime("ctlb", "c#tau #Lambda_{b} [mm]",0.02, 0.30);
  RooDataSet lt_data("lt_data", "#Lambda_{b} lifetime dataset", lttree, lb_lifetime);
 
  RooRealVar tau("tau","tau",-27.0,-40.,-2.0) ;
  RooExponential lt_model("lt_model","Lifetime model",lb_lifetime,tau);
  lt_model.fitTo(lt_data);

  RooPlot* ltframe = lb_lifetime.frame(Name("ltframe"),Title("#Lambda_{b} lifetime - data"));
  ltframe->SetTitle("");
  lt_data.plotOn(ltframe,Binning(20)) ;
  lt_model.plotOn(ltframe);
  ltframe->Draw();

  cout << "c tau = "<<-10.0/tau.getVal()<<" mm"<<endl;

  txt = new TLatex(0.5,0.7,("c#tau: " + roundToString(-10000.0/tau.getVal(),0) + "#pm"+roundToString(tau.getError()/tau.getVal()*10000.0/tau.getVal(),0)+ " #mum").c_str());
  txt->SetNDC(true);
  txt->SetTextSize(0.06);
  ltframe->addObject(txt);
  txt = new TLatex(0.5,0.6,("Entries: " + roundToString(lttree->GetEntries(),0)).c_str());
  txt->SetNDC(true);
  txt->SetTextSize(0.06);
  ltframe->addObject(txt);
  ltframe->Draw();
  canvas->SaveAs("Lifetime_Data.pdf");
}


