#include <iostream>
#include <fstream>

#include "TTree.h"
#include "TCanvas.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include <string>
#include "TLatex.h"
#include "TText.h"
#include <sstream>
#include <iomanip>

using namespace std;
using namespace RooFit;
TCanvas *canvas1, *canvas2;

template <typename T>
std::string toString(T i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}

std::string roundToString(double v, std::streamsize precision)
{   
    ostringstream oss;
    oss << std::setprecision(precision) << std::fixed << v;
    return oss.str();
}

void MassFit(){
  TLatex* txt;
  canvas1= new TCanvas("mass","Lambda_b mass");
  cout <<"Canvas opened"<<endl;
  string tmp;
  int a;
  double Lb_mass, L0_lifetime, Lb_lifetime;
  //std::string strTitle = "MuOnia L #approx40nb^{-1}";
  //ifstream infile("mlb1011081859.out");
  //ifstream infile("mlb1011081859cut25s.out");
  std::string strTitle = "MC #Lambda_{b}";
  ifstream infile("mlbMCsig1012091538.out");
  TTree *tree=new TTree("LBmass","Mass of Lambda_b");;
  tree->Branch("mass", &Lb_mass, "mass/D");
  tree->Branch("l0_lifetime", &L0_lifetime, "l0_lifetime/D");
  tree->Branch("lb_lifetime", &Lb_lifetime, "lb_lifetime/D");

  while(!infile.eof()){
    infile >> a >> Lb_mass >> L0_lifetime >> Lb_lifetime;
    Lb_lifetime*=10.0;
    tree->Fill();
  }
  infile.close();

  RooRealVar mass("mass", "J/#psi #Lambda mass [GeV]", 5.22, 6.045);
  RooDataSet data("data", "Lambda_b dataset", tree, mass);

  RooRealVar lb_lifetime("lb_lifetime", "c#tau #Lambda_{b} [mm]",0.25, 3.5);
  RooDataSet lt_data("lt_data", "#Lambda_{b} lifetime dataset", tree, lb_lifetime);

  //   mass fit

  RooRealVar mean("mean","mean",5.62,5.6,5.7) ;
  RooRealVar sigma("sigma","sigma",0.02,0.01,0.04) ;
  RooGaussian sig("sig","signal p.d.f.",mass,mean,sigma) ;

  RooRealVar c0("c0","coefficient #0", 1.0,-15,15) ;
  RooRealVar c1("c1","coefficient #1", -0.3,-15,15) ;
  RooPolynomial bkg("bkg","background p.d.f.",mass,RooArgList(c0,c1),0) ;

  RooRealVar nsig("nsig","signal fraction", 55.,0.,250.) ;
  RooRealVar nbkg("nbkg","Background fraction", 520.,0.,1050.) ;

  RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg)) ;
  model.fitTo(data,Extended(kTRUE));

  RooPlot* xframe = mass.frame(Name("xframe"),Title(("Mass of #Lambda_{b} - "+strTitle).c_str())) ;
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

  {
    const double txtPosLeft = 5.7;
    //const double txtPosTop = 50;
    //const double txtLineSpace = 5;
    const double txtPosTop = 100;
    const double txtLineSpace = 20;
    txt = new TLatex(txtPosLeft,txtPosTop-0*txtLineSpace,("n Signal: " + roundToString(nsigWindow,0)).c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    txt = new TLatex(txtPosLeft,txtPosTop-1*txtLineSpace,("n Bgr: " + roundToString(nbkgWindow,1)).c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    txt = new TLatex(txtPosLeft,txtPosTop-2*txtLineSpace,("S/#sqrt{S+B}: " + roundToString(nsigWindow/sqrt(nsigWindow+nbkgWindow),0)).c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    txt = new TLatex(txtPosLeft,txtPosTop-3*txtLineSpace,("mass: " + roundToString(m,3) + " GeV/c^{2}").c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    txt = new TLatex(txtPosLeft,txtPosTop-4*txtLineSpace,("width: " + roundToString(s,3) + " GeV/c^{2}").c_str());
    txt->SetTextSize(0.04);
    xframe->addObject(txt);
    xframe->Draw();
  }


  //  lifetime fit

  RooRealVar tau("tau","tau",-2.3,-5.,-0.5) ;
  RooExponential lt_model("lt_model","Lifetime model",lb_lifetime,tau);
  lt_model.fitTo(lt_data);

  canvas2= new TCanvas("lblt","Lambda_b lifetime");
  RooPlot* ltframe = lb_lifetime.frame(Name("ltframe"),Title(("Lifetime of #Lambda_{b} - "+strTitle).c_str()));
  lt_data.plotOn(ltframe,Binning(20)) ;
  lt_model.plotOn(ltframe);
  ltframe->Draw();

  cout << "c tau = "<<-1/tau.getVal()<<" mm"<<endl;

  txt = new TLatex(2,20,("c#tau: " + roundToString(1000*-1/tau.getVal(),0) + " #mum").c_str());
  txt->SetTextSize(0.04);
  ltframe->addObject(txt);
  txt = new TLatex(2,18.5,"(Very preliminary!)");
  txt->SetTextSize(0.04);
  ltframe->addObject(txt);
  ltframe->Draw();
}


