#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFractionFitter.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "THStack.h"

#include "Cuts.C"

using std::cout;
using std::endl;
using std::string;

void doShapeFit1d()
{
    const string path = "../data/";
    TFile *f_data = TFile::Open((path+"run393__run397.root").c_str());
    //TFile *f_sig = TFile::Open((path+"run398__run399.root").c_str());
    TFile *f_sig = TFile::Open((path+"run404.root").c_str());
    TFile *f_bgr_Jp = TFile::Open((path+"run400.root").c_str());
    TFile *f_bgr_Bs = TFile::Open((path+"run401.root").c_str());
    TFile *f_bgr_Bp = TFile::Open((path+"run402.root").c_str());
    TFile *f_bgr_Lb = TFile::Open((path+"run403.root").c_str());

    TTree *t_data = (TTree*)f_data->Get("events");
    TTree *t_sig = (TTree*)f_sig->Get("events");
    TTree *t_bgr_Jp = (TTree*)f_bgr_Jp->Get("events");
    TTree *t_bgr_Bs = (TTree*)f_bgr_Bs->Get("events");
    TTree *t_bgr_Bp = (TTree*)f_bgr_Bp->Get("events");
    TTree *t_bgr_Lb = (TTree*)f_bgr_Lb->Get("events");

    Cuts cut;
    cut.selectCut("B005", "acc05B0", "muSoft");
    string curCut = cut.getCut();

    TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 1200);
    canvas->Divide(2,2);

    const string binstring = "(50,4.5,6)";
    canvas->cd(1);
    t_data->Draw(("mB0>>h_data"+binstring).c_str(), curCut.c_str());
    t_sig->Draw(("mB0>>h_sig"+binstring).c_str(), curCut.c_str(), "same");
    t_bgr_Jp->Draw(("mB0>>h_bgr_Jp"+binstring).c_str(), curCut.c_str(), "same");
    t_bgr_Bs->Draw(("mB0>>h_bgr_Bs"+binstring).c_str(), curCut.c_str(), "same");
    t_bgr_Bp->Draw(("mB0>>h_bgr_Bp"+binstring).c_str(), curCut.c_str(), "same");
    t_bgr_Lb->Draw(("mB0>>h_bgr_Lb"+binstring).c_str(), curCut.c_str(), "same");

    TH1F *h_data = (TH1F*)gDirectory->GetList()->FindObject("h_data");
    TH1F *h_sig = (TH1F*)gDirectory->GetList()->FindObject("h_sig");
    TH1F *h_bgr_Jp = (TH1F*)gDirectory->GetList()->FindObject("h_bgr_Jp");
    TH1F *h_bgr_Bs = (TH1F*)gDirectory->GetList()->FindObject("h_bgr_Bs");
    TH1F *h_bgr_Bp = (TH1F*)gDirectory->GetList()->FindObject("h_bgr_Bp");
    TH1F *h_bgr_Lb = (TH1F*)gDirectory->GetList()->FindObject("h_bgr_Lb");

    const int nmc = 5;
    TObjArray *mc = new TObjArray(nmc);
    mc->Add(h_sig->Clone("h_sig_fit"));
    mc->Add(h_bgr_Jp->Clone("h_bgr_Jp_fit"));
    mc->Add(h_bgr_Bs->Clone("h_bgr_Bs_fit"));
    mc->Add(h_bgr_Bp->Clone("h_bgr_Bp_fit"));
    mc->Add(h_bgr_Lb->Clone("h_bgr_Lb_fit"));

    TFractionFitter* fit = new TFractionFitter(h_data, mc);

    double lo = 0, hi = 100;
    fit->Constrain(0, lo, hi);
    fit->Constrain(1, lo, hi);
    fit->Constrain(2, lo, hi);
    fit->Constrain(3, lo, hi);
    fit->Constrain(4, lo, hi);

    Int_t status = fit->Fit();
    cout << "fit status: " << status << endl;
    Color_t col[nmc] = { 2, 3, 4, 6, 7};
    if (status == 0)
    {                       // check on fit status
	TH1F* result = (TH1F*) fit->GetPlot();
	canvas->cd(2);
	h_data->Draw("Ep");
	result->Draw("same");

	canvas->cd(3);
	h_data->Draw();
	double max = h_data->GetMaximum();

	canvas->cd(4);
	double w[nmc], wE[nmc];
	THStack *hstack = new THStack("hstack","Stacked");
	for (int i = 0; i!=nmc; i++)
	{
	    fit->GetResult(i, w[i], wE[i]);
	    cout << "res " << i << ": " << w[i] << endl;
	    //fit->GetMCPrediction(i)->Draw(i!=0 ? "same" : "");
	    fit->GetMCPrediction(i)->SetFillColor(col[i]);
	    hstack->Add(fit->GetMCPrediction(i));
	    /*
	    ((TH1F*)mc->At(i))->Scale(7*w[i]);    
	    ((TH1F*)mc->At(i))->SetMinimum(0);
	    ((TH1F*)mc->At(i))->SetMaximum(max);
	    ((TH1F*)mc->At(i))->SetLineColor(col[i]);
	    ((TH1F*)mc->At(i))->Draw(i!=0 ? "same" : "");
	    */
	}
	hstack->Draw();
    }
}

