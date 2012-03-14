#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooTruthModel.h"

#include "utils.h"
#include "setTDRStyle_modified.C"


using std::cout;
using std::endl;
using std::string;
using namespace RooFit;

int doCtFitGenReco(TTree *tree, string title, bool isB0, double tLo = -2e-12, double tHi = 16e-12)
{
    setTDRStyle();
    const int nBins(50);

    // Observables in TTree
    RooRealVar t("t", "t", tLo, tHi);
    RooRealVar tE("tE", "tE", 0, 20e-12);

    // Data
    RooDataSet data("data", "dataset", RooArgSet(t,tE), Import(*tree));
    int curEntries = data.numEntries();
    cout << "Dataset contains " << curEntries << " entries." << endl;

    // Fit parameters
    //const double tauLb(1.391e-12), tauB0(1.525e-12);
    const double tauLb(1.229e-12), tauB0(1.536e-12); // corresponds to the evt.pdl table used in production
    const double tauTruthVal = (isB0 ? tauB0 : tauLb);
    RooRealVar tau("tau","tau", tauTruthVal, 0e-12, (tHi < 2*tauTruthVal ? 2*tauTruthVal : tHi));

    RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries*1.4) ;

    // Model building
    // -- decay
    RooRealVar bias("bias","bias",0,-10,10);
    RooRealVar sigma("sigma","per-event error scale factor",1,0.1,10);
    RooGaussModel gm("gm1","gauss model scaled by per-event error",t,bias,sigma,tE) ;
    RooTruthModel idealres("idealres","Ideal resolution model",t);
    RooDecay decay("decay", "decay", t, tau, gm, RooDecay::SingleSided);

    // final model
    RooAddPdf model("model", "model", RooArgSet(decay), RooArgList(nsig));

    // do the fit
    model.fitTo(data);

    // now make a model using the truth lifetime
    RooRealVar tauTruth("tauTruth", "tau truth", tauTruthVal);
    RooDecay decayTruth("decayTruth", "decay truth", t, tauTruth, idealres, RooDecay::SingleSided);
    RooAddPdf modelTruth("modelTruth", "model truth", RooArgSet(decayTruth), RooArgList(nsig));

    // some plotting
    RooPlot * frame_t = t.frame(Title(title.c_str()));
    data.plotOn(frame_t, Binning(nBins));
    modelTruth.plotOn(frame_t,ProjWData(data),LineStyle(kDashed),LineColor(kGreen));
    model.plotOn(frame_t,ProjWData(data));

    // write the results to the canvas
    const double txtPosLeft   = .40;
    const double txtPosTop    = .80;
    const double txtLineSpace = .05;
    const double textsize     = .04;

    frame_t->addObject(writeTLatex("#tau: " + roundToString(tau.getVal()*1e12, 3) + " #pm " + roundToString(tau.getError()*1e12, 3) + " ps", txtPosLeft, txtPosTop-0*txtLineSpace, textsize));
    frame_t->addObject(writeTLatex("n Signal: " + roundToString(nsig.getVal(), 0) + " #pm " + roundToString(nsig.getError(), 0), txtPosLeft, txtPosTop-1*txtLineSpace, textsize));
    frame_t->addObject(writeTLatex("(errors statistical only)", txtPosLeft, txtPosTop-2*txtLineSpace, textsize));

    // draw the plots
    TCanvas *canvas = new TCanvas("fitres","fit results",800,800);
    canvas->Divide(1,1);
    canvas->cd(1);
    frame_t->Draw();
    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.12);
    gPad->SetTopMargin(0.10);

    // make a legend
    const double legPosLeft = .20;
    const double legPosTop = .20;
    const double legWidth = .30;
    const double legHeight = .07;
    TLegend *legend = new TLegend(legPosLeft, legPosTop, legPosLeft+legWidth, legPosLeft+legHeight);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    legend->AddEntry(frame_t->findObject(isB0 ? "model_Norm[t3dB0]" : "model_Norm[t3dLb]"), "fit", "l");
    legend->AddEntry(frame_t->findObject(isB0 ? "modelTruth_Norm[t3dB0]" : "modelTruth_Norm[t3dLb]"), ("truth (" + roundToString(tauTruthVal*1e12,3) + " ps)").c_str(), "l");
    legend->Draw();

    // get out object names (only needed to catch names of objects for legend during development)
    // taken from http://root.cern.ch/phpBB3/viewtopic.php?p=31694
    /*
    cout << "Frame objects:\n";
    for (int i=0; i<frame_t->numItems(); i++) {
	const string obj_name = frame_t->nameOf(i);
	if (obj_name=="") continue;
	cout << i << ".: " << obj_name << endl;
	//cout << Form("%d. '%s'\n",i,obj_name.Data());
    }
    */
//	0.: h_data
//	1.: model_Norm[t3dLb]
//	2.: modelTruth_Norm[t3dLb]


    // we reached the end without an error
    return 0;
}

int doCtFitGenReco(string filename, string title, bool isB0, double tLo = 0e-12, double tHi = 16e-12)
{
    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
	return -1;
    }
    TTree *tree = (TTree*)f->Get("fittree");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
	cout << "This script requires a very reduced tree called fittree" << endl;
	return -2;
    }
    return doCtFitGenReco(tree, title, isB0, tLo, tHi);
}

// Variante mit Cut
int doCtFitGenReco(string filename, string cutstring, string title, bool isB0, double tLo = 0e-12, double tHi = 16e-12)
{
    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
	return -1;
    }
    TTree *tree = (TTree*)f->Get("genevents");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
	return -2;
    }
    // do a preselection
    TTree* treeCutted = 0;
    TFile* tmpfile = new TFile("tmpDoFitGen.root","RECREATE");
    if (tmpfile == 0)
    {
       cout << "temporary file creation failed - exiting" << endl;
       return -3;
    }
    cout << "tmpfile: " << tmpfile << endl;
    cout << "Preselecting events..." << endl;
    cout << "Cutstring: " << cutstring << endl;

    treeCutted = tree->CopyTree(cutstring.c_str());
    cout << "now " << treeCutted->GetEntries() << " of " << tree->GetEntries()<< endl;
    delete tree;
    tmpfile->Write();

    return doCtFitGenReco(treeCutted, title, isB0, tLo, tHi);
}

int doSomeCtFitPlots(string run)
{
    string filename = "../data/" + run + ".root";
    string pdfFilebase = run;
    string titlePref = "B^{0} lifetime gen data";
    string titleFull = titlePref; // default behaviour
    if (run=="run249") titleFull = titlePref + " - private production";
    if (run=="run250") titleFull = titlePref + " - Summer11 new run";
    if (run=="run251") titleFull = titlePref + " - Summer11 run ~month ago";
    doCtFitGenReco(filename, titleFull, true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_nocut.pdf").c_str());

    // Cuts auf B0
    doCtFitGenReco(filename, "ptB0>10", "B^{0} lifetime gen data - p_{T}(B^{0})>10", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_ptB0_gt10.pdf").c_str());

    return 0;

    doCtFitGenReco(filename, "pB0>15", "B^{0} lifetime gen data - p(B^{0})>15", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_pB0_gt15.pdf").c_str());

    doCtFitGenReco(filename, "TMath::Abs(etaB0)<1", "B^{0} lifetime gen data - |#eta(B^{0})|<1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_etaB0_lt1.pdf").c_str());

    doCtFitGenReco(filename, "TMath::Abs(etaB0)>1", "B^{0} lifetime gen data - |#eta(B^{0})|>1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_etaB0_gt1.pdf").c_str());

    // Cuts auf Muonen
    doCtFitGenReco(filename, "ptmu1>7&&ptmu2>7", "B^{0} lifetime gen data - p_{T}(#mu_{1,2})>7", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_ptmu1_gt7.pdf").c_str());

    doCtFitGenReco(filename, "ptmu1>3&&ptmu2>3", "B^{0} lifetime gen data - p_{T}(#mu_{1,2})>3", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_ptmu1_gt3.pdf").c_str());

    doCtFitGenReco(filename, "pmu1>4&&pmu2>4", "B^{0} lifetime gen data - p(#mu_{1,2})>4", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_pmu1_gt4.pdf").c_str());

    doCtFitGenReco(filename, "TMath::Abs(etamu1)<1&&TMath::Abs(etamu2)<1", "B^{0} lifetime gen data - |#eta(#mu_{1,2})|<1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_etamu1_lt1.pdf").c_str());

    doCtFitGenReco(filename, "TMath::Abs(etamu1)>1&&TMath::Abs(etamu2)>1", "B^{0} lifetime gen data - |#eta(#mu_{1,2})|>1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_etamu1_gt1.pdf").c_str());

    return 0;
}

