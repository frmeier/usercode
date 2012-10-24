#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"

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

#include "Cuts.C"

#include "utils.h"
#include "setTDRStyle_modified.C"

bool isRecoTree(true);
bool useEventsTree(false);

using std::cout;
using std::endl;
using std::string;
using namespace RooFit;

int doCtFitGen(TTree *tree, string title, bool isB0, double tLo = 0e-12, double tHi = 16e-12)
{
    setTDRStyle();
    const int nBins(150);

    // Observables in TTree
    //RooRealVar t((isB0 ? "t3dB0" : "t3dLb"), "t", tLo, tHi);
    string prtcl = "bc";
    string treevar;
    if (isRecoTree)
    {	// it is from b0Reader
	if (useEventsTree)
	    treevar = "ct" + prtcl + "truth";
	else
	    treevar = "ct" + prtcl;
    }
    else
    {	// it is from b0GenReader
	treevar = "t3d" + prtcl;
    }
    cout << "Using leaf " << treevar << endl;
    RooRealVar t(treevar.c_str(), "t", tLo, tHi);

    // Data
    RooDataSet data("data", "dataset", RooArgSet(t), Import(*tree));
    int curEntries = data.numEntries();
    cout << "Dataset contains " << curEntries << " entries." << endl;

    // Fit parameters
    //const double tauLb(1.391e-12), tauB0(1.525e-12);
    const double tauLb(1.424e-12), tauB0(1.536e-12);
    //const double tauLb(1.229e-12), tauB0(1.536e-12); // corresponds to the evt.pdl table used in production
    const double tauTruthVal = (isB0 ? tauB0 : tauLb);
    RooRealVar tau("tau","tau", tauTruthVal, 0e-12, (tHi < 2*tauTruthVal ? 2*tauTruthVal : tHi));

    RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries*1.4) ;

    // Model building
    // -- decay
    RooTruthModel idealres("idealres","Ideal resolution model",t);
    RooDecay decay("decay", "decay", t, tau, idealres, RooDecay::SingleSided);

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
    cout << "Result: " << roundToString(tau.getVal()*1e12, 3) << " +/- " << roundToString(tau.getError()*1e12, 3) << " ps  N: " << roundToString(nsig.getVal(), 0) << endl;

    gStyle->SetTitleColor(0);

    // draw the plots
    //TCanvas *canvas = new TCanvas("fitres","fit results",1600,800);
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
    if (useEventsTree)
    {
	legend->AddEntry(frame_t->findObject("model_Norm[ctbctruth]"), "fit", "l");
	legend->AddEntry(frame_t->findObject("modelTruth_Norm[ctbctruth]"), ("truth (" + roundToString(tauTruthVal*1e12,3) + " ps)").c_str(), "l");
    }
    else
    {
	legend->AddEntry(frame_t->findObject("model_Norm[ctbc]"), "fit", "l");
	legend->AddEntry(frame_t->findObject("modelTruth_Norm[ctbc]"), ("truth (" + roundToString(tauTruthVal*1e12,3) + " ps)").c_str(), "l");
    }
    legend->SetBorderSize(0);
    legend->Draw();

    // get out object names (only needed to catch names of objects for legend during development)
    // taken from http://root.cern.ch/phpBB3/viewtopic.php?p=31694
    cout << "Frame objects:\n";
    for (int i=0; i<frame_t->numItems(); i++) {
	const string obj_name = frame_t->nameOf(i);
	if (obj_name=="") continue;
	cout << i << ".: " << obj_name << endl;
	//cout << Form("%d. '%s'\n",i,obj_name.Data());
    }
//	0.: h_data
//	1.: model_Norm[t3dLb]
//	2.: modelTruth_Norm[t3dLb]


    // we reached the end without an error
    return 0;
}

int doCtFitGen(string filename, string title, bool isB0, double tLo = 0e-12, double tHi = 16e-12)
{
    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
	return -1;
    }
    TTree *tree = (TTree*)f->Get(useEventsTree ? "events" : "genevents");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
	return -2;
    }
    return doCtFitGen(tree, title, isB0, tLo, tHi);
    tree->Delete();
    delete tree;
    f->Delete();
    delete f;
}

// Variante mit Cut
int doCtFitGen(string filename, string cutstring, string title, bool isB0, double tLo = 0e-12, double tHi = 16e-12)
{
    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
	return -1;
    }
    TTree *tree = (TTree*)f->Get(useEventsTree ? "events" : "genevents");
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

    return doCtFitGen(treeCutted, title, isB0, tLo, tHi);
}

void setUseEventsTree(bool b) { useEventsTree = b; }

int doSomeCtFitPlots(string run)
{
    string filename = "../data/" + run + ".root";
    string pdfFilebase = run;
    string titlePref = "B^{0} lifetime gen data";
    string titleFull = titlePref; // default behaviour
    if (run=="run249") titleFull = titlePref + " - private production";
    if (run=="run250") titleFull = titlePref + " - Summer11 new run";
    if (run=="run251") titleFull = titlePref + " - Summer11 run ~month ago";
    doCtFitGen(filename, titleFull, true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_nocut.pdf").c_str());

    // Cuts auf B0
    doCtFitGen(filename, "ptB0>10", "B^{0} lifetime gen data - p_{T}(B^{0})>10", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_ptB0_gt10.pdf").c_str());

    return 0;

    doCtFitGen(filename, "pB0>15", "B^{0} lifetime gen data - p(B^{0})>15", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_pB0_gt15.pdf").c_str());

    doCtFitGen(filename, "TMath::Abs(etaB0)<1", "B^{0} lifetime gen data - |#eta(B^{0})|<1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_etaB0_lt1.pdf").c_str());

    doCtFitGen(filename, "TMath::Abs(etaB0)>1", "B^{0} lifetime gen data - |#eta(B^{0})|>1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_etaB0_gt1.pdf").c_str());

    // Cuts auf Muonen
    doCtFitGen(filename, "ptmu1>7&&ptmu2>7", "B^{0} lifetime gen data - p_{T}(#mu_{1,2})>7", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_ptmu1_gt7.pdf").c_str());

    doCtFitGen(filename, "ptmu1>3&&ptmu2>3", "B^{0} lifetime gen data - p_{T}(#mu_{1,2})>3", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_ptmu1_gt3.pdf").c_str());

    doCtFitGen(filename, "pmu1>4&&pmu2>4", "B^{0} lifetime gen data - p(#mu_{1,2})>4", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_pmu1_gt4.pdf").c_str());

    doCtFitGen(filename, "TMath::Abs(etamu1)<1&&TMath::Abs(etamu2)<1", "B^{0} lifetime gen data - |#eta(#mu_{1,2})|<1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_etamu1_lt1.pdf").c_str());

    doCtFitGen(filename, "TMath::Abs(etamu1)>1&&TMath::Abs(etamu2)>1", "B^{0} lifetime gen data - |#eta(#mu_{1,2})|>1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_etamu1_gt1.pdf").c_str());

    return 0;
}

void doCtFitGenCutSequenceB0_bins(double ptpp)
{
    const string filename = "../data/run460_472.root";
    setUseEventsTree(true);

    Cuts cutB0;
    cutB0.selectCut("B006", "acc06B0", "muSoft", "HLT_jpsiBarrel", "HLT_matched");
    cutB0.removeOneCut("rptha1");
    cutB0.removeOneCut("rptha2");
    const string anaAll = "isMCmatch==1&&isSig==1";
    const string anaPass = anaAll+"&&"+cutB0.getCut();

	    //doCtFitGen(filename, anaAll+"&&TMath::Abs(retamu1)<1.25&&TMath::Abs(retamu2)<1.25&&TMath::Abs(retaha1)<1.25&&TMath::Abs(retaha2)<1.25", "B^{0} reco all matched sig, |#eta_{i}|<2.5", true, 0e-12, 16e-12);
	    //return;
    int i = 0;
//    for (double ptpp = 0.9; ptpp <= 2.5; ptpp+=.1)
	for (double ptpm = 0.5; ptpm <= 2.5; ptpm+=.1)
	{
	    const string addCut = "&&(rqha1>0?rptha1:rptha2)>" + toString(ptpp) + "&&(rqha1<0?rptha1:rptha2)>" + toString(ptpm);
	    cout << "addCut: " << addCut << endl;
	    cout << "ptcut: pi+: " << ptpp << " pt-: " << ptpm << endl;
	    doCtFitGen(filename, anaPass+addCut, "B^{0} reco all matched sig, p_{T}(#pi^{+})>"+toString(ptpp)+" p_{T}(#pi^{-})>"+toString(ptpm), true, 0e-12, 16e-12);
	    gPad->SaveAs(("doCtFitGen_B0_ptbins" + toString(i) +".pdf").c_str());
	}
}

void doCtFitGenCutSequenceB0(const string filename)
{
    //TCanvas *c = new TCanvas("c", "c", 400, 400);
    //const string filename = "../data/run365.root";
    //const string filename = "../data/run337.root";
    //TFile *f = TFile::Open("../data/run337.root");
    //TTree *genevents = (TTree*)f->Get("genevents");
    //TTree *events = (TTree*)f->Get("events");
    const string pdfFilebase = "doCtFitGenCutSequenceB0";

    // ------------------------------------------------------------------------- step 0: all
    setUseEventsTree(false);
    doCtFitGen(filename, "", "B^{0} all", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_all.pdf").c_str());

    //doCtFitGen(filename, "hasCand>0", "B^{0} cand", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_cand.pdf").c_str());

    /*
    doCtFitGen(filename, "vrrs>2", "B^{0} cand vrrs>2", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candeffstep2.pdf").c_str());

    doCtFitGen(filename, "vrrs>2&&(vrrs>6?getrn()<.3:1)>0", "B^{0} cand effstep vrrs>2", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candeffstep3.pdf").c_str());

    doCtFitGen(filename, "(vrrs>6?getrn()<.3:1)>0", "B^{0} cand effstep", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candeffstep4.pdf").c_str());
    */

    //doCtFitGen(filename, "getProb(vrrs)<getrn()", "B^{0} cand effenvelope", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_candeffenvelope.pdf").c_str());

    //doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3", "B^{0} simple acc", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_candacc.pdf").c_str());

    //doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3&&(vrrs>6?getrn()<.3:1)>0", "B^{0} simple acc effstep", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_candacceffstep.pdf").c_str());

    //doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3&&(vrrs>6?getrn()<.1:1)>0", "B^{0} simple acc effstep .1", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_candacceffstep2.pdf").c_str());

    //doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3&&getProb(vrrs)<getrn()", "B^{0} simple acc effenvel", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_candacceffenvel.pdf").c_str());

    /*
    doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3&&(ptha1>ptha2?ptha1:ptha2)>1.8&&(ptha1<ptha2?ptha1:ptha2)>0.5&&TMath::Abs(etaha1)<1.3&&TMath::Abs(etaha2)<1.3", "B^{0} acc hacut", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candaccha.pdf").c_str());

    doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3&&(ptha1>ptha2?ptha1:ptha2)>1.8&&(ptha1<ptha2?ptha1:ptha2)>0.5&&TMath::Abs(etaha1)<1.3&&TMath::Abs(etaha2)<1.3&&(vrrs>6?getrn()<.3:1)>0", "B^{0} acc hacut effstep .3", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candacchastep1.pdf").c_str());

    doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3&&(ptha1>ptha2?ptha1:ptha2)>1.8&&(ptha1<ptha2?ptha1:ptha2)>0.5&&TMath::Abs(etaha1)<1.3&&TMath::Abs(etaha2)<1.3&&(vrrs>6?getrn()<.1:1)>0", "B^{0} acc hacut effstep .1", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candacchastep2.pdf").c_str());

    doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3&&(ptha1>ptha2?ptha1:ptha2)>1.8&&(ptha1<ptha2?ptha1:ptha2)>0.5&&TMath::Abs(etaha1)<1.3&&TMath::Abs(etaha2)<1.3&&getProb(vrrs)<getrn()", "B^{0} acc hacut effenvel", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candacchaenvel.pdf").c_str());
    */

    //doCtFitGen(filename, "ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<1.3&&TMath::Abs(etamu2)<1.3&&(ptha1>ptha2?ptha1:ptha2)>1.8&&(ptha1<ptha2?ptha1:ptha2)>0.5&&TMath::Abs(etaha1)<1.3&&TMath::Abs(etaha2)<1.3&&getProbZ(vzrs)<getrn()", "B^{0} acc hacut effenvelz", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_candacchaenvelz.pdf").c_str());


    doCtFitGen(filename, "hasCand>0&&(isMCmatch&1)==1", "B^{0} cand match", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candmatch.pdf").c_str());

    //doCtFitGen(filename, "hasCand>0&&(isMCmatch&1)==0", "B^{0} cand NOmatch", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_candnotmatch.pdf").c_str());

    /*
    // ------------------------------------------------------------------------- step 1: acceptance
    // Muon acceptance
    const string acc_mu1_1 = "TMath::Abs(etamu1)<1.2&&ptmu1>3.5";
    const string acc_mu1_2 = "TMath::Abs(etamu1)>=1.2&&TMath::Abs(etamu1)<1.6&&ptmu1>8-3.75*TMath::Abs(etamu1)";
    const string acc_mu1_3 = "TMath::Abs(etamu1)>=1.6&&TMath::Abs(etamu1)<2.4&&ptmu1>2.0";
    const string acc_mu2_1 = "TMath::Abs(etamu2)<1.2&&ptmu2>3.5";
    const string acc_mu2_2 = "TMath::Abs(etamu2)>=1.2&&TMath::Abs(etamu2)<1.6&&ptmu2>8-3.75*TMath::Abs(etamu2)";
    const string acc_mu2_3 = "TMath::Abs(etamu2)>=1.6&&TMath::Abs(etamu2)<2.4&&ptmu2>2.0";
    const string acc_mu = "(("+acc_mu1_1+")||("+acc_mu1_2+")||("+acc_mu1_3+"))&&(("+acc_mu2_1+")||("+acc_mu2_2+")||("+acc_mu2_3+"))";

    doCtFitGen(filename, acc_mu, "B^{0} Acceptance", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_01_muAcc.pdf").c_str());

    // ------------------------------------------------------------------------- step 2: Jpsi Kandidat
    const string jpCandAll = acc_mu;
    const string jpCandPass = "hasJpCand==1";

    doCtFitGen(filename, jpCandAll+"&&"+jpCandPass, "+J/#psi candidate", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_05_muKsAccTrgJpCand.pdf").c_str());

    //doCtFitGen(filename, jpCandPass, "B^{0} JpCand", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_05_isol_JpCand.pdf").c_str());

    // ------------------------------------------------------------------------- step 3: Ks Kandidat
    const string KsCandAll = jpCandAll + "&&" + jpCandPass;
    const string KsCandPass = "hasKsCand==1";

    doCtFitGen(filename, KsCandAll+"&&"+KsCandPass, "B^{0} +K_{s} candidate", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_06_muKsAccTrgJpCandKsCand.pdf").c_str());

    //doCtFitGen(filename, KsCandPass, "B^{0} KsCand", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_06_isol_KsCand.pdf").c_str());

    // ------------------------------------------------------------------------- step 4: B0 Kandidat
    const string lbCandAll = KsCandAll + "&&" + KsCandPass;
    const string lbCandPass = "hasCand==1";

    doCtFitGen(filename, lbCandAll+"&&"+lbCandPass, "B^{0} + B^{0} candidate unmatched", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_07_muKsAccTrgJpCandKsCandB0Cand.pdf").c_str());

    //doCtFitGen(filename, lbCandPass, "B^{0} B0Cand", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_07_isolB0Cand.pdf").c_str());

    doCtFitGen(filename, lbCandAll+"&&"+lbCandPass+"&&isMCmatch==1", "B^{0} B^{0} candidate MCmatched", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_07_isolB0CandMatch.pdf").c_str());

    // ------------------------------------------------------------------------- step 5: Trigger
    const string trigAll = lbCandAll + "&&" + lbCandPass;
    const string trigPass = "HLTokBarrelJpsi==1";
    const string trigPassAll = trigAll + "&&" + trigPass;

    doCtFitGen(filename, trigPassAll, "B^{0} + barrel trigger", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_04_muKsAccTrg.pdf").c_str());

    // ------------------------------------------------------------------------- step 6: Triggermatching
    const string hltMatchAll = lbCandAll + "&&" + lbCandPass;
    const string hltMatchPass = "HLTmatch==1";

    doCtFitGen(filename, hltMatchAll+"&&"+hltMatchPass, "B^{0} muKs Acc Trg JpCand KsCand B0Cand TrgMatch", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_08_muKsAccTrgJpCandKsCandB0CandHltmatch.pdf").c_str());

    // everything up to now:
    const string hltMatchPassAll = hltMatchAll + "&&" + hltMatchPass;

    */

    // =====================================================================================================
    // ------------------------------------------------------------------------- step 7: analysis cuts
    setUseEventsTree(true);

    Cuts cutB0;
    cutB0.selectCut("B008", "acc06B0", "muSoft");
    const string anaAll = "isMCmatch==1&&isSig==1";
    const string anaPass = anaAll+"&&"+cutB0.getCut();
    const string anaPassBarrel = anaPass+"&&"+"HLTmatch==1&&HLTokBarrelJpsi==1";
    const string anaPassDispl = anaPass+"&&"+"HLTmatch==1&&HLTokDisplJpsi==1";

    //doCtFitGen(filename, "", "B^{0} reco all", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_09_recoAll.pdf").c_str());

    /*
    for (int i=0; i!=5; i++)
    {
	doCtFitGen(filename, "vrrs>"+toString(i), "B^{0} reco all, vrrs>"+toString(i), true, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_09_recoAll_rvrrs_"+toString(i)+".pdf").c_str());
    }
    */

    //doCtFitGen(filename, "isMCmatch==1", "B^{0} reco all matched", true, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_09_recoAllMatch.pdf").c_str());

    /*
    doCtFitGen(filename, anaAll, "B^{0} reco all matched sig", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_09_recoAllMatchSig.pdf").c_str());

    */

    doCtFitGen(filename, anaPass, "B^{0} analysis cuts", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacut.pdf").c_str());

    doCtFitGen(filename, anaPassBarrel, "B^{0} analysis cuts, barrel trigger", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutBarrel.pdf").c_str());

    for (int i=0; i!=5; i++)
    {
	doCtFitGen(filename, anaPassBarrel+"&&vrrs>"+toString(i), "B^{0} analysis cuts, barrel trigger, vrrs>"+toString(i), true, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_10_nacutBarrel_rvrrs_"+toString(i)+".pdf").c_str());
    }

    doCtFitGen(filename, anaPassDispl, "B^{0} analysis cuts, displaced trigger", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutDispl.pdf").c_str());


    /*
    Cuts cutMuSoft;
    cutMuSoft.selectCut("muSoft");
    const string anaMuSoft = cutMuSoft.getCut();
    doTEffPlotSetRecoB0(c, events, "eff_analMuSoft", "isol. soft #mu selection efficiency", "B^{0}", 25, anaAll, anaMuSoft);

    Cuts cutAccB0;
    cutAccB0.selectCut("acc04B0");
    const string anaAccB0 = cutAccB0.getCut();
    doTEffPlotSetRecoB0(c, events, "eff_analAcc", "isol. reco acceptance efficiency", "B^{0}", 25, anaAll, anaAccB0);

    // ------------------------------------------------------------------------- finish
    cout << "If you want to get one pdf then run the following command (requires a qorking copy of pdftk):" << endl;
    cout << "pdftk B0{eff_muAcc,eff_KsAcc,eff_muKsAcc,eff_HLTbarrEff,eff_muKsAcc_HLTbarrEff,eff_JpCand,eff_KsCand,eff_B0Cand,eff_hltMatch,eff_upto_hltMatch,eff_anal,eff_analMuSoft,eff_analAcc}_{p,d3,ct}.pdf"
	 << " cat output B0effplots.pdf" << endl;
    */
}

void doCtFitGenCutSequenceB0Final(const string filename)
{
    //TCanvas *c = new TCanvas("c", "c", 400, 400);
    //const string filename = "../data/run365.root";
    //const string filename = "../data/run337.root";
    //TFile *f = TFile::Open("../data/run337.root");
    //TTree *genevents = (TTree*)f->Get("genevents");
    //TTree *events = (TTree*)f->Get("events");
    const string pdfFilebase = "doCtFitGenCutSequenceB0";

    // ------------------------------------------------------------------------- step 0: all
    setUseEventsTree(false);
    doCtFitGen(filename, "", "B^{0} all", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_all.pdf").c_str());

    doCtFitGen(filename, "hasCand>0&&(isMCmatch&1)==1", "B^{0} cand match", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candmatch.pdf").c_str());

    // =====================================================================================================
    // ------------------------------------------------------------------------- step 7: analysis cuts
    setUseEventsTree(true);

    Cuts cutB0;
    cutB0.selectCut("B008", "acc06B0", "muSoft");
    const string anaAll = "isMCmatch==1&&isSig==1";
    const string anaPass = anaAll+"&&"+cutB0.getCut();
    const string anaPassBarrel = anaPass+"&&"+"HLTmatch==1&&HLTokBarrelJpsi==1";
    const string anaPassDispl = anaPass+"&&"+"HLTmatch==1&&HLTokDisplJpsi==1";

    doCtFitGen(filename, anaPass, "B^{0} analysis cuts", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacut.pdf").c_str());

    doCtFitGen(filename, anaPassBarrel, "B^{0} analysis cuts, barrel trigger", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutBarrel.pdf").c_str());

    doCtFitGen(filename, anaPassDispl, "B^{0} analysis cuts, displaced trigger", true, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutDispl.pdf").c_str());
}

void doCtFitGenCutSequenceLb(const string filename)
{
    //TCanvas *c = new TCanvas("c", "c", 400, 400);
    //const string filename = "../data/run365.root";
    //const string filename = "../data/run337.root";
    //TFile *f = TFile::Open("../data/run337.root");
    //TTree *genevents = (TTree*)f->Get("genevents");
    //TTree *events = (TTree*)f->Get("events");
    const string pdfFilebase = "doCtFitGenCutSequenceLb";

    // ------------------------------------------------------------------------- step 0: all
    setUseEventsTree(false);
    doCtFitGen(filename, "", "#Lambda_{b} all", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_all.pdf").c_str());

    //doCtFitGen(filename, "hasCand>0", "#Lambda_{b} cand", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_cand.pdf").c_str());

    doCtFitGen(filename, "hasCand>0&&(isMCmatch&1)==1", "#Lambda_{b} cand match", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candmatch.pdf").c_str());

    //doCtFitGen(filename, "hasCand>0&&(isMCmatch&1)==0", "#Lambda_{b} cand NOmatch", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_00_candnotmatch.pdf").c_str());

    /*
    // ------------------------------------------------------------------------- step 1: acceptance
    // Muon acceptance
    const string acc_mu1_1 = "TMath::Abs(etamu1)<1.2&&ptmu1>3.5";
    const string acc_mu1_2 = "TMath::Abs(etamu1)>=1.2&&TMath::Abs(etamu1)<1.6&&ptmu1>8-3.75*TMath::Abs(etamu1)";
    const string acc_mu1_3 = "TMath::Abs(etamu1)>=1.6&&TMath::Abs(etamu1)<2.4&&ptmu1>2.0";
    const string acc_mu2_1 = "TMath::Abs(etamu2)<1.2&&ptmu2>3.5";
    const string acc_mu2_2 = "TMath::Abs(etamu2)>=1.2&&TMath::Abs(etamu2)<1.6&&ptmu2>8-3.75*TMath::Abs(etamu2)";
    const string acc_mu2_3 = "TMath::Abs(etamu2)>=1.6&&TMath::Abs(etamu2)<2.4&&ptmu2>2.0";
    const string acc_mu = "(("+acc_mu1_1+")||("+acc_mu1_2+")||("+acc_mu1_3+"))&&(("+acc_mu2_1+")||("+acc_mu2_2+")||("+acc_mu2_3+"))";

    doCtFitGen(filename, acc_mu, "#Lambda_{b} Acceptance", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_01_muAcc.pdf").c_str());

    // ------------------------------------------------------------------------- step 2: Jpsi Kandidat
    const string jpCandAll = acc_mu;
    const string jpCandPass = "hasJpCand==1";

    doCtFitGen(filename, jpCandAll+"&&"+jpCandPass, "+J/#psi candidate", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_05_muKsAccTrgJpCand.pdf").c_str());

    //doCtFitGen(filename, jpCandPass, "#Lambda_{b} JpCand", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_05_isol_JpCand.pdf").c_str());

    // ------------------------------------------------------------------------- step 3: L0 Kandidat
    const string L0CandAll = jpCandAll + "&&" + jpCandPass;
    const string L0CandPass = "hasL0Cand==1";

    doCtFitGen(filename, L0CandAll+"&&"+L0CandPass, "#Lambda_{b} +#Lambda^{0} candidate", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_06_muL0AccTrgJpCandL0Cand.pdf").c_str());

    //doCtFitGen(filename, L0CandPass, "#Lambda_{b} L0Cand", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_06_isol_L0Cand.pdf").c_str());

    // ------------------------------------------------------------------------- step 4: B0 Kandidat
    const string lbCandAll = L0CandAll + "&&" + L0CandPass;
    const string lbCandPass = "hasCand==1";

    doCtFitGen(filename, lbCandAll+"&&"+lbCandPass, "#Lambda_{b} + #Lambda_{b} candidate unmatched", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_07_muL0AccTrgJpCandL0CandB0Cand.pdf").c_str());

    //doCtFitGen(filename, lbCandPass, "#Lambda_{b} B0Cand", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_07_isolB0Cand.pdf").c_str());

    doCtFitGen(filename, lbCandAll+"&&"+lbCandPass+"&&isMCmatch==1", "#Lambda_{b} #Lambda_{b} candidate MCmatched", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_07_isolB0CandMatch.pdf").c_str());

    // ------------------------------------------------------------------------- step 5: Trigger
    const string trigAll = lbCandAll + "&&" + lbCandPass;
    const string trigPass = "HLTokBarrelJpsi==1";
    const string trigPassAll = trigAll + "&&" + trigPass;

    doCtFitGen(filename, trigPassAll, "#Lambda_{b} + barrel trigger", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_04_muL0AccTrg.pdf").c_str());

    // ------------------------------------------------------------------------- step 6: Triggermatching
    const string hltMatchAll = lbCandAll + "&&" + lbCandPass;
    const string hltMatchPass = "HLTmatch==1";

    doCtFitGen(filename, hltMatchAll+"&&"+hltMatchPass, "#Lambda_{b} muL0 Acc Trg JpCand L0Cand B0Cand TrgMatch", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_08_muL0AccTrgJpCandL0CandB0CandHltmatch.pdf").c_str());

    // everything up to now:
    const string hltMatchPassAll = hltMatchAll + "&&" + hltMatchPass;

    */

    // =====================================================================================================
    // ------------------------------------------------------------------------- step 7: analysis cuts
    setUseEventsTree(true);

    Cuts cutLb;
    cutLb.selectCut("lb12", "acc06Lb", "muSoft");
    const string anaAll = "isMCmatch==1&&isSig==1";
    const string anaPass = anaAll+"&&"+cutLb.getCut();
    const string anaPassBarrel = anaPass+"&&"+"HLTmatch==1&&HLTokBarrelJpsi==1";
    const string anaPassDispl = anaPass+"&&"+"HLTmatch==1&&HLTokDisplJpsi==1";

    /*
    for (int i=0; i!=5; i++)
    {
	doCtFitGen(filename, "vrrs>"+toString(i), "#Lambda_{b}+#bar{#Lambda}_{b} reco all, vrrs>"+toString(i), false, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_09_recoAll_rvrrs_"+toString(i)+".pdf").c_str());

	doCtFitGen(filename, "rqha1>0&&vrrs>"+toString(i), "#Lambda_{b} reco all, vrrs>"+toString(i), false, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_09_recoAll_Lb_rvrrs_"+toString(i)+".pdf").c_str());

	doCtFitGen(filename, "rqha1<0&&vrrs>"+toString(i), "#bar{#Lambda}_{b} reco all, vrrs>"+toString(i), false, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_09_recoAll_Lbbar_rvrrs_"+toString(i)+".pdf").c_str());
    }
    */

    /*
    // Study with increasing cut on vrrs
    for (int i=0; i!=5; i++)
    {
	doCtFitGen(filename, anaPassBarrel+"&&vrrs>"+toString(i), "#Lambda_{b}+#bar{#Lambda}_{b} reco cuts, barrel trigger, vrrs>"+toString(i), false, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_10_anacutBarrel_rvrrs_"+toString(i)+".pdf").c_str());

	doCtFitGen(filename, anaPassBarrel+"&&rqha1>0&&vrrs>"+toString(i), "#Lambda_{b} reco cuts, barrel trigger, vrrs>"+toString(i), false, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_10_anacutBarrel_Lb_rvrrs_"+toString(i)+".pdf").c_str());

	doCtFitGen(filename, anaPassBarrel+"&&rqha1<0&&vrrs>"+toString(i), "#bar{#Lambda}_{b} reco cuts, barrel trigger, vrrs>"+toString(i), false, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_10_anacutBarrel_Lbbar_rvrrs_"+toString(i)+".pdf").c_str());
    }
    */

    doCtFitGen(filename, "rqha1<0", "#bar{#Lambda}_{b} reco", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_0.pdf").c_str());

    /*
    doCtFitGen(filename, anaAll+"&&rqha1<0", "#bar{#Lambda}_{b} reco, mcmatch, sig", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_27.pdf").c_str());

    doCtFitGen(filename, anaPass+"&&rqha1<0", "#bar{#Lambda}_{b} reco, mcmatch, sig, anacut", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_28.pdf").c_str());

    doCtFitGen(filename, anaPass+"&&rqha1<0", "#bar{#Lambda}_{b} reco, mcmatch, sig, acc", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_29.pdf").c_str());

    doCtFitGen(filename, anaPass+"&&rqha1<0", "#bar{#Lambda}_{b} reco, mcmatch, sig, softMu", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_30.pdf").c_str());
    */

    doCtFitGen(filename, anaPass+"&&rqha1<0", "#bar{#Lambda}_{b} reco, anacut ohne trg", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_31.pdf").c_str());

    std::vector<string> cutlist;
    cutlist.push_back("rptmu1");
    cutlist.push_back("rptmu2");
    cutlist.push_back("mrs");
    cutlist.push_back("Kshypo");
    cutlist.push_back("rptha1");
    cutlist.push_back("rptha2");
    cutlist.push_back("probrs");
    cutlist.push_back("alphars");
    cutlist.push_back("d3rs");
    cutlist.push_back("d3rs/d3Ers");
    cutlist.push_back("vrrs");

    int counter = 0;
    for (std::vector<string>::const_iterator it = cutlist.begin(); it!=cutlist.end(); it++)
    {
	Cuts cut;
	cut.selectCut("lb13", "acc06Lb", "muSoft");
	cut.removeOneCut(*it);
	doCtFitGen(filename, cut.getCut()+"&&rqha1<0", "#bar{#Lambda}_{b} reco, acc, soft, anacut ohne "+*it, false, 0e-12, 16e-12);
	gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_Remove"+toString(counter)+".pdf").c_str());
	counter++;
    }

    //doCtFitGen(filename, anaPassBarrel+"&&rqha1<0", "#bar{#Lambda}_{b} reco, anacut mit trg", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_32.pdf").c_str());

    /*
       // Studien für Bias durch Cuts
    doCtFitGen(filename, "rqha1<0&&HLTokBarrelJpsi==1", "#bar{#Lambda}_{b} reco, barrel", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_1.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&HLTokBarrelJpsi==1&&HLTmatch==1", "#bar{#Lambda}_{b} reco, barrel match", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_2.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&rptmu1>3", "#bar{#Lambda}_{b} reco, , rptmu1>3", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_3.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&rptmu2>3", "#bar{#Lambda}_{b} reco, , rptmu2>3", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_4.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&mrs>1.110683&&mrs<1.120683", "#bar{#Lambda}_{b} reco, ", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_5.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&(Kshypo<0.485614||Kshypo>0.509614)", "#bar{#Lambda}_{b} reco, ", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_6.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&rptha1>1.8", "#bar{#Lambda}_{b} reco, rptha1>1.8", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_7.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&rptha2>.5", "#bar{#Lambda}_{b} reco, rptha2>.5", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_8.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&probrs>.15", "#bar{#Lambda}_{b} reco, probrs>.15", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_9.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&alphars<.012", "#bar{#Lambda}_{b} reco, alphars<.012", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_10.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&d3rs>3", "#bar{#Lambda}_{b} reco, d3rs>3", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_11.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&d3rs/d3Ers>15", "#bar{#Lambda}_{b} reco, d3rs/d3Ers>15", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_12.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&vrrs>3", "#bar{#Lambda}_{b} reco, vrrs>3", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_13.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&TMath::Abs(etajp)<1.25", "#bar{#Lambda}_{b} reco, |etajp|<1.25", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_14.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>6.5", "#bar{#Lambda}_{b} reco, ptjp>6.5", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_15.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>10", "#bar{#Lambda}_{b} reco, ptjp>10", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_16.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>13", "#bar{#Lambda}_{b} reco, ptjp>13", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_17.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>15", "#bar{#Lambda}_{b} reco, ptjp>15", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_18.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>6.5&&TMath::Abs(etajp)<1.25", "#bar{#Lambda}_{b} reco, ptjp>6.5 + |etajp|<1.25", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_19.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>10&&TMath::Abs(etajp)<1.25", "#bar{#Lambda}_{b} reco, ptjp>10 + |etajp|<1.25", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_20.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>13&&TMath::Abs(etajp)<1.25", "#bar{#Lambda}_{b} reco, ptjp>13 + |etajp|<1.25", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_21.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&isCowboy==1", "#bar{#Lambda}_{b} reco, cowboy", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_22.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&isCowboy==0", "#bar{#Lambda}_{b} reco, seagull", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_23.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&maxdocajp<.5", "#bar{#Lambda}_{b} reco, maxdocajp<.5", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_24.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>10&&TMath::Abs(etajp)<1.25&&isCowboy==0", "#bar{#Lambda}_{b} reco, pt>10, eta<1.25, seagull", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_25.pdf").c_str());

    doCtFitGen(filename, "rqha1<0&&ptjp>13&&TMath::Abs(etajp)<1.25&&isCowboy==0", "#bar{#Lambda}_{b} reco, pt>13, eta<1.25, seagull", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_20_anacutTest_Lbbar_26.pdf").c_str());
    */


    //doCtFitGen(filename, "", "#Lambda_{b} reco all", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_09_recoAll.pdf").c_str());

    /*
    doCtFitGen(filename, anaAll, "#Lambda_{b} reco all matched sig", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_09_recoAllMatchSig.pdf").c_str());

    */
    doCtFitGen(filename, anaPass, "#Lambda_{b} analysis cuts", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacut.pdf").c_str());

    doCtFitGen(filename, anaPassBarrel, "#Lambda_{b} analysis cuts, barrel trigger", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutBarrel.pdf").c_str());

    doCtFitGen(filename, anaPassDispl, "#Lambda_{b} analysis cuts, displaced trigger", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutDispl.pdf").c_str());

    //doCtFitGen(filename, anaPass, "#Lambda_{b} reco HLT anaCut", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_11_recoHLTAnaCut.pdf").c_str());

    //doCtFitGen(filename, anaPass, "#Lambda_{b} reco HLT anaCut match", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_12_recoHLTAnaCutMatch.pdf").c_str());

    //doCtFitGen(filename, anaPass+"&&rpt1m>3&&rpt2m>3", "#Lambda_{b} reco HLT anaCut + ptmu>3", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_12_recoHLTAnaCut.pdf").c_str());

    //doCtFitGen(filename, anaPass+"&&(PvLip2==9999||PvLip2/PvLipE2>4)", "#Lambda_{b} reco HLT anaCut + 2ndPVSign>4", false, 0e-12, 16e-12);
    //gPad->SaveAs((pdfFilebase + "_13_recoHLTAnaCut.pdf").c_str());


    /*
    Cuts cutMuSoft;
    cutMuSoft.selectCut("muSoft");
    const string anaMuSoft = cutMuSoft.getCut();
    doTEffPlotSetRecoB0(c, events, "eff_analMuSoft", "isol. soft #mu selection efficiency", "#Lambda_{b}", 25, anaAll, anaMuSoft);

    Cuts cutAccB0;
    cutAccB0.selectCut("acc04B0");
    const string anaAccB0 = cutAccB0.getCut();
    doTEffPlotSetRecoB0(c, events, "eff_analAcc", "isol. reco acceptance efficiency", "#Lambda_{b}", 25, anaAll, anaAccB0);

    // ------------------------------------------------------------------------- finish
    cout << "If you want to get one pdf then run the following command (requires a qorking copy of pdftk):" << endl;
    cout << "pdftk B0{eff_muAcc,eff_KsAcc,eff_muKsAcc,eff_HLTbarrEff,eff_muKsAcc_HLTbarrEff,eff_JpCand,eff_KsCand,eff_B0Cand,eff_hltMatch,eff_upto_hltMatch,eff_anal,eff_analMuSoft,eff_analAcc}_{p,d3,ct}.pdf"
	 << " cat output B0effplots.pdf" << endl;
    */
}

void doCtFitGenCutSequenceLbFinal(const string filename)
{
    //TCanvas *c = new TCanvas("c", "c", 400, 400);
    //const string filename = "../data/run365.root";
    //const string filename = "../data/run337.root";
    //TFile *f = TFile::Open("../data/run337.root");
    //TTree *genevents = (TTree*)f->Get("genevents");
    //TTree *events = (TTree*)f->Get("events");
    const string pdfFilebase = "doCtFitGenCutSequenceLb";

    // ------------------------------------------------------------------------- step 0: all
    setUseEventsTree(false);
    doCtFitGen(filename, "", "#Lambda_{b} all", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_all.pdf").c_str());

    doCtFitGen(filename, "qha1>0", "#Lambda_{b} only", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_all_lb.pdf").c_str());

    doCtFitGen(filename, "qha1<0", "#bar{#Lambda}_{b} only", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_all_lbbar.pdf").c_str());

    doCtFitGen(filename, "hasCand>0&&(isMCmatch&1)==1", "#Lambda_{b} cand match", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candmatch.pdf").c_str());

    doCtFitGen(filename, "hasCand>0&&(isMCmatch&1)==1&&qha1>0", "#Lambda_{b} only cand match", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candmatch_lb.pdf").c_str());

    doCtFitGen(filename, "hasCand>0&&(isMCmatch&1)==1&&qha1<0", "#bar{#Lambda}_{b} only cand match", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_00_candmatch_lbbar.pdf").c_str());

    // =====================================================================================================
    // ------------------------------------------------------------------------- step 7: analysis cuts
    setUseEventsTree(true);

    Cuts cutLb;
    cutLb.selectCut("lb14", "acc06Lb", "muSoft");
    const string anaAll = "isMCmatch==1&&isSig==1";
    const string anaPass = anaAll+"&&"+cutLb.getCut();
    const string anaPassBarrel = anaPass+"&&"+"HLTmatch==1&&HLTokBarrelJpsi==1";
    const string anaPassDispl = anaPass+"&&"+"HLTmatch==1&&HLTokDisplJpsi==1";

    doCtFitGen(filename, anaPass, "#Lambda_{b} analysis cuts", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacut.pdf").c_str());

    doCtFitGen(filename, anaPass+"&&rqha1>0", "#Lambda_{b} only analysis cuts", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacut_lb.pdf").c_str());

    doCtFitGen(filename, anaPass+"&&rqha1<0", "#bar{#Lambda}_{b} only analysis cuts", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacut_lbbar.pdf").c_str());

    doCtFitGen(filename, anaPassBarrel, "#Lambda_{b} analysis cuts, barrel trigger", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutBarrel.pdf").c_str());

    doCtFitGen(filename, anaPassBarrel+"&&rqha1>0", "#Lambda_{b} only analysis cuts, barrel trigger", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutBarrel_lb.pdf").c_str());

    doCtFitGen(filename, anaPassBarrel+"&&rqha1<0", "#bar{#Lambda}_{b} only analysis cuts, barrel trigger", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutBarrel_lbbar.pdf").c_str());

    doCtFitGen(filename, anaPassDispl, "#Lambda_{b} analysis cuts, displaced trigger", false, 0e-12, 16e-12);
    gPad->SaveAs((pdfFilebase + "_10_anacutDispl.pdf").c_str());

}

