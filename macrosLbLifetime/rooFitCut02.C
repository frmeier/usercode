#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include "TMath.h"
#include <string>
#include "TLatex.h"
#include <sstream>
#include <iomanip>
#include "setTDRStyle_modified.C"

#include "doMassFitLb01.C"

#include "cut.h"
#include "cuts.C" // implies cut.C

using std::cout;
using std::endl;
using RooFit::Extended;
using RooFit::Name;
using RooFit::Title;
using RooFit::Binning;
using RooFit::Components;
using RooFit::LineStyle;

TCanvas *canvas;

struct addtlCuts
{
    addtlCuts(std::string mycut, std::string n) : cut(mycut), name(n) {};
    std::string cut;
    std::string name;
};

struct Fitresults
{
    std::vector<double> sig;
    std::vector<double> bgr;
    std::vector<double> soversqrtsb;
    std::vector<double> soverb;
    std::vector<double> mass;
    std::vector<double> width;
    std::vector<string> title;
};

void plotResult(TPad* pad, std::string histoName, const std::vector<double> &data, const std::vector<std::string> &names, std::string title)
{
    if (data.size() != names.size())
    {
        cout << "plotResult: data.size() != names.size() - aborting" << endl;
        throw ("plotResult: data.size() != names.size() - aborting");
    }
    const int nBins = data.size();
    TH1F *h = new TH1F(histoName.c_str(), title.c_str(), nBins, 0, nBins);
    for (unsigned int i=0; i!=nBins; i++)
    {
        if(!TMath::IsNaN(data[i])) h->SetBinContent(i+1,data[i]);
        h->GetXaxis()->SetBinLabel(i+1,names[i].c_str());
    }
    if (h->GetMinimum() > 0) h->SetMinimum(0);
    h->Draw("P0");
}

void plotResultErrors(TPad* pad, std::string histoName, const std::vector<double> &data, const std::vector<double> & errors, const std::vector<std::string> &names, std::string title)
{
    if (data.size() != names.size())
    {
        cout << "plotResult: data.size() != names.size() - aborting" << endl;
        throw ("plotResult: data.size() != names.size() - aborting");
    }
    const int nBins = data.size();
    TH1F *h = new TH1F(histoName.c_str(), title.c_str(), nBins, 0, nBins);
    for (unsigned int i=0; i!=nBins; i++)
    {
        if(!TMath::IsNaN(data[i])) h->SetBinContent(i+1,data[i]);
        if(!TMath::IsNaN(errors[i])) h->SetBinError(i+1,errors[i]);
        h->GetXaxis()->SetBinLabel(i+1,names[i].c_str());
    }
    //if (h->GetMinimum() > 0) h->SetMinimum(0);
    h->Draw("P0E");
}


void rooFitCut (std::string filename, std::string cutname, double lumi = 0, bool doPublicationGrade = false)
{
    const bool prelim(true);
    const std::string outpdf = "rooFitCut";
    const bool doL1(false);
    const bool doHLT(false);
    const bool doL1_HLT(false);
    const bool doMCtruth(false);
    const bool doOthers(false);
    const bool doEta(false);
    const bool doEtaMu(false);
    const bool doPt(false);
    const bool doRun(false);
    const bool doTrgEta(false);
    const bool doTrgMu(false);
    const bool doB0(false);
    const bool doIsSig(false); // does just this - useful for general plots for presentation of data (plot0) and MC (plot1)
    const bool doFirst2011(false);
    const bool do2011Av2(true);
    const bool do2011Av2_4(true);

    const bool noTitle(doPublicationGrade);
    const bool officialPlot(false);

    std::string exclCut("");

    // container for the results
    Fitresults fitresults;

    // Select cuts from cuts.C
    Cuts cut;
    //cut.selectCut(cutname,"acc03_16","HLT_matched_01");
    //cut.selectCut(cutname,"acc03","HLT_matched_2011_02");
    cut.selectCut(cutname,"acc03");
    //cut.selectCut(cutname);

    cout << "This is rooFitCut using the following cuts as a basis:" << endl;
    cout << cut.getCut() << endl;

    // Open file
    TFile* file = TFile::Open(filename.c_str());
    if (file ==0)
    {
        cout << "File " << filename << " not found. Exiting." << endl;
        return;
    }
    cout << "File " << filename << " succesfully opened." << endl;

    // Get tree
    std::string treename = "events";
    TTree *tree = (TTree*) file->Get(treename.c_str());
    if (tree == 0)
    {
        cout << "Unable to get TTree " << treename << " from file. Exiting." << endl;
        return;
    }
    cout << "TTree found with " << tree->GetEntries() << " entries." << endl;

    // here we define some additional cuts
    std::vector<addtlCuts> addtlCutsVec;
    addtlCutsVec.push_back(addtlCuts("","all data"));
    if (doL1)
    {
        addtlCutsVec.push_back(addtlCuts("L1TMu0==1","L1_SingleMu0"));
        addtlCutsVec.push_back(addtlCuts("L1TMu3==1","L1_SingleMu3"));
        addtlCutsVec.push_back(addtlCuts("L1TMu5==1","L1_SingleMu5"));
        addtlCutsVec.push_back(addtlCuts("L1TMu7==1","L1_SingleMu7"));
        addtlCutsVec.push_back(addtlCuts("L1TMu10==1","L1_SingleMu10"));
        addtlCutsVec.push_back(addtlCuts("L1TMu14==1","L1_SingleMu14"));
        addtlCutsVec.push_back(addtlCuts("L1TMu20==1","L1_SingleMu20"));
        addtlCutsVec.push_back(addtlCuts("L1TDMu0==1","L1_DoubleMuOpen"));
        addtlCutsVec.push_back(addtlCuts("L1TDMu3==1","L1_DoubleMu3"));
        addtlCutsVec.push_back(addtlCuts("((L1TMu0==1)||(L1TDMu0==1))","SMu0||DMu0"));
        addtlCutsVec.push_back(addtlCuts("((L1TMu3==1)||(L1TDMu3==1))","SMu3||DMu3"));
    }
    if (doL1_HLT)
    {
	const std::string cutL1 = "((L1TMu3==1)||(L1TDMu3==1))&&";
	cout << "Ask for additional L1 trigger(s): " << cutL1 << endl;
	//const std::string cutL1 = "((L1TMu0==1)||(L1TDMu0==1))&&";
	//cut.cs.addCut(new cutConst("((L1TMu3==1)||(L1TDMu3==1))")); cut.parvec.push_back(0);
        addtlCutsVec.push_back(addtlCuts(cutL1+"HLTqrk","HLTQuarkonium"));
        //addtlCutsVec.push_back(addtlCuts("HLTqrkLS","HLTQuarkoniumLS"));
        //addtlCutsVec.push_back(addtlCuts("HLTDMuOp","HLTDoubleMuOpen"));
        //addtlCutsVec.push_back(addtlCuts("HLTDMu0","HLTDoubleMu0"));
        //addtlCutsVec.push_back(addtlCuts("HLTDMu3","HLTDoubleMu3"));
        //addtlCutsVec.push_back(addtlCuts("HLTDMu5","HLTDoubleMu5"));
        //addtlCutsVec.push_back(addtlCuts("HLTMu3t3jp","HLTMu3Trck3Jpsi"));
        //addtlCutsVec.push_back(addtlCuts("HLTMu3t5jp","HLTMu3Trck5Jpsi"));
        //addtlCutsVec.push_back(addtlCuts("HLTMu5t0jp","HLTMu5Trck0Jpsi"));
        addtlCutsVec.push_back(addtlCuts(cutL1+"HLTMu0jp","HLTMu0OSTJpsi"));
        addtlCutsVec.push_back(addtlCuts(cutL1+"HLTMu0jpT","HLTMu0OSTJpsiTight"));
        //addtlCutsVec.push_back(addtlCuts("HLTMu3jp","HLTMu3OSTJpsi"));
        //addtlCutsVec.push_back(addtlCuts("HLTMu3jpT","HLTMu3OSTJpsiTight"));
        //addtlCutsVec.push_back(addtlCuts("HLTMu5jp","HLTMu5OSTJpsi"));
        //addtlCutsVec.push_back(addtlCuts("HLTMu5jpT","HLTMu5OSTJpsiTight"));
        //addtlCutsVec.push_back(addtlCuts("HLTL1DMu0","HLTL1DoubleMu0"));
        //addtlCutsVec.push_back(addtlCuts("HLTL2DMu0","HLTL2DoubleMu0"));
        //addtlCutsVec.push_back(addtlCuts("HLTL2Mu0","HLTL2Mu0"));
	// composite
        //addtlCutsVec.push_back(addtlCuts("(HLTqrk||HLTDMu3||HLTMu0jp||HLTMu0jpT||HLTMu3jp||HLTMu3jpT)","HLT 1,5,10-13"));
        //addtlCutsVec.push_back(addtlCuts("(HLTqrk||HLTDMu3)","HLTQuarkonium||HLTDoubleMu3"));
        //addtlCutsVec.push_back(addtlCuts(cutL1+"(HLTqrk||HLTMu0jp||HLTMu0jpT)","HLTQuarkonium||HLTMu0OSTJpsi||HLTMu0OSTJpsiTight"));
        addtlCutsVec.push_back(addtlCuts(cutL1+"(HLTMu0jp||HLTMu0jpT)","HLTMu0OSTJpsi||HLTMu0OSTJpsiTight"));
        addtlCutsVec.push_back(addtlCuts("(HLTMu0jp||HLTMu0jpT)","HLTMu0OSTJpsi||HLTMu0OSTJpsiTight w/o L1"));
        //addtlCutsVec.push_back(addtlCuts("(HLTMu0jp||HLTMu0jpT||HLTMu3jp||HLTMu3jpT)","HLT 10-13"));
    }
    if (doHLT)
    {
	std::string rStart(""), rEnd("");
	//rStart = "140058"; rEnd = "144114";
	//rStart = "146428"; rEnd = "147116";
	//rStart = "147196"; rEnd = "148058";
	//rStart = "148819"; rEnd = "149182";
	//rStart = "149291"; rEnd = "149294";
	std::string addCut = rStart.size()>0 ? "run>=" + rStart + "&&run<=" + rEnd : "";
	std::string addCutT = " runs " + rStart + " to " + rEnd;
        //addtlCutsVec.push_back(addtlCuts("HLTL2Mu0","HLTL2Mu0"));
	if (addCut.size()>0) addtlCutsVec.push_back(addtlCuts(addCut,addCutT));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu0TkMu0jp==1","HLT_Mu0_TkMu0_Jpsi"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu0TkMu0jpNC==1","HLT_Mu0_TkMu0_Jpsi_NoCharge"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu3TkMu0jp==1","HLT_Mu3_TkMu0_Jpsi"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu3TkMu0jpNC==1","HLT_Mu3_TkMu0_Jpsi_NoCharge"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu5TkMu0jp==1","HLT_Mu5_TkMu0_Jpsi"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu3t3jp==1","HLT_Mu3_Track3_Jpsi"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu3t3jp2==1","HLT_Mu3_Track3_Jpsi_v2"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu3t3jp3==1","HLT_Mu3_Track3_Jpsi_v3"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTqrk==1","HLT_DoubleMu0_Quarkonium_v1"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu0jp==1","HLT_Mu0_TkMu0_OST_Jpsi"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu0jpT1==1","HLT_Mu0_TkMu0_OST_Jpsi_Tight_v1"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu0jpT2==1","HLT_Mu0_TkMu0_OST_Jpsi_Tight_v2"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTMu0jpT3==1","HLT_Mu0_TkMu0_OST_Jpsi_Tight_v3"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTok==1","HLTok"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"((run>=140058&&run<=144114&&HLTMu0TkMu0jp==1)||(run>=146428&&run<=148058&&HLTMu0jp==1)||(run>=148819&&run<=149182&&HLTMu0jpT2==1)||(run>=149291&&run<=149294&&HLTMu0jpT3==1))","HLT combined"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTmatch==1","HLTmatch"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"HLTmatch==1&&HLTok==1","HLTmatch and HLTok"+(addCut.size()>0?" ("+addCutT+")":"")));
    }
    if (doMCtruth)
    {
        addtlCutsVec.push_back(addtlCuts("isMCmatch==1","MCmatch"));
        addtlCutsVec.push_back(addtlCuts("isMCmatch==0","noMCmatch"));
        addtlCutsVec.push_back(addtlCuts("isSig==1","sig"));
        addtlCutsVec.push_back(addtlCuts("isSig==0","nosig"));
        addtlCutsVec.push_back(addtlCuts("isMCmatch==1&&isSig==1","MCmatch + sig"));
        addtlCutsVec.push_back(addtlCuts("isMCmatch==0&&isSig==1","noMCmatch + sig"));
        addtlCutsVec.push_back(addtlCuts("isMCmatch==1&&isSig==0","MCmatch + nosig"));
    }
    if (doEta)
    {
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(etalb)<1","|#eta(#Lambda_b)|<1"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(etalb)>1","|#eta(#Lambda_b)|>1"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(etalb)<1.6","|#eta(#Lambda_b)|<1.6"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(etalb)>1.6","|#eta(#Lambda_b)|>1.6"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(etalb)>1&&TMath::Abs(etalb)<1.6","1<|#eta(#Lambda_b)|<1.6"));
    }
    if (doEtaMu)
    {
        addtlCutsVec.push_back(addtlCuts("isSig==1","sig"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)<1.0&&TMath::Abs(Seta2m)<1.0","|#eta(#mu_{1,2})|<1.0"));
        addtlCutsVec.push_back(addtlCuts("isSig==1&&TMath::Abs(Seta1m)<1.0&&TMath::Abs(Seta2m)<1.0","|#eta(#mu_{1,2})|<1.0"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)<1.6&&TMath::Abs(Seta2m)<1.6","|#eta(#mu_{1,2})|<1.6"));
        addtlCutsVec.push_back(addtlCuts("isSig==1&&TMath::Abs(Seta1m)<1.6&&TMath::Abs(Seta2m)<1.6","|#eta(#mu_{1,2})|<1.6"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)<2.0&&TMath::Abs(Seta2m)<2.0","|#eta(#mu_{1,2})|<2.0"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)>1&&TMath::Abs(Seta2m)>1","|#eta(#mu_{1,2})|>1.0"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)>1.6&&TMath::Abs(Seta2m)>1.6","|#eta(#mu_{1,2})|>1.6"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)>2&&TMath::Abs(Seta2m)>2","|#eta(#mu_{1,2})|>2.0"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)>1&&TMath::Abs(Seta1m)<2&&TMath::Abs(Seta2m)>1&&TMath::Abs(Seta2m)<2","1<|#eta(#mu_{1,2})|<2"));
        addtlCutsVec.push_back(addtlCuts("(TMath::Abs(Seta1m)<1.0&&TMath::Abs(Seta2m)>1.0)||(TMath::Abs(Seta2m)<1.0&&TMath::Abs(Seta1m)>1.0)","|#eta(#mu_i)|<1.0, |#eta(#mu_j)|>1.0"));
    }
    if (doPt)
    {
        addtlCutsVec.push_back(addtlCuts("ptlb>1","p_{T}(#Lambda_b)>1"));
        addtlCutsVec.push_back(addtlCuts("ptlb>5","p_{T}(#Lambda_b)>5"));
        addtlCutsVec.push_back(addtlCuts("ptlb>10","p_{T}(#Lambda_b)>10"));
        addtlCutsVec.push_back(addtlCuts("ptlb>15","p_{T}(#Lambda_b)>15"));
        addtlCutsVec.push_back(addtlCuts("ptlb>20","p_{T}(#Lambda_b)>20"));
    }
    if (doRun)
    {
	std::string etacut = "1.0";
	std::string addCut = "TMath::Abs(Seta1m)<"+etacut+"&&TMath::Abs(Seta2m)<"+etacut;
	std::string addCutT = "|#eta(#mu_{1,2})|<"+etacut;
	if (addCut.size()>0) addtlCutsVec.push_back(addtlCuts(addCut,addCutT));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"run>=140058&&run<=144114","runs [140558,144114]"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"run>=146428&&run<=147116","runs [146428,147116]"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"run>=147196&&run<=148058","runs [147196,148058]"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"run>=148819&&run<=149182","runs [148819,149182]"+(addCut.size()>0?" ("+addCutT+")":"")));
        addtlCutsVec.push_back(addtlCuts((addCut.size()>0?addCut+"&&":"")+"run>=149291&&run<=149294","runs [149291,149294]"+(addCut.size()>0?" ("+addCutT+")":"")));
    }
    if (doTrgEta)
    {
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1","HLT passed and matched"));
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1&&TMath::Abs(Seta1m)<1.6&&TMath::Abs(Seta2m)<1.6","trigger matched, |#eta(#mu_{1,2})|<1.6"));
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1&&TMath::Abs(Seta1m)<1.0&&TMath::Abs(Seta2m)<1.0","trigger matched, |#eta(#mu_{1,2})|<1.0"));
    }
    if (doTrgMu) // needs a cutset without a cut on the muid
    {
	const int muId = 2;
	std::cout << "muIdStr: " << muIdStr(2) << endl;
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1","HLT passed and matched"));
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1&&"+muIdStr(2),"ditto, muId==2 (global)"));
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1&&"+muIdStr(4),"ditto, muId==4 (tracker)"));
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1&&"+muIdStr(6),"ditto, muId==4&6 (global and tracker)"));
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1&&"+muIdStr(16),"ditto, muId==16 (tracker arbitrated)"));
    }
    if (doB0)
    {
	exclCut = "Kshypo";
	addtlCutsVec.push_back(addtlCuts("isSig","+isSig"));
	//Cuts cutsAcc, cutsHLT;
	//cutsAcc.selectCut("acc03");
	//cutsHLT.selectCut("HLT_matched_01");
	//addtlCutsVec.push_back(addtlCuts(cutsAcc.getCut(),"+acc03"));
	//addtlCutsVec.push_back(addtlCuts(cutsHLT.getCut(),"+HLT_matched_01"));
	//addtlCutsVec.push_back(addtlCuts(cutsAcc.getCut()+"&"+cutsHLT.getCut(),"+acc03+HLT_matched_01"));
	addtlCutsVec.push_back(addtlCuts(cut.getOneCut(exclCut),"Ks window excluded"));
	addtlCutsVec.push_back(addtlCuts("isSig&"+cut.getOneCut(exclCut),"isSig + Ks window excluded"));
	//addtlCutsVec.push_back(addtlCuts(cut.getOneCut("Kshypo")+"&"+cutsAcc.getCut(),"+Kshypo+acc03"));
	//addtlCutsVec.push_back(addtlCuts(cut.getOneCut("Kshypo")+"&"+cutsHLT.getCut(),"+Kshypo+HLT_matched_01"));
	//addtlCutsVec.push_back(addtlCuts(cut.getOneCut("Kshypo")+"&"+cutsAcc.getCut()+"&"+cutsHLT.getCut(),"+Kshypo+acc03+HLT_matched_01"));
    }
    if (doIsSig)
    {
	addtlCutsVec.push_back(addtlCuts("isSig","+isSig"));
    }
    if (doFirst2011)
    {
	addtlCutsVec.push_back(addtlCuts("HLTqrk==1","HLT_DoubleMu3_Quarkonium_v1"));
	addtlCutsVec.push_back(addtlCuts("HLTqrk==1&&HLTmatch==1","HLT_DoubleMu3_Quarkonium_v1 matched"));
	addtlCutsVec.push_back(addtlCuts("HLTDMu3jp==1","HLT_DoubleMu3_Jpsi_v1"));
	addtlCutsVec.push_back(addtlCuts("HLTDMu3jp==1&&HLTmatch==1","HLT_DoubleMu3_Jpsi_v1 matched"));
	const std::string runcut = "(run==161311||run==161312)";
	const std::string runcutTitle = "Runs 161311/2";
	addtlCutsVec.push_back(addtlCuts(runcut,runcutTitle));
	addtlCutsVec.push_back(addtlCuts(runcut+"&&HLTqrk==1","HLT_DoubleMu3_Quarkonium_v1 "+runcutTitle));
	addtlCutsVec.push_back(addtlCuts(runcut+"&&HLTqrk==1&&HLTmatch==1","HLT_DoubleMu3_Quarkonium_v1 matched "+runcutTitle));
	addtlCutsVec.push_back(addtlCuts(runcut+"&&HLTDMu3jp==1","HLT_DoubleMu3_Jpsi_v1 "+runcutTitle));
	addtlCutsVec.push_back(addtlCuts(runcut+"&&HLTDMu3jp==1&&HLTmatch==1","HLT_DoubleMu3_Jpsi_v1 matched "+runcutTitle));
    }
    if (do2011Av2)
    {
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1","all data, trigger matched"));

	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5BarJp==1","matched HLT_Dimuon6p5_Barrel_Jpsi_v1"));
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5JpDis==1","matched HLT_Dimuon6p5_Jpsi_Displaced_v1"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5Jp==1","matched HLT_Dimuon6p5_Jpsi_v1"));
	std::string etacut = "1.5";
	std::string addCut = "TMath::Abs(Seta1m)<"+etacut+"&&TMath::Abs(Seta2m)<"+etacut;
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5JpDis==1&&"+addCut,"matched HLT_Dimuon6p5_Jpsi_Displaced_v1, |#eta(#mu_{i})|<"+etacut));

	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5BarJp==1&&HLTDMu6p5JpDis==1","matched Barrel_Jpsi AND Jpsi_Displaced"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5BarJp==1&&HLTDMu6p5Jp==1","matched Barrel_Jpsi AND Mu6p5Jp"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5JpDis==1&&HLTDMu6p5Jp==1","matched Jpsi_Displaced AND Mu6p5Jp"));

	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&(HLTDMu6p5BarJp==1||HLTDMu6p5JpDis==1)","matched Barrel_Jpsi OR Jpsi_Displaced"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&(HLTDMu6p5BarJp==1||HLTDMu6p5Jp==1)","matched Barrel_Jpsi OR Mu6p5Jp"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&(HLTDMu6p5JpDis==1||HLTDMu6p5Jp==1)","matched Jpsi_Displaced OR Mu6p5Jp"));

	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5BarJp==1&&HLTDMu6p5JpDis==1&&HLTDMu6p5Jp==1","matched all 3 ANDed"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&(HLTDMu6p5BarJp==1||HLTDMu6p5JpDis==1||HLTDMu6p5Jp==1)","matched all 3 ORed"));

	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTMu5L2Mu2Jpsi==1","matched HLT_Mu5_L2Mu2_Jpsi_v3"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTMu5Tr2Jpsi==1","matched HLT_Mu5_Track2_Jpsi_v2"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTMu5Tr7Jpsi==1","matched HLT_Mu7_Track7_Jpsi_v3"));
    }
    if (do2011Av2_4)
    {
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1","all data, trigger matched"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5BarJp==1","matched HLT_Dimuon6p5_Barrel_Jpsi_v1"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5JpDis==1","matched HLT_Dimuon6p5_Jpsi_Displaced_v1"));
	//addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu6p5Jp==1","matched HLT_Dimuon6p5_Jpsi_v1"));
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu10BarJp==1","matched HLT_Dimuon10_Barrel_Jpsi_v1"));
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu7JpDis==1","matched HLT_Dimuon7_Jpsi_Displaced_v1"));
	std::string etacut = "1.0";
	std::string addCut = "TMath::Abs(Seta1m)<"+etacut+"&&TMath::Abs(Seta2m)<"+etacut;
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu7JpDis==1&&"+addCut,"matched HLT_Dimuon7_Jpsi_Displaced_v1, |#eta(#mu_{i})|<"+etacut));
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu10BarJp==1&&"+addCut,"matched HLT_Dimuon10_Barrel_Jpsi_v1, |#eta(#mu_{i})|<"+etacut));
    }

    //addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)<1.0&&TMath::Abs(Seta2m)<1.0","|#eta(#mu_{1,2})|<1.0"));
    // initialise the TCanvas
    setTDRStyle();
    const int maxCanvasCols = doPublicationGrade ? 1 : 3;
    const int maxCanvasRows = doPublicationGrade ? 1 : 5;
    const int nAdditionalPlots = 5;
    const int nPlots = addtlCutsVec.size() + nAdditionalPlots;
    int nCanvasCols = maxCanvasCols;
    int nCanvasRows = maxCanvasRows;
    if (nPlots < maxCanvasCols*maxCanvasRows)
    {
        if(nPlots == 1)
            nCanvasCols = nCanvasCols = 1;
        else
        {
            nCanvasCols = maxCanvasCols;
            nCanvasRows = TMath::Ceil(nPlots / nCanvasCols);
        }
    }
    int nCanvasCd = nCanvasCols*nCanvasRows;
    //canvas = new TCanvas("canvas","Roofit",nCanvasCols*500,nCanvasRows*350);
    canvas = new TCanvas("canvas","Roofit",nCanvasCols*500,nCanvasRows*300);
    canvas->Divide(nCanvasCols,nCanvasRows);

    // check if we can make a preselection
    TFile* tmpfile = new TFile("tmp.root","RECREATE");
    if (tmpfile == 0)
    {
	cout << "temporary file creation failed - exiting" << endl;
	return;
    }
    cout << "tmpfile: " << tmpfile << endl;
    TTree* treepresel = 0;
    bool usePreseltree(false);
    if(addtlCutsVec.size()>0 && addtlCutsVec[0].cut=="")
    {
        cout << "Preselecting events..." << endl;
	if (0 == exclCut.size())
	    treepresel = tree->CopyTree(cut.getCut().c_str());
	else
	    treepresel = tree->CopyTree(cut.getCutExceptOne(exclCut).c_str());
        cout << "now " << treepresel->GetEntries() << " of " << tree->GetEntries()<< endl;
        delete tree;
        usePreseltree = true; // needed to detach the preseltree from the original one
    }

    int canvasCdCounter(0), canvasPageCounter(0);
    for (unsigned int i=0; i!=addtlCutsVec.size(); i++)
    {
        canvasCdCounter++;
        if(canvasCdCounter > nCanvasCd)
        {
            canvas->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
            canvas->SaveAs((outpdf + toString(canvasPageCounter) + ".png").c_str());
            canvasPageCounter++;
            canvasCdCounter = 1;
            canvas->Clear("D");
        }
        canvas->cd(canvasCdCounter);

        // Apply cut selection
	std::string curCut;
	if (0 == exclCut.size())
	    curCut = addtlCutsVec[i].cut=="" ? cut.getCut() : (cut.getCut() + "&&" + addtlCutsVec[i].cut);
	else
	    curCut = addtlCutsVec[i].cut=="" ? cut.getCutExceptOne(exclCut) : (cut.getCutExceptOne(exclCut) + "&&" + addtlCutsVec[i].cut);
	cout << "Using additional cut: " << addtlCutsVec[i].cut << " name: " << addtlCutsVec[i].name << endl;
        TTree* subtree;
        if(usePreseltree)
        {
            subtree = treepresel->CopyTree(curCut.c_str());
            cout << "subtree: " << subtree->GetEntries() << " entries of " << treepresel->GetEntries()<< endl;
        }
        else
        {
            subtree = tree->CopyTree(curCut.c_str());
            cout << "subtree: " << subtree->GetEntries() << " entries of " << tree->GetEntries()<< endl;
        }
	const double curEntries = subtree->GetEntries();
	cout << "curEntries: " << curEntries << endl;

	doMassFitLb01_fitresults res = doMassFitLb01(subtree, lumi, addtlCutsVec[i].name, noTitle, prelim, officialPlot);

        // collect results
        fitresults.sig.push_back(res.sig);
        fitresults.bgr.push_back(res.bgr);
        fitresults.soversqrtsb.push_back(res.soversqrtsb);
        fitresults.soverb.push_back(res.soverb);
        fitresults.mass.push_back(res.mass);
        fitresults.width.push_back(res.width);
        fitresults.title.push_back(addtlCutsVec[i].name);
    }

    // draw result histos
    for (unsigned int i=0; i!= nAdditionalPlots; i++)
    {
        canvasCdCounter++;
        if(canvasCdCounter > nCanvasCd)
        {
            canvas->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
            canvas->SaveAs((outpdf + toString(canvasPageCounter) + ".png").c_str());
            canvasPageCounter++;
            canvasCdCounter = 1;
            canvas->Clear("D");
        }
        TPad *pad = (TPad*)canvas->cd(canvasCdCounter);
        if (0 == i) plotResult(pad, "hsig", fitresults.sig, fitresults.title, "Signal");
        if (1 == i) plotResult(pad, "hbgr", fitresults.bgr, fitresults.title, "Background");
        if (2 == i) plotResult(pad, "hssqtsb", fitresults.soversqrtsb, fitresults.title, "S/#sqrt(S+B)");
        if (3 == i) plotResult(pad, "hsoverb", fitresults.soverb, fitresults.title, "S/B");
        if (4 == i) plotResultErrors(pad, "hmass", fitresults.mass, fitresults.width, fitresults.title, "Mass");
    }

    // finalize page
    canvas->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
    canvas->SaveAs((outpdf + toString(canvasPageCounter) + ".png").c_str());
    delete canvas;
}

