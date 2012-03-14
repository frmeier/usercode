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

//#include "cut.C"
#include "cuts.C" // implies cut.C

using std::cout;
using std::endl;
using RooFit::Extended;
using RooFit::Name;
using RooFit::Title;
using RooFit::Binning;
using RooFit::Components;
using RooFit::LineStyle;

TCanvas *c;

struct addtlCuts
{
    addtlCuts(std::string c, std::string n) : cut(c), name(n) {};
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

std::string roundToString(double v, std::streamsize precision)
{
    ostringstream oss;
    oss << std::setprecision(precision) << std::fixed << v;
    return oss.str();
}

void setTLatexInit(TLatex *txt)
{
    txt->SetTextSize(0.04);
    txt->SetNDC(true);
}

TLatex* writeTLatex(std::string text, double x, double y)
{
    TLatex* txt = new TLatex(x,y,text.c_str());
    txt->SetTextSize(0.04);
    txt->SetNDC(true);
    return txt;
}

std::string muIdStr(const int &id)
{
    const std::string sid = toString(id);
    return "(rid1m&"+sid+")=="+sid+"&&(rid2m&"+sid+")=="+sid;
}

void rooFitCut (std::string filename, std::string cutname)
{
    const std::string outpdf = "rooFitCut";
    const bool doL1(false);
    const bool doHLT(false);
    const bool doHLTbrief(true);
    const bool doL1_HLT(false);
    const bool doMCtruth(false);
    const bool doOthers(false);
    const bool doEta(false);
    const bool doEtaMu(false);
    const bool doPt(false);
    const bool doRun(false);
    const bool doTrgEta(false);
    const bool doTrgMu(false);

    // container for the results
    Fitresults fitresults;

    // Select cuts from cuts.C
    Cuts cut;
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
    if (doHLTbrief) // meanwhile not so brief anymore ;)
    {
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTok==1","HLTmatch and HLTok"));
        addtlCutsVec.push_back(addtlCuts("HLTmatch==1","HLTmatch"));
        addtlCutsVec.push_back(addtlCuts("HLTDMu0==1","HLTDoubleMu0"));
        addtlCutsVec.push_back(addtlCuts("HLTDMu0==1&&HLTmatch==1","HLTDoubleMu0 and HLTmatched (caution)"));
        addtlCutsVec.push_back(addtlCuts("HLTDMu3==1","HLTDoubleMu3"));
        addtlCutsVec.push_back(addtlCuts("HLTDMu3==1&&HLTmatch==1","HLTDoubleMu3 and HLTmatched (caution)"));
        addtlCutsVec.push_back(addtlCuts("HLTqrk==1","HLTDoubleMu0_Quarkonium_v1"));
        addtlCutsVec.push_back(addtlCuts("HLTqrk==1&&HLTmatch==1","HLTDoubleMu0_Quarkonium_v1 and HLTmatched (caution)"));
        addtlCutsVec.push_back(addtlCuts("((HLTDMu0==1&&run<=147196)||(HLTqrk==1&&run>=147196))","HLTDoubleMu0 (-147196), later HLTQrk"));
        addtlCutsVec.push_back(addtlCuts("((HLTDMu0==1&&run<=147196)||(HLTqrk==1&&run>=147196))&&HLTmatch==1","HLTDoubleMu0/HLTQrk and HLTmatched"));
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
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)<1.0&&TMath::Abs(Seta2m)<1.0","|#eta(#mu_{1,2})|<1.0"));
        addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)<1.6&&TMath::Abs(Seta2m)<1.6","|#eta(#mu_{1,2})|<1.6"));
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


    // initialise the TCanvas
    setTDRStyle();
    const int maxCanvasCols = 1;
    const int maxCanvasRows = 1;
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
    c = new TCanvas("c","Roofit",nCanvasCols*500,nCanvasRows*350);
    c->Divide(nCanvasCols,nCanvasRows);

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
        treepresel = tree->CopyTree(cut.getCut().c_str());
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
            c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
            canvasPageCounter++;
            canvasCdCounter = 1;
            c->Clear("D");
        }
        c->cd(canvasCdCounter);

        // Apply cut selection
        std::string curCut = addtlCutsVec[i].cut=="" ? cut.getCut() : (cut.getCut() + "&&" + addtlCutsVec[i].cut);
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

	// some configs
	const double valLo(5.22), valHi(6.045);
	//const double valLo(5.35), valHi(5.90);
	const int nBins(33);
	//const int nBins(22);

        RooRealVar mass("mlb", "J/#psi #Lambda mass [GeV]", valLo, valHi);
        RooDataSet data("data", "Lambda_b dataset", subtree, mass);

        //RooRealVar lb_lifetime("lb_lifetime", "c#tau #Lambda_{b} [mm]",0.25, 3.5);
        //RooDataSet lt_data("lt_data", "#Lambda_{b} lifetime dataset", tree, lb_lifetime);

        //   mass fit
        RooRealVar mean("mean","mean",5.62,5.60,5.64) ;
        RooRealVar sigma("sigma","sigma",0.02,0.01,0.03) ;
        RooGaussian sig("sig","signal p.d.f.",mass,mean,sigma) ;

        RooRealVar c0("c0","coefficient #0", 1.0,-1,2) ;
        RooRealVar c1("c1","coefficient #1", 0.1,-1,1) ;
        RooPolynomial bkg("bkg","background p.d.f.",mass,RooArgList(c0,c1),0) ;

        RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries) ;
        RooRealVar nbkg("nbkg","Background fraction", .01* curEntries,0.,curEntries*1.4) ;

        RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg)) ;
        model.fitTo(data,Extended(kTRUE));

        RooPlot* xframe = mass.frame(Name("xframe"),Title(("Mass of #Lambda_{b} - "+ addtlCutsVec[i].name).c_str())) ;
        data.plotOn(xframe,Binning(nBins)) ;
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

        // write the results to the canvas
        TLatex *txt;
        const double txtPosLeft = .7;
        const double txtPosTop = .80;
        const double txtLineSpace = .05;
        xframe->addObject(writeTLatex("n Signal: " + roundToString(nsigWindow,1),txtPosLeft,txtPosTop-0*txtLineSpace));
        xframe->addObject(writeTLatex("n Bgr: " + roundToString(nbkgWindow,1), txtPosLeft,txtPosTop-1*txtLineSpace));
        xframe->addObject(writeTLatex("S/#sqrt{S+B}: " + roundToString(nsigWindow/sqrt(nsigWindow+nbkgWindow),1),txtPosLeft,txtPosTop-2*txtLineSpace));
        xframe->addObject(writeTLatex("S/B: " + roundToString(nsigWindow/nbkgWindow,2),txtPosLeft,txtPosTop-3*txtLineSpace));
        xframe->addObject(writeTLatex("mass: " + roundToString(m,3) + " GeV/c^{2}", txtPosLeft,txtPosTop-4*txtLineSpace));
        xframe->addObject(writeTLatex("width: " + roundToString(s,3) + " GeV/c^{2}",txtPosLeft,txtPosTop-5*txtLineSpace));
        xframe->Draw();
        //c->Update();

        // collect results
        fitresults.sig.push_back(nsigWindow);
        fitresults.bgr.push_back(nbkgWindow);
        fitresults.soversqrtsb.push_back(nsigWindow/sqrt(nsigWindow+nbkgWindow));
        fitresults.soverb.push_back(nsigWindow/nbkgWindow);
        fitresults.mass.push_back(m);
        fitresults.width.push_back(s);
        fitresults.title.push_back(addtlCutsVec[i].name);
    }

    // draw result histos
    for (unsigned int i=0; i!= nAdditionalPlots; i++)
    {
        canvasCdCounter++;
        if(canvasCdCounter > nCanvasCd)
        {
            c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
            canvasPageCounter++;
            canvasCdCounter = 1;
            c->Clear("D");
        }
        TPad *pad = (TPad*)c->cd(canvasCdCounter);
        if (0 == i) plotResult(pad, "hsig", fitresults.sig, fitresults.title, "Signal");
        if (1 == i) plotResult(pad, "hbgr", fitresults.bgr, fitresults.title, "Background");
        if (2 == i) plotResult(pad, "hssqtsb", fitresults.soversqrtsb, fitresults.title, "S/#sqrt(S+B)");
        if (3 == i) plotResult(pad, "hsoverb", fitresults.soverb, fitresults.title, "S/B");
        if (4 == i) plotResultErrors(pad, "hmass", fitresults.mass, fitresults.width, fitresults.title, "Mass");
    }

    // finalize page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
    delete c;
}

