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
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "TMath.h"
#include <string>
#include "TLatex.h"
#include <sstream>
#include <iomanip>
#include "setTDRStyle_modified.C"

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
using namespace RooFit;

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
    for (int i=0; i!=nBins; i++)
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
    for (int i=0; i!=nBins; i++)
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

void rooFitCt01 (std::string filename, std::string cutname, bool doPublicationGrade = false)
{
    const std::string outpdf = "rooFitCt";
    const bool doIsSig(false); // does just this - useful for general plots for presentation of data (plot0) and MC (plot1)
    const bool doData2011(true);

    const bool noTitle(doPublicationGrade);

    std::string exclCut("");

    // container for the results
    Fitresults fitresults;

    // Select cuts from cuts.C
    Cuts cut;
    //cut.selectCut(cutname,"acc03","HLT_matched_01");
    cut.selectCut(cutname,"acc03","HLT_matched_2011_02");
    //cut.selectCut(cutname);
    //cut.selectCut("nocuts");

    cout << "This is rooFitCt01 using the following cuts as a basis:" << endl;
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
    //addtlCutsVec.push_back(addtlCuts("","all data"));
    if (doIsSig)
    {
	//addtlCutsVec.push_back(addtlCuts("isSig","+isSig"));
	//addtlCutsVec.push_back(addtlCuts("isSig==1&&isMCmatch==1&&ct3dlb>1e-12","isSig&&isMCmatch"));
	addtlCutsVec.push_back(addtlCuts("isSig==1&&isMCmatch==1&&TMath::Abs(ct3dlbE/ct3dlb)<2","isSig&&isMCmatch"));
	//addtlCutsVec.push_back(addtlCuts("isSig==1&&isMCmatch==1&&ct3dlb>-2e-12&&ct3dlb<100e-12&&ct3dlbE<100e-12","isSig&&isMCmatch"));
    }
    if (doData2011)
    {
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu7JpDis==1","matched HLT_Dimuon7_Jpsi_Displaced_v1"));
	std::string etacut = "1.0";
	std::string addCut = "TMath::Abs(Seta1m)<"+etacut+"&&TMath::Abs(Seta2m)<"+etacut;
	addtlCutsVec.push_back(addtlCuts("HLTmatch==1&&HLTDMu7JpDis==1&&"+addCut,"matched HLT_Dimuon7_Jpsi_Displaced_v1, |#eta(#mu_{i})|<"+etacut));
    }

    //addtlCutsVec.push_back(addtlCuts("TMath::Abs(Seta1m)<1.0&&TMath::Abs(Seta2m)<1.0","|#eta(#mu_{1,2})|<1.0"));
    // initialise the TCanvas
    setTDRStyle();
    const int maxCanvasCols = doPublicationGrade ? 1 : 2;
    const int maxCanvasRows = doPublicationGrade ? 1 : 2;
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
            nCanvasRows = (int)TMath::Ceil(nPlots / nCanvasCols);
        }
    }
    int nCanvasCd = nCanvasCols*nCanvasRows;
    //c = new TCanvas("c","Roofit",nCanvasCols*500,nCanvasRows*350);
    c = new TCanvas("c","Roofit",nCanvasCols*800,nCanvasRows*800);
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
            c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
            c->SaveAs((outpdf + toString(canvasPageCounter) + ".png").c_str());
            canvasPageCounter++;
            canvasCdCounter = 1;
            c->Clear("D");
        }
        c->cd(canvasCdCounter);

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

	// some configs
	const int nBins(15*int(curEntries/100));

	// Observables
	RooRealVar t("ct3dlb","measured lifetime #Lambda_{b}",-2e-12, 13e-12);
	RooRealVar dt("ct3dlbE","per-event error on t",0.01e-12,1000e-12);

	// Build a gaussian resolution model scaled by the per-event error = gauss(dt,bias,sigma*dterr)
	const double biasMax(10);
	RooRealVar bias("bias","bias",0,-biasMax,biasMax);
	const double sigmaFactor(4);
	RooRealVar sigma("sigma","per-event error scale factor",1,1/sigmaFactor,1*sigmaFactor);
	RooGaussModel gm("gm1","gauss model scaled by per-event error",t,bias,sigma,dt);

	// Construct decay(dt) (x) gauss1(dt|dterr)
	RooRealVar tau("tau","tau",1.391e-12,.5e-12,5e-12);
	RooDecay decay_gm("decay_gm","decay",t,tau,gm,RooDecay::SingleSided);

	// Getting data from tree
        RooDataSet data("data", "#Lambda_{b} lifetime dataset", RooArgSet(t,dt), RooFit::Import(*subtree));

	// Now do the fit
	decay_gm.fitTo(data, ConditionalObservables(dt));

	t.Print();
	dt.Print();
	tau.Print();
	bias.Print();
	sigma.Print();

	// Plotting
	RooPlot* frame = t.frame(Title("Projection of decay(d|dt) on t"));
	data.plotOn(frame,RooFit::Binning(nBins));
	decay_gm.plotOn(frame,ProjWData(data));
	
	frame->Draw();
    }

    // finalize page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".png").c_str());
    //delete c;
    tmpfile->Write();
}

