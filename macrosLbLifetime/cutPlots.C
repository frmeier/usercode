#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLine.h"
#include "TPaveText.h"
#include <string>
#include "TLatex.h"
#include <sstream>
#include <iomanip>
#include <map>
#include <string>
#include "setTDRStyle_modified.C"

//#include "cut.C"
#include "cuts.C" // implies cut.C
#include "canvaspager.h"

using std::cout;
using std::endl;

TCanvas *c;

std::map<std::string,std::string> initMap()
{
    std::map<std::string,std::string> namemap;
    namemap["ml0"] = "m(#Lambda^{0})";
    namemap["mjp"] = "m(J/#psi)";
    namemap["ptjp"] = "p_{T}(J/#psi)";
    namemap["rptpr"] = "p_{T}(p)";
    namemap["rptpi"] = "p_{T}(#pi)";
    namemap["rpt1m"] = "p_{T}(#mu_{1})";
    namemap["rpt2m"] = "p_{T}(#mu_{2})";
    namemap["ptl0"] = "p_{T}(#Lambda^{0})";
    namemap["probjp"] = "#chi^{2} probability of J/#psi vertex";
    namemap["probl0"] = "#chi^{2} probability of #Lambda^{0} vertex";
    namemap["problb"] = "#chi^{2} probability of #Lambda_{b} vertex";
    namemap["probpr"] = "#chi^{2} probability of p track";
    namemap["probpi"] = "#chi^{2} probability of #pi track";
    namemap["Kshypo"] = "veto on m(K_{s}) assuming #pi#pi";
    namemap["alphal0"] = "#alpha(#Lambda^{0})";
    namemap["alphalb"] = "#alpha(#Lambda_{b})";
    namemap["d3l0"] = "flight length of #Lambda^{0}";
    namemap["d3l0/d3El0"] = "significance of #Lambda^{0} vertex";
    namemap["ptgangDRl0"] = "Pointing angle of #Lambda^{0} as #DeltaR";
    return namemap;
}


struct addtlCuts
{
    addtlCuts(std::string c, std::string n) : cut(c), name(n) {};
    std::string cut;
    std::string name;
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

void setPadParams(TPad *pad)
{
    pad->SetTopMargin(0.1);
    pad->SetBottomMargin(0.10);
    pad->SetLeftMargin(0.10);
    pad->SetRightMargin(0.10);
}

void cutPlots(std::string filename, std::string cutname, double factor = 4, bool doPres = false)
{
    const std::string outpdf = "cutPlots";
    const std::string drawstring = "mlb";
    const std::string histoFormat = "(20,5,6)";
    const bool doPlotExcl(false);
    const bool doPlotInv(false);
    const bool doPlotDistr(true);

    std::map<std::string,std::string> namemap = initMap();

    // Select cuts from cuts.C
    Cuts cut;
    //cut.selectCut(cutname);
    cut.selectCut(cutname,"acc03","HLT_matched_2011_02");

    cout << "This is cutPlots using the following cuts as a basis:" << endl;
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

    // initialise the TCanvas
    setTDRStyle();
    const int maxCanvasCols = doPres ? 1 : 3;
    const int maxCanvasRows = doPres ? 1 : 5;
    const int nAdditionalPlotRows = 1;
    const int nPlotRows = cut.cs.size();
    int nCanvasCols = maxCanvasCols;
    int nCanvasRows = nPlotRows > maxCanvasRows ? maxCanvasRows : nPlotRows;
    int nCanvasCd = nCanvasCols*nCanvasRows;
    const int nPads = nCanvasCd;
    c = new TCanvas("c","cutPlots",nCanvasCols*500,nCanvasRows*350);
    c->Divide(nCanvasCols,nCanvasRows);

    int canvasCdCounter(0), canvasPageCounter(0);
    for (unsigned int i=0; i!=cut.cs.size(); i++)
    {
	if(doPlotExcl)
	{
	    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	    TPad *pad = (TPad*)c->cd(canvasCdCounter);
	    const std::string curCutExceptOne = cut.cs.getCutExceptOne(cut.parvec, i);
	    const std::string curCutName = cut.getOneCut(i);
	    const std::string curHistName = "hexc" + toString(i);
	    const std::string curdrawstring = drawstring + ">>" + curHistName + histoFormat;
	    //cout << i << ": " << curCutExceptOne.c_str() << endl;
	    tree->Draw(curdrawstring.c_str(),curCutExceptOne.c_str());
	    TH1F *h = (TH1F*)gDirectory->Get(curHistName.c_str());
	    h->SetTitle((curCutName + " excluded").c_str());
	    h->Draw();
	}

	if (doPlotInv)
	{
	    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	    c->cd(canvasCdCounter);
	    const std::string curCutOneInverted = cut.cs.getCutInvertN(cut.parvec, i);
	    const std::string curCutName = cut.getOneCut(i);
	    const std::string curHistName = "hinv" + toString(i);
	    const std::string curdrawstring = drawstring + ">>" + curHistName + histoFormat;
	    //cout << i << ": " << curCutOneInverted.c_str() << endl;
	    tree->Draw(curdrawstring.c_str(),curCutOneInverted.c_str());
	    TH1F *h = (TH1F*)gDirectory->Get(curHistName.c_str());
	    h->SetTitle((curCutName + " inverted").c_str());
	    h->Draw();
	}


	if (doPlotDistr)
	{
	    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	    TPad* pad = (TPad*)c->cd(canvasCdCounter);
	    setPadParams(pad);
	    const std::string curCutExceptOne = cut.cs.getCutExceptOne(cut.parvec, i);
	    const std::string curVarName = cut.cs.getCutName(i);
	    const std::string curHistName = "hvar" + toString(i);
	    std::string curHistoFormat = "(25,0," + toString(factor*cut.parvec[i]) + ")";
	    const std::string curCutClassName = cut.cs.getCutClassName(i);
	    if(curCutClassName=="cutSymWindow" || curCutClassName=="cutSymWindowVeto") 
	    {
		const double val = cut.parvec[i];
		const double lo = cut.cs.getOneCutValue(i,-val*factor);
		const double hi = cut.cs.getOneCutValue(i,+val*factor);
		curHistoFormat = "(25," + toString(lo) + "," + toString(hi) + ")";
	    }
	    const std::string curdrawstring = curVarName + ">>" + curHistName + curHistoFormat;
	    cout << curdrawstring << endl;
	    tree->Draw(curdrawstring.c_str(),curCutExceptOne.c_str());
	    TH1F *h = (TH1F*)gDirectory->Get(curHistName.c_str());
	    h->SetTitle("");
	    h->Draw();
	    TPaveText *ptTitle = new TPaveText(0.25,0.91,0.80,0.96,"LNDC");
	    if(namemap.find(curVarName.c_str())!=namemap.end())
		ptTitle->AddText(namemap[curVarName].c_str());
	    else
		ptTitle->AddText(curVarName.c_str());
	    ptTitle->SetBorderSize(0);
	    ptTitle->SetFillColor(0);
	    ptTitle->SetTextSize(0.05);
	    ptTitle->Draw();
	    if(curCutClassName=="cutSymWindow" || curCutClassName=="cutSymWindowVeto") 
	    {
		const double val = cut.parvec[i];
		const double lo = cut.cs.getOneCutValue(i,-val);
		const double hi = cut.cs.getOneCutValue(i,+val);
		TLine *l1 = new TLine(lo,0,lo,h->GetMaximum());
		l1->SetLineColor(2);
		l1->Draw();
		TLine *l2 = new TLine(hi,0,hi,h->GetMaximum());
		l2->SetLineColor(2);
		l2->Draw();
	    }
	    else
	    {
		TLine *l = new TLine(cut.parvec[i],0,cut.parvec[i],h->GetMaximum());
		l->SetLineColor(2);
		l->Draw();
	    }
	}
    }

    // finalize page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
    delete c;
}

