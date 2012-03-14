#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TLine.h"
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

TCanvas *c;

struct addtlCuts
{
    addtlCuts(std::string c, std::string n) : cut(c), name(n) {};
    std::string cut;
    std::string name;
};

struct Fitresults
{
    std::vector<double> x;
    std::vector<double> sig, sigE;
    std::vector<double> bgr, bgrE;
    std::vector<double> soversqrtsb, soversqrtsbE;
    std::vector<double> soverb, soverbE;
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

void plotResultErrors(TPad* pad, std::string histoName, const std::vector<double> &x,
	const std::vector<double> &data, const std::vector<double> & errors,
	const std::vector<std::string> &names, std::string title,
	bool forceMin = false, double linepos = -9999)
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
	//h->GetXaxis()->SetNdivisions(505);
    }
    if (forceMin) if (h->GetMinimum() > 0) h->SetMinimum(0);
    h->Draw("P0E");
    if(linepos!=-9999)
    {
	cout << "Draw a line at " << linepos << endl;
	cout << "Draw " << h->GetXaxis()->GetXmin() << " " << h->GetXaxis()->GetXmax() << endl;
	// draw a vertical line at the given x position
	const double fractionOfHeight = .96;
	const double curYMax = h->GetMaximum();
	const double curYMin = h->GetMinimum();
	const double height = curYMax-curYMin ;
	const double lineMin = curYMin + (1-fractionOfHeight) * height;
	const double lineMax = curYMax - (1-fractionOfHeight) * height;
	TLine *l = new TLine(linepos,lineMin,linepos,lineMax);
	l->SetLineColor(2);
	l->Draw();
    }
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

void rooFitCutScan01 (std::string file1, std::string file2, std::string cutToVary, int nVariations, double lo, double hi, bool doPublicationGrade = false)
{
    const std::string outpdf = "rooFitCutScan";

    const bool noTitle(doPublicationGrade);

    // container for the results
    Fitresults fitresults;

    // Select cuts from cuts.C
    Cuts cut;
    cut.selectCut("an04","acc03_16","HLT_matched_01");

    cout << "This is rooFitCutScan01 using the following cuts as a basis:" << endl;
    cout << cut.getCut() << endl;
    cout << "----- BEGIN value block -----" << endl;
    cut.printCutLaTeX_values("cuts:anal:");
    cout << "----- END value block -----" << endl;

    // Open file
    TFile* file = TFile::Open(file1.c_str());
    if (file ==0)
    {
        cout << "File " << file1 << " not found. Exiting." << endl;
        return;
    }
    cout << "File " << file1 << " succesfully opened." << endl;

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
    const int maxCanvasCols = doPublicationGrade ? 1 : 4;
    const int maxCanvasRows = doPublicationGrade ? 1 : 5;
    const int nAdditionalPlots = 5;
    const int nPlots = nVariations + nAdditionalPlots;
    int nCanvasCols = maxCanvasCols;
    int nCanvasRows = maxCanvasRows;
    if (nPlots < maxCanvasCols*maxCanvasRows)
    {
        if(nPlots == 1)
            nCanvasCols = nCanvasCols = 1;
        else
        {
            nCanvasCols = maxCanvasCols;
            nCanvasRows = TMath::Ceil((double)nPlots / nCanvasCols);
        }
    }
    int nCanvasCd = nCanvasCols*nCanvasRows;
    c = new TCanvas("c","Roofit",nCanvasCols*500,nCanvasRows*350);
    c->Divide(nCanvasCols,nCanvasRows);

    // now create a temporary file to speed up the process later
    TFile* tmpfile = new TFile("tmp.root","RECREATE");
    if (tmpfile == 0)
    {
	cout << "temporary file creation failed - exiting" << endl;
	return;
    }
    cout << "tmpfile: " << tmpfile << endl;
    TTree* treepresel = 0;

    cout << "Preselecting events... Excluding " << cutToVary << endl;
    treepresel = tree->CopyTree(cut.getCutExceptOne(cutToVary).c_str());
    cout << "Cutstring is: " << cut.getCutExceptOne(cutToVary) << endl;
    cout << "now " << treepresel->GetEntries() << " of " << tree->GetEntries()<< endl;
    delete tree;

    // now comes the main loop
    int canvasCdCounter(0), canvasPageCounter(0);
    for (unsigned int i=0; i!=nVariations+1; i++)
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
	const double curCutVal = lo+i*(hi-lo)/nVariations;
	const std::string curCut = cut.getOneCut(cutToVary,curCutVal);
	const std::string curTitle = cutToVary + ": " + toString(curCutVal);
	cout << curTitle << " - " << curCut << endl;

        TTree* subtree;
        subtree = treepresel->CopyTree(curCut.c_str());
        cout << "subtree: " << subtree->GetEntries() << " entries of " << treepresel->GetEntries()<< endl;
	const double curEntries = subtree->GetEntries();
	cout << "curEntries: " << curEntries << endl;

	// some configs
	const double valLo(5.22), valHi(6.045);
	//const double valLo(5.35), valHi(5.90);
	const int nBins(33);
	//const int nBins(22);

        RooRealVar mass("mlb", "J/#psi #Lambda mass [GeV/c^{2}]", valLo, valHi);
        RooDataSet data("data", "Lambda_b dataset", subtree, mass);

        //RooRealVar lb_lifetime("lb_lifetime", "c#tau #Lambda_{b} [mm]",0.25, 3.5);
        //RooDataSet lt_data("lt_data", "#Lambda_{b} lifetime dataset", tree, lb_lifetime);

        //   mass fit
        RooRealVar mean("mean","mean",5.62,5.60,5.64) ;
        RooRealVar sigma("sigma","sigma",0.014,0.01,0.03) ;
        RooGaussian sig("sig","signal p.d.f.",mass,mean,sigma) ;

        RooRealVar c0("c0","coefficient #0", 1.0,-1,2) ;
        RooRealVar c1("c1","coefficient #1", 0.1,-1,1) ;
        //RooPolynomial bkg("bkg","background p.d.f.",mass,RooArgList(c0,c1),0) ;
        RooChebychev bkg("bkg","background p.d.f.",mass,RooArgSet(c0)) ;

        RooRealVar nsig("nsig","signal fraction", .5*curEntries,0.,curEntries) ;
        RooRealVar nbkg("nbkg","Background fraction", .01* curEntries,0.,curEntries*1.4) ;

        RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg)) ;
        model.fitTo(data,Extended(kTRUE));

        RooPlot* xframe = mass.frame(Name("xframe"),Title(("Mass of #Lambda_{b} - "+ curTitle).c_str())) ;
        if(noTitle) xframe->SetTitle("");
        data.plotOn(xframe,Binning(nBins)) ;
        model.plotOn(xframe);
        model.plotOn(xframe,Components("bkg"),LineStyle(kDashed));
	
	xframe->SetYTitle("Events / (0.025 GeV/c^{2})");
        xframe->Draw();

        Double_t m=mean.getVal();
	Double_t m_err=mean.getError();
        Double_t s1=sigma.getVal();
        Double_t s1_err=sigma.getError();

        mass.setRange("window",m-2.5*s1,m+2.5*s1);
        RooAbsReal* fracSigRange = sig.createIntegral(mass,mass,"window");
        Double_t nsigWindow = fracSigRange->getVal()*nsig.getVal();
	Double_t nsig_err = nsig.getError()*fracSigRange->getVal();
        RooAbsReal* fracBGRange = bkg.createIntegral(mass,mass,"window");
        Double_t nbkgWindow = nbkg.getVal()*fracBGRange->getVal();
	Double_t nbkg_err = nbkg.getError()*fracBGRange->getVal();

        cout << "n_Signal     =  "<<nsigWindow<<endl;
        cout << "n_Background =  "<<nbkgWindow<<endl;
        cout << "S/Sqrt(S+B)  =  "<<nsigWindow/sqrt(nsigWindow+nbkgWindow)<<endl;

        // write the results to the canvas
        TLatex *txt;
        const double txtPosLeft = .6;
        const double txtPosTop = .80;
        const double txtLineSpace = .05;

	double sb_err = sqrt(nsig_err*nsig_err/(nsigWindow*nsigWindow)+nbkg_err*nbkg_err/(nbkgWindow*nbkgWindow))*nsigWindow/nbkgWindow;
	double s_over_sqrt_err = sqrt(nsig_err*nsig_err/(nsigWindow*nsigWindow)+(nbkg_err*nbkg_err+nsig_err*nsig_err)/(4.0*(nbkgWindow+nsigWindow)*(nbkgWindow+nsigWindow)))*nsigWindow/sqrt(nsigWindow+nbkgWindow);
        xframe->addObject(writeTLatex("n Signal: " + roundToString(nsigWindow,1)+" #pm "+roundToString(nsig_err,1),txtPosLeft,txtPosTop-0*txtLineSpace));
        xframe->addObject(writeTLatex("n Bgr: " + roundToString(nbkgWindow,1)+" #pm "+roundToString(nbkg_err,1), txtPosLeft,txtPosTop-1*txtLineSpace));
        xframe->addObject(writeTLatex("S/#sqrt{S+B}: " + roundToString(nsigWindow/sqrt(nsigWindow+nbkgWindow),1) + " #pm " + roundToString(s_over_sqrt_err,1) ,txtPosLeft,txtPosTop-2*txtLineSpace));
        xframe->addObject(writeTLatex("S/B: " + roundToString(nsigWindow/nbkgWindow,2) + " #pm " + roundToString(sb_err,2),txtPosLeft,txtPosTop-3*txtLineSpace));
        xframe->addObject(writeTLatex("mass: " + roundToString(m,3) + " #pm " + roundToString(m_err,3)+" GeV/c^{2}", txtPosLeft,txtPosTop-4*txtLineSpace));
        xframe->addObject(writeTLatex("width: " + roundToString(s1,3) + " #pm " + roundToString(s1_err,3)+ " GeV/c^{2}",txtPosLeft,txtPosTop-5*txtLineSpace));
        xframe->Draw();
        //c->Update();

        // collect results
	fitresults.x.push_back(curCutVal);
        fitresults.sig.push_back(nsigWindow);
        fitresults.sigE.push_back(nsig_err);
        fitresults.bgr.push_back(nbkgWindow);
        fitresults.bgrE.push_back(nbkg_err);
        fitresults.soversqrtsb.push_back(nsigWindow/sqrt(nsigWindow+nbkgWindow));
        fitresults.soversqrtsbE.push_back(s_over_sqrt_err);
        fitresults.soverb.push_back(nsigWindow/nbkgWindow);
        fitresults.soverbE.push_back(sb_err);
        fitresults.mass.push_back(m);
        fitresults.width.push_back(s1);
        fitresults.title.push_back(toString(curCutVal));
    }

    // draw result histos
    for (unsigned int i=0; i!= nAdditionalPlots; i++)
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
        TPad *pad = (TPad*)c->cd(canvasCdCounter);
        if (0 == i) plotResultErrors(pad, "hsig", fitresults.x, fitresults.sig, fitresults.sigE, fitresults.title, "Signal", true, cut.getOneCutValue(cutToVary));
        if (1 == i) plotResultErrors(pad, "hbgr", fitresults.x, fitresults.bgr, fitresults.bgrE, fitresults.title, "Background", true);
        if (2 == i) plotResultErrors(pad, "hssqtsb", fitresults.x, fitresults.soversqrtsb, fitresults.soversqrtsbE, fitresults.title, "S/#sqrt(S+B)", true);
        if (3 == i) plotResultErrors(pad, "hsoverb", fitresults.x, fitresults.soverb, fitresults.soverbE, fitresults.title, "S/B", true);
        if (4 == i) plotResultErrors(pad, "hmass", fitresults.x, fitresults.mass, fitresults.width, fitresults.title, "Mass");
    }

    // finalize page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".png").c_str());
    delete c;

    delete treepresel;
    cout << "this Cut: " << cutToVary << " " << cut.getOneCutValue(cutToVary)
	 << " nBins: " << fitresults.x.size()
	 << " range: " << fitresults.x[0] << " " << fitresults.x[fitresults.x.size()-1] << endl;
    tmpfile->Write();
}

