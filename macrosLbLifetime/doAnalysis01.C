#include <iostream>
#include <memory>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TDirectory.h"
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
#include <vector>
#include "TLatex.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include "TSystem.h"
#include "TEventList.h"

#include "setTDRStyle_modified.C"

#include "Canvaspager.h"
#include "doMassFitLb01.C"
#include "doMassFitB001.C"
#include "doMassFitJp01.C"
#include "doSidebandPlot.C"
#include "TriggerPlot.C"
#include "do2dPlot.C"
#include "do1dPlot.C"
#include "do1dPlotGaus.C"
#include "do2dProfilePlot.C"
#include "HtmlReport.C"
#include "JsonRunList.h"
#include "doEfficiencyPlotFitterLb.C"
#include "doEfficiencyPlotFitterB0.C"
#include "doEfficiencyPlotFitterJp.C"
#include "doCtFitB0_01.C"
#include "doStackedPlot.C"

#include "Cut.h"
#include "Cuts.h" // implies cut.C
#include "Datafile.h"

using std::cout;
using std::endl;
using RooFit::Extended;
using RooFit::Name;
using RooFit::Title;
using RooFit::Binning;
using RooFit::Components;
using RooFit::LineStyle;

TCanvas *canvas;
HtmlReport *htrep;

bool flgNoTitle, flgPrelim, flgOfficialPlot, flgDoHtmlReport;

// ===================================================================================
void writePickFile(TTree *t, string cut, string filename)
{
    t->Draw(">>lst", cut.c_str());
    TEventList *lst;
    lst = (TEventList*)gDirectory->Get("lst");
    int run, ls;
    unsigned int event;
    t->SetBranchAddress("run",&run);
    t->SetBranchAddress("LS",&ls);
    t->SetBranchAddress("event",&event);
    std::ofstream ofs(filename.c_str());
    for(int i = 0; i!=lst->GetN(); i++)
    {
	t->GetEntry(lst->GetEntry(i));
	ofs << run << ":" << ls << ":" << event << endl;
    }
    ofs.close();
    cout << "Wrote " << lst->GetN() << " events in " << filename << endl;
}

// ===================================================================================
void writeJsonFile(TTree *t, string cut, string filename)
{
    JsonRunList jsnlist;
    t->Draw(">>lst", cut.c_str());
    TEventList *lst;
    lst = (TEventList*)gDirectory->Get("lst");
    int run, ls;
    t->SetBranchAddress("run",&run);
    t->SetBranchAddress("LS",&ls);
    for(int i = 0; i!=lst->GetN(); i++)
    {
	t->GetEntry(lst->GetEntry(i));
	jsnlist.insert(run,ls);
    }
    delete lst;
    std::ofstream ofs(filename.c_str());
    ofs << jsnlist;
    ofs.close();
    cout << "Wrote " << jsnlist.getSize() << " LS in " << filename << endl;
}

// ===================================================================================
void massPlotsPtYbind(TTree *tree, Canvaspager &cp, Cuts curCut, string addTitle, double lumi)
{
    if (flgDoHtmlReport) htrep->addTableCell("Binned in p<sub>T</sub>:");
    std::auto_ptr<TTree> subtree;

    // do some bin selection
    cp.cdNext();
    Cuts cutPt_10_15;
    cutPt_10_15.selectCut("ptlb_10_15");
    subtree.reset(tree->CopyTree((curCut+cutPt_10_15).getCut().c_str()));
    doMassFitLb01_fitresults resMassFitLb_ptlb_10_15 = doMassFitLb01(subtree.get(), lumi, "10<p_{T}(#Lambda_{b})<15 " + addTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"10 &lt; p<sub>T</sub> &lt 15");

    cp.cdNext();
    Cuts cutPt_15_20;
    cutPt_15_20.selectCut("ptlb_15_20");
    subtree.reset(tree->CopyTree((curCut+cutPt_15_20).getCut().c_str()));
    doMassFitLb01_fitresults resMassFitLb_ptlb_15_20 = doMassFitLb01(subtree.get(), lumi, "15<p_{T}(#Lambda_{b})<20 " + addTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"15 &lt; p<sub>T</sub> &lt 20");

    cp.cdNext();
    Cuts cutPt_20_infty;
    cutPt_20_infty.selectCut("ptlb_20_infty");
    subtree.reset(tree->CopyTree((curCut+cutPt_20_infty).getCut().c_str()));
    doMassFitLb01_fitresults resMassFitLb_ptlb_20_infty = doMassFitLb01(subtree.get(), lumi, "20<p_{T}(#Lambda_{b}) " + addTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"20 &lt; p<sub>T</sub>");

    if (flgDoHtmlReport) htrep->addTableCell("Binned in rapidity <em>y</em>:");
    cp.cdNext();
    Cuts cutY_00_05;
    cutY_00_05.selectCut("ylb_00_05");
    subtree.reset(tree->CopyTree((curCut+cutY_00_05).getCut().c_str()));
    doMassFitLb01_fitresults resMassFitLb_ylb_00_05 = doMassFitLb01(subtree.get(), lumi, "|y(#Lambda_{b})|<0.5 " + addTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"|y| &lt; 0.5");

    cp.cdNext();
    Cuts cutY_05_12;
    cutY_05_12.selectCut("ylb_05_12");
    subtree.reset(tree->CopyTree((curCut+cutY_05_12).getCut().c_str()));
    doMassFitLb01_fitresults resMassFitLb_ylb_05_12 = doMassFitLb01(subtree.get(), lumi, "0.5#leq|y(#Lambda_{b})|<1.2 " + addTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"0.5 &lt; |y| &lt; 1.2");

    cp.cdNext();
    Cuts cutY_12_22;
    cutY_12_22.selectCut("ylb_12_22");
    subtree.reset(tree->CopyTree((curCut+cutY_12_22).getCut().c_str()));
    doMassFitLb01_fitresults resMassFitLb_ylb_12_22 = doMassFitLb01(subtree.get(), lumi, "1.2#leq|y(#Lambda_{b})|<2.2 " + addTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"1.2 &lt; |y| &lt; 2.2");

    if (flgDoHtmlReport)
    {
	htrep->beginTable();
	htrep->addTableRow("Bin","# sig","# bgr");
	htrep->addTableRow("10 &lt; p<sub>T</sub> &lt 15",roundToString(resMassFitLb_ptlb_10_15.sig,1),roundToString(resMassFitLb_ptlb_10_15.bgr,1));
	htrep->addTableRow("15 &lt; p<sub>T</sub> &lt 20",roundToString(resMassFitLb_ptlb_15_20.sig,1),roundToString(resMassFitLb_ptlb_15_20.bgr,1));
	htrep->addTableRow("20 &lt; p<sub>T</sub>",roundToString(resMassFitLb_ptlb_20_infty.sig,1),roundToString(resMassFitLb_ptlb_20_infty.bgr,1));
	htrep->addTableRow("|y|<0.5 ",roundToString(resMassFitLb_ylb_00_05.sig,1),roundToString(resMassFitLb_ylb_00_05.bgr,1));
	htrep->addTableRow("0.5 &lt; |y| &lt; 1.2",roundToString(resMassFitLb_ylb_05_12.sig,1),roundToString(resMassFitLb_ylb_05_12.bgr,1));
	htrep->addTableRow("1.2 &lt; |y| &lt; 2.2",roundToString(resMassFitLb_ylb_12_22.sig,1),roundToString(resMassFitLb_ylb_12_22.bgr,1));
    }
}

// ===================================================================================
void cutsPlotsSidebandLb(TTree *tree, Canvaspager &cp, Cuts curCut, string addTitle, string addHname, doMassFitLb01_fitresults resMassFit)
{
    if (flgDoHtmlReport) htrep->addP(addTitle + ":", true);
    cout << "=======================================" << endl;
    cout << "Doing cut parameter plots for Lb.... " << addTitle << endl;

    const double sigWindow(2.0), sideWindowClose(3.0), sideWindowFar(6.0);
    const double mass = resMassFit.mass;
    const double width = resMassFit.width;
    string sigCut = "mlb>" + toString(mass-sigWindow*width) + "&&mlb<" + toString(mass+sigWindow*width);
    string sideCut = "(mlb>=" + toString(mass+sideWindowClose*width) + "&& mlb <" + toString(mass+sideWindowFar*width) + ")";
    //string sideCut = "((mlb>" + toString(mass-sideWindowFar*width) + "&& mlb <=" + toString(mass-sideWindowClose*width) + ")||" +
	//		"(mlb>=" + toString(mass+sideWindowClose*width) + "&& mlb <" + toString(mass+sideWindowFar*width) + "))";
    //string sideCut = "(mlb<=" + toString(mass-sideWindow*width) + "||mlb>=" + toString(mass+sigWindow*width) + ")";

    /*
    cp.cdNext();
    doSidebandPlot(tree, "ct3dlb_all"+addHname, "ct3dlb", "L_{b} lifetime 3d "+addTitle, "#tau_{3d}(#Lambda_{b})", "s", 20, 0, 9e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ct3dlb_sig"+addHname, "ct3dlb/ct3dlbE", "L_{b} lifetime significance 3d "+addTitle, "#tau_{3d}/#varepsilon (#Lambda_{b})", "", 20, 0, 10, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlotLogY(tree, "ct3dlb_alllog"+addHname, "ct3dlb", "L_{b} lifetime 3d (log)  "+addTitle, "#tau_{3d}(#Lambda_{b})", "s", 30, 0, 9e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ct3dlb_red"+addHname, "ct3dlb-3*ct3dlbE", "L_{b} reduced lifetime 3d "+addTitle, "#tau_{3d}(#Lambda_{b})", "s", 20, 0, 9e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ct3dlbE"+addHname, "ct3dlbE", "L_{b} lifetime uncertainty 3d "+addTitle, "#varepsilon_{#tau_{3d}}(#Lambda_{b})", "s", 30, 0, 3e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ctxylb_all"+addHname, "ctxylb", "L_{b} lifetime 2d "+addTitle, "#tau_{xy}(#Lambda_{b})", "s", 20, 0, 9e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ctxylb_sig"+addHname, "ctxylb/ctxylbE", "L_{b} lifetime 2d significance "+addTitle, "#tau_{xy}/#varepsilon (#Lambda_{b})", "", 20, 0, 10, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ctxylbE"+addHname, "ctxylbE", "L_{b} lifetime 2d uncertainty "+addTitle, "#varepsilon_{#tau_{xy}}(#Lambda_{b})", "s", 60, 0, 3e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "d3lb_all"+addHname, "d3lb", "#Lambda_{b} flightlength 3d"+addTitle, "d_{3d}(#Lambda_{b})", "cm", 20, 0, 1, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlotLogY(tree, "d3lb_alllog"+addHname, "d3lb", "#Lambda_{b} flightlength 3d (log) "+addTitle, "d_{3d}(#Lambda_{b})", "cm", 20, 0, 1, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "d3Elb_all"+addHname, "d3Elb", "#Lambda_{b} flightlength uncertainty 3d "+addTitle, "#varepsilon_{d_{3d}(#Lambda_{b})}", "cm", 40, 0, 0.3, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "dxyElb_all"+addHname, "dxyElb", "#Lambda_{b} flightlength uncertainty xy "+addTitle, "#varepsilon_{d_{xy}(#Lambda_{b})}", "cm", 40, 0, 0.3, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "d3lb_sig"+addHname, "d3lb/d3Elb", "#Lambda_{b} flightlength significance 3d "+addTitle, "d_{3d}(#Lambda_{b}) significance", "", 20, 0, 5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "mjp"+addHname, "mjp", "Jpsi mass "+addTitle, "m(J/#psi)", "GeV/c^{2}", 20, 2.9, 3.3, curCut.getCutExceptOne("mjp"), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ml0"+addHname, "ml0", "#Lambda^{0} mass "+addTitle, "m(#Lambda^{0})", "GeV/c^{2}", 20, 1.091, 1.141, curCut.getCutExceptOne("ml0"), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "maxdocalb"+addHname, "maxdocalb", "max doca Lb "+addTitle, "max d_{ca}(#Lambda_{b})", "cm", 25, 0, .25, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ptlb"+addHname, "ptlb", "pt lb "+addTitle, "p_{T}(#Lambda_{b})", "GeV/c", 25, 0, 50, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "etalb"+addHname, "TMath::Abs(etalb)", "eta lb "+addTitle, "|#eta(#Lambda_{b})|", "", 20, 0, 2.5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ptgangDRl0"+addHname, "ptgangDRl0", "pointing angle l0 as dR "+addTitle, "dR_{pointing}(#Lambda^{0})", "", 50, 0, .1, curCut.getCutExceptOne("ptgangDRl0"), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "dRmumu"+addHname, "dRmumu", "dR between muons "+addTitle, "dR_(#mu#mu})", "", 50, 0, 1, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");
    */

    cp.cdNext();
    doSidebandPlot(tree, "PvLip"+addHname, "PvLip", "longit I.P. to best PV "+addTitle, "PV lip", "cm", 25, 0, .2, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip2"+addHname, "PvLip2", "longit I.P. to 2nd best PV "+addTitle, "PV lip", "cm", 25, 0, 2, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLipSig"+addHname, "PvLip/PvLipE", "longit I.P. to best PV "+addTitle, "PV lip significance", "", 25, 0, 5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip2Sig"+addHname, "PvLip2/PvLipE2", "longit I.P. to 2nd best PV "+addTitle, "PV lip significance", "", 25, 0, 5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip2Sig_2"+addHname, "PvLip2/PvLipE2", "longit I.P. to 2nd best PV "+addTitle, "PV lip significance", "", 25, 0, 50, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip2Sig_inv"+addHname, "PvLipE2/PvLip2", "longit I.P. to 2nd best PV "+addTitle, "PV lip inv. significance", "", 25, 0, 1, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "rpt1m"+addHname, "rpt1m", "p_{T}(#mu,1) "+addTitle, "p_{T}(#mu,1)", "GeV/c", 50, 0, 20, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "rpt2m"+addHname, "rpt2m", "p_{T}(#mu,2) "+addTitle, "p_{T}(#mu,2)", "GeV/c", 50, 0, 20, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "rptpi"+addHname, "rptpi", "p_{T}(#pi) "+addTitle, "p_{T}(#pi)", "GeV/c", 80, 0, 5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "rptpr"+addHname, "rptpr", "p_{T}(p) "+addTitle, "p_{T}(p)", "GeV/c", 80, 0, 20, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ptl0"+addHname, "ptl0", "p_{T}(#Lambda^{0}) "+addTitle, "p_{T}(#Lambda^{0})", "GeV/c", 60, 0, 30, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ptjp"+addHname, "ptjp", "p_{T}(J/#psi) "+addTitle, "p_{T}(J/#psi)", "GeV/c", 60, 0, 30, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");
}

void cutPlotSidebandDataMCLb(TTree *datatree, TTree *mctree, Canvaspager &cp, std::string varToPlot, bool exclCut, int nBins, Double_t lo, Double_t hi, Cuts cutData, Cuts cutMC, string addTitle, string addHname, string axisName, string axisUnit, const doMassFitLb01_fitresults &fitres)
{
    const std::string cutstringData = exclCut ? cutData.getCutExceptOne(varToPlot) : cutData.getCut();
    const std::string cutstringMC   = exclCut ? cutMC.getCutExceptOne(varToPlot)   : cutMC.getCut();

    const double mass = fitres.mass;
    const double twoSigmas = fitres.twoSigmas;
    const double threeSigmas = fitres.threeSigmas;

    string sigCut = "mbc>" + toString(mass-twoSigmas) + "&&mbc<" + toString(mass+twoSigmas);
    string sideCut = "(mbc>=" + toString(mass+threeSigmas) + "&&mbc<" + toString(6) + ")";

    cp.cdNext("Lb_sidebandsubtracted_"+addHname);

    // get cut lines
    double drawLine1(-9999), drawLine2(-9999);
    if (exclCut)
    {
	cutSet::sizeType cutno = cutData.cs.getCutPos(varToPlot);
	const string curCutClassName = cutData.cs.getCutClassName(cutno);

	if(curCutClassName=="cutSymWindow")
	{
	    const double val = cutData.parvec[cutno];
	    drawLine1 = cutData.cs.getOneCutValue(cutno,-val);
	    drawLine2 = cutData.cs.getOneCutValue(cutno,+val);
	}
	else if(curCutClassName=="cutSymWindowVeto")
	{
	    const double val = cutData.parvec[cutno];
	    drawLine1 = cutData.cs.getOneCutValue(cutno,+val);
	    drawLine2 = cutData.cs.getOneCutValue(cutno,-val);
	}
	else if(curCutClassName=="cutBoundUpper")
	{
	    drawLine2 = cutData.parvec[cutno];
	}
	else
	{
	    drawLine1 = cutData.parvec[cutno];
	}
    }

    doSidebandPlot(datatree, mctree, addHname, varToPlot, addTitle, axisName, axisUnit, nBins, lo, hi, cutstringData, cutstringMC, sigCut, sideCut, fitres.sig, fitres.bgr, drawLine1, drawLine2);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(), "", cp.getCurPdf());
    if (flgDoHtmlReport) htrep->endTableRow();
}

doMassFitLb01_fitresults cutPlotSidebandDataMCLb(TTree *datatree, TTree *mctree, Canvaspager &cp, std::string varToPlot, bool exclCut, int nBins, Double_t lo, Double_t hi, Cuts cutData, Cuts cutMC, string addTitle, string addHname, string axisName, string axisUnit)
{
    const std::string cutstringData = exclCut ? cutData.getCutExceptOne(varToPlot) : cutData.getCut();

    std::auto_ptr<TTree> subtree(datatree->CopyTree(cutstringData.c_str()));
    cp.cdNext("B0_massplot4sidebandsubtracted_"+addHname);
    doMassFitLb01_fitresults fitres = doMassFitLb01_DG(subtree.get(), 0, addTitle, flgNoTitle, flgPrelim, flgOfficialPlot, true);
    const double mass = fitres.mass;
    const double twoSigmas = fitres.twoSigmas;
    const double threeSigmas = fitres.threeSigmas;
    const string annotation = "2&sigma; range: " + toString(mass-twoSigmas) + "-" + toString(mass+twoSigmas)
	+ " 3&sigma; range: " + toString(mass-threeSigmas) + "-" + toString(mass+threeSigmas);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(), annotation, cp.getCurPdf());

    cutPlotSidebandDataMCLb(datatree, mctree, cp, varToPlot, exclCut, nBins, lo, hi, cutData, cutMC, addTitle, addHname, axisName, axisUnit, fitres);
    return fitres;
}

void cutPlotSidebandDataMCB0(TTree *datatree, TTree *mctree, Canvaspager &cp, std::string varToPlot, bool exclCut, int nBins, Double_t lo, Double_t hi, Cuts cutData, Cuts cutMC, string addTitle, string addHname, string axisName, string axisUnit, const doMassFitB001_fitresults &fitres)
{
    const std::string cutstringData = exclCut ? cutData.getCutExceptOne(varToPlot) : cutData.getCut();
    const std::string cutstringMC   = exclCut ? cutMC.getCutExceptOne(varToPlot)   : cutMC.getCut();

    const double mass = fitres.mass;
    const double twoSigmas = fitres.twoSigmas;
    const double threeSigmas = fitres.threeSigmas;

    string sigCut = "mbc>" + toString(mass-twoSigmas) + "&&mbc<" + toString(mass+twoSigmas);
    string sideCut = "(mbc>=" + toString(mass+threeSigmas) + "&&mbc<" + toString(6) + ")";

    cp.cdNext("B0_sidebandsubtracted_"+addHname);

    // get cut lines
    double drawLine1(-9999), drawLine2(-9999);
    if (exclCut)
    {
	cutSet::sizeType cutno = cutData.cs.getCutPos(varToPlot);
	const string curCutClassName = cutData.cs.getCutClassName(cutno);

	if(curCutClassName=="cutSymWindow")
	{
	    const double val = cutData.parvec[cutno];
	    drawLine1 = cutData.cs.getOneCutValue(cutno,-val);
	    drawLine2 = cutData.cs.getOneCutValue(cutno,+val);
	}
	else if(curCutClassName=="cutSymWindowVeto")
	{
	    const double val = cutData.parvec[cutno];
	    drawLine1 = cutData.cs.getOneCutValue(cutno,+val);
	    drawLine2 = cutData.cs.getOneCutValue(cutno,-val);
	}
	else if(curCutClassName=="cutBoundUpper")
	{
	    drawLine2 = cutData.parvec[cutno];
	}
	else
	{
	    drawLine1 = cutData.parvec[cutno];
	}
    }

    doSidebandPlot(datatree, mctree, addHname, varToPlot, addTitle, axisName, axisUnit, nBins, lo, hi, cutstringData, cutstringMC, sigCut, sideCut, fitres.sig, fitres.bgr, drawLine1, drawLine2);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(), "", cp.getCurPdf());
    if (flgDoHtmlReport) htrep->endTableRow();
}

doMassFitB001_fitresults cutPlotSidebandDataMCB0(TTree *datatree, TTree *mctree, Canvaspager &cp, std::string varToPlot, bool exclCut, int nBins, Double_t lo, Double_t hi, Cuts cutData, Cuts cutMC, string addTitle, string addHname, string axisName, string axisUnit)
{
    const std::string cutstringData = exclCut ? cutData.getCutExceptOne(varToPlot) : cutData.getCut();

    std::auto_ptr<TTree> subtree(datatree->CopyTree(cutstringData.c_str()));
    cp.cdNext("B0_massplot4sidebandsubtracted_"+addHname);
    doMassFitB001_fitresults fitres = doMassFitB001_DG(subtree.get(), 0, addTitle, flgNoTitle, flgPrelim, flgOfficialPlot, true);
    const double mass = fitres.mass;
    const double twoSigmas = fitres.twoSigmas;
    const double threeSigmas = fitres.threeSigmas;
    const string annotation = "2&sigma; range: " + toString(mass-twoSigmas) + "-" + toString(mass+twoSigmas)
	+ " 3&sigma; range: " + toString(mass-threeSigmas) + "-" + toString(mass+threeSigmas);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(), annotation, cp.getCurPdf());

    cutPlotSidebandDataMCB0(datatree, mctree, cp, varToPlot, exclCut, nBins, lo, hi, cutData, cutMC, addTitle, addHname, axisName, axisUnit, fitres);
    return fitres;
}

void cutsPlotsSidebandLbLbBar(TTree *tree, Canvaspager &cp, Cuts curCut, string addTitle, string addHname, doMassFitLb01_fitresults resMassFit)
{
    if (flgDoHtmlReport) htrep->addP(addTitle + ":", true);
    cout << "=======================================" << endl;
    cout << "Doing cut parameter plots for LbLbBar.... " << addTitle << endl;

    const double sigWindow(2.5), sideWindowClose(3.5), sideWindowFar(5.5);
    const double mass = resMassFit.mass;
    const double width = resMassFit.width;
    string sigCut = "mlb>" + toString(mass-sigWindow*width) + "&&mlb<" + toString(mass+sigWindow*width);
    string sideCut = "((mlb>" + toString(mass-sideWindowFar*width) + "&& mlb <=" + toString(mass-sideWindowClose*width) + ")||" +
			"(mlb>=" + toString(mass+sideWindowClose*width) + "&& mlb <" + toString(mass+sideWindowFar*width) + "))";

    cp.cdNext();
    doSidebandPlot(tree, "mjp_lb"+addHname, "mjp", "Jpsi mass from #Lambda_{b} "+addTitle, "m(J/#psi)", "GeV/c^{2}", 20, 2.9, 3.3, curCut.getCutExceptOne("mjp")+"&&rqpr>0", sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "mjp_lbbar"+addHname, "mjp", "Jpsi mass from #bar{#Lambda}_{b} "+addTitle, "m(J/#psi)", "GeV/c^{2}", 20, 2.9, 3.3, curCut.getCutExceptOne("mjp")+"&&rqpr<0", sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ml0_lb"+addHname, "ml0", "#Lambda^{0} mass from #Lambda_{b} "+addTitle, "m(#Lambda^{0})", "GeV/c^{2}", 20, 1.091, 1.141, curCut.getCutExceptOne("ml0")+"&&rqpr>0", sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ml0_lbbar"+addHname, "ml0", "#Lambda^{0} mass from #bar{#Lambda}_{b} "+addTitle, "m(#Lambda^{0})", "GeV/c^{2}", 20, 1.091, 1.141, curCut.getCutExceptOne("ml0")+"&&rqpr<0", sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");
}

// ===================================================================================
void cutsPlotsSidebandB0(TTree *tree, Canvaspager &cp, Cuts curCut, string addTitle, string addHname, doMassFitB001_fitresults resMassFit)
{
    if (flgDoHtmlReport) htrep->addP(addTitle + ":", true);
    cout << "=======================================" << endl;
    cout << "Doing cut parameter plots for B0.... " << addTitle << endl;

    const double sigWindow(2.5), sideWindowClose(3.5), sideWindowFar(5.5);
    const double mass = resMassFit.mass;
    const double width = resMassFit.width;
    string sigCut = "mB0>" + toString(mass-sigWindow*width) + "&&mB0<" + toString(mass+sigWindow*width);
    string sideCut = "((mB0>" + toString(mass-sideWindowFar*width) + "&& mB0 <=" + toString(mass-sideWindowClose*width) + ")||" +
			"(mB0>=" + toString(mass+sideWindowClose*width) + "&& mB0 <" + toString(mass+sideWindowFar*width) + "))";
    //string sideCut = "(mB0<=" + toString(mass-sideWindow*width) + "||mB0>=" + toString(mass+sigWindow*width) + ")";

    cp.cdNext();
    doSidebandPlot(tree, "ct3dB0_all"+addHname, "ct3dB0", "B^{0} lifetime 3d "+addTitle, "#tau_{3d}(B^{0})", "s", 20, 0, 9e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlotLogY(tree, "ct3dB0_alllog"+addHname, "ct3dB0", "B^{0} lifetime 3d (log) "+addTitle, "#tau_{3d}(B^{0})", "s", 30, 0, 9e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ct3dB0_red"+addHname, "ct3dB0-3*ct3dB0E", "B^{0} reduced lifetime 3d "+addTitle, "#tau_{3d}(B^{0})", "s", 20, 0, 9e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ct3dB0E"+addHname, "ct3dB0E", "B^{0} lifetime uncertainty 3d "+addTitle, "#varepsilon_{#tau_{3d}}(B^{0})", "s", 30, 0, 3e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ctxyB0_all"+addHname, "ctxyB0", "B^{0} lifetime 2d "+addTitle, "#tau_{xy}(B^{0})", "s", 20, 0, 9e-12, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "d3B0_all"+addHname, "d3B0", "B^{0} flightlength 3d"+addTitle, "d_{3d}(B^{0})", "cm", 20, 0, 1, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlotLogY(tree, "d3B0_alllog"+addHname, "d3B0", "B^{0} flightlength 3d (log) "+addTitle, "d_{3d}(B^{0})", "cm", 20, 0, 1, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "d3EB0_all"+addHname, "d3EB0", "B^{0} flightlength uncertainty 3d"+addTitle, "#varepsilon_{d_{3d}(B^{0})}", "cm", 20, 0, 0.5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "d3B0_sig"+addHname, "d3B0/d3EB0", "B^{0} flightlength significance 3d "+addTitle, "d_{3d}(B^{0}) significance", "", 20, 0, 5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "mjp"+addHname, "mjp", "Jpsi mass "+addTitle, "m(J/#psi)", "GeV/c^{2}", 20, 2.9, 3.3, curCut.getCutExceptOne("mjp"), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "mKs"+addHname, "mKs", "K_{s} mass "+addTitle, "m(K_{s})", "GeV/c^{2}", 20, 0.4476, 0.5476, curCut.getCutExceptOne("mKs"), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "maxdocaB0"+addHname, "maxdocaB0", "max doca B0 "+addTitle, "max d_{ca}(B^{0})", "cm", 25, 0, .25, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ptB0"+addHname, "ptB0", "pt B0 "+addTitle, "p_{T}(B^{0})", "GeV/c", 25, 0, 50, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "etaB0"+addHname, "TMath::Abs(etaB0)", "eta B0 "+addTitle, "|#eta(B^{0})|", "", 20, 0, 2.5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "ptgangDRKs"+addHname, "ptgangDRKs", "pointing angle Ks as dR "+addTitle, "dR_{pointing}(K_{s})", "", 50, 0, .1, curCut.getCutExceptOne("ptgangDRKs"), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "dRmumu"+addHname, "dRmumu", "dR between muons "+addTitle, "dR_(#mu#mu})", "", 50, 0, 1, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip"+addHname, "PvLip", "longit I.P. to best PV "+addTitle, "PV lip", "cm", 25, 0, .2, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip2"+addHname, "PvLip2", "longit I.P. to 2nd best PV "+addTitle, "PV lip", "cm", 25, 0, 2, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLipSig"+addHname, "PvLip/PvLipE", "longit I.P. to best PV "+addTitle, "PV lip significance", "", 25, 0, 5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip2Sig"+addHname, "PvLip2/PvLipE2", "longit I.P. to 2nd best PV "+addTitle, "PV lip significance", "", 25, 0, 5, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip2Sig_2"+addHname, "PvLip2/PvLipE2", "longit I.P. to 2nd best PV "+addTitle, "PV lip significance", "", 25, 0, 50, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");

    cp.cdNext();
    doSidebandPlot(tree, "PvLip2Sig_inv"+addHname, "PvLipE2/PvLip2", "longit I.P. to 2nd best PV "+addTitle, "PV lip inv. significance", "", 25, 0, 1, curCut.getCut(), sigCut, sideCut, resMassFit.sig, resMassFit.bgr);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"");
}

// ===================================================================================
void doAnalysis01(std::string mainTitle, bool doPublicationGrade = false)
{
    //if (filename.size() > 0) cout << "Filename option currently not used - configure files in code." << endl;
    flgPrelim = true;			    // true: prints CMS preliminary on certain plots
    flgNoTitle = true;//doPublicationGrade;	    // suppresses titles on histograms - ok for report, bad for overview
    flgOfficialPlot = false;                 // steers some other stuff
    const string outpath = "doAnalysis01plots" + toString(gSystem->Now().AsString());
    if (gSystem->MakeDirectory(outpath.c_str()) != 0)
    {
	cout << "Unable to create directory for output... Aborting." << endl;
	return;
    }
    const std::string outpdf = "doAnalysis01plots"; // stem of the output files
    std::string addCut("");			    // ad-hoc cut, string appended to the ones from cut selection
    //std::string addCut("HLTokDisplJpsi==1&&HLTmatch==1");			    // ad-hoc cut, string appended to the ones from cut selection
    //std::string addCut("HLTokDisplJpsi&&HLTmatch==1");			    // ad-hoc cut, string appended to the ones from cut selection
    //std::string addCut("HLTokBarrelJpsi==1&&HLTmatch==1");			    // ad-hoc cut, string appended to the ones from cut selection

    // Start HTML report, makes only sense if we create output with single plots per file, so doPublicationGrade is required
    flgDoHtmlReport = false;
    if (doPublicationGrade)
    {
	flgDoHtmlReport = true;
	htrep = new HtmlReport(outpath + "/index.html", "doAnalyis01 output" + (mainTitle.size() > 0 ? (" - " + mainTitle) : "" ));
	htrep->addH("Report",'1');
	if (mainTitle.size() > 0) htrep->addP(mainTitle);
	cout << "Creating HTML report as well." << endl;
    }

    // Selection of plot groups
    const bool doLb(true); // plots for Lb
    const bool doB0(true); // plots for Lb
    const bool doMC(true); // do MC plots as well
    const bool doMassPlots(true); // mass plots including fit
    const bool doMassPlotsBinned(false); // mass plots including fit
    const bool doMassPlotsHLT(false); // mass plots including fit
    const bool doMassPlotsPtYbins(false); // mass plots including fit
    const bool doLbLbBarPlots(false); // several plots regarding Lb and LbBar
    const bool doMCtruthPlots(false); // MC truth plots
    const bool doCutsPlots(false); // cut plots
    const bool doCuts2dPlots(false); // cut plots 2d
    const bool doSidebandPlots(true); // cut plots using improved sidebandsubtraction
    const bool doTriggerPlots(false); // trigger sanity plot
    const bool doEfficiencyPlotFitterPlots(false); // make efficiency plots
    const bool doJsonFiles(false); // generate json files on events selected in mass plots
    const bool doPickFiles(false); // generate pick files on events selected in mass plots
    const bool doLifetimeFit(false); // lifetime fits
    const bool doTPjpsi(false); // T&P trigger efficiencies
    const bool doMCchannelPlots(false);
    const bool doB0tTruthPlots(false); // checks sanity of lifetime reconstruction
    const bool doLbtTruthPlots(false); // checks sanity of lifetime reconstruction
    const bool doBgrChannelPlots(false); // background plots for Lb

    // Select cuts from cuts.C
    Cuts cutPreselLb, cutPreselB0; // these cuts are applied to the preselection, nothing can be wider than this
    //cutPreselLb.selectCut("acc04Lb");
    //cutPreselB0.selectCut("acc04B0");
    Cuts cutAnalLb, cutAnalLbMC, cutAnalB0, cutAnalB0MC; // these are the cuts on top the preselected ones, may also disable certain ones
    const string triggerSel = "HLT_jpsiBarrel";
    //const string triggerSel = "HLT_jpsiDispl";
    cutAnalLb.selectCut("lb12", "acc06Lb", "muSoft", "HLT_matched", triggerSel);
    cutAnalLbMC.selectCut("lb12", "acc06Lb", "muSoft", "HLT_matched", triggerSel);
    cutAnalB0.selectCut("B006", "acc06B0", "muSoft", "HLT_matched", triggerSel);
    cutAnalB0MC.selectCut("B006", "acc06B0", "muSoft", "HLT_matched", triggerSel);

    cout << "=======================================================================" << endl;
    cout << "This is rooFitCut using the following cuts as a basis:" << endl;
    cout << "Preselected cuts:" << endl;
    cout << cutPreselLb.getCut() << endl;
    cout << "Analysis cuts:" << endl;
    cout << cutAnalLb.getCut() << endl;
    cout << "=======================================================================" << endl;

    // Open file
    /*
    TFile* file = TFile::Open(filename.c_str());
    if (file ==0)
    {
        cout << "File " << filename << " not found. Exiting." << endl;
        return;
    }
    cout << "File " << filename << " succesfully opened." << endl;
    */

    string path = "../data/";
    //string path = "/scratch/frmeier/";
    Datafiles filesLb;
    //filesLb.add(Datafile(path+"run173.root","Bla",62498308.488+85771229.232,889059+1063909,0,0));
    //filesLb.add(Datafile(path+"run174.root","Bla",182316289.168,2712130,0,0));
    //filesLb.add(Datafile(path+"run175.root","Bla",221490476.336,5720870,0,0));
    //filesLb.add(Datafile(path+"run176.root","Bla",202502392.026,3088807,0,0));

    //filesLb.add(Datafile(path+"run194.root","Bla",515837510.122,15135557+14825014,0,0));
    //filesLb.add(Datafile(path+"run195.root","Bla",203016656.418,0,0,0));
    //filesLb.add(Datafile(path+"run196.root","Bla",166731765.211,0,0,0));
    //filesLb.add(Datafile(path+"run197.root","Bla",236158511.176,0,0,0));

    //filesLb.add(Datafile(path+"run265.root","Run2011A_v4",952.7e6,0,0,0));
    //filesLb.add(Datafile(path+"run266.root","Run2011A_v5",411.9e6,0,0,0));
    //filesLb.add(Datafile(path+"run267.root","Run2011A_v6",125.8e6,0,0,0));
    //filesLb.add(Datafile(path+"run268.root","Run2011A_v6",540.9e6,0,0,0));
    //filesLb.add(Datafile(path+"run269.root","Run2011B_v1",1556e6,0,0,0));

    //filesLb.add(Datafile(path+"run344.root","Run2011 incomplete!!!!",4913e6,0,0,0));

    //filesLb.add(Datafile(path+"run350.root","Run2011A_MayReReco_v1 ",236e6,0,0,0));
    //filesLb.add(Datafile(path+"run351.root","Run2011A_PromptReco_v4 ",958e6,0,0,0));
    //filesLb.add(Datafile(path+"run352.root","Run2011A_PromptReco_v5 ",412e6,0,0,0));
    //filesLb.add(Datafile(path+"run353.root","Run2011A_PromptReco_v6 ",660e6,0,0,0));
    //filesLb.add(Datafile(path+"run354.root","Run2011B_PromptReco_v1 ",2481e6,0,0,0));

    //filesLb.add(Datafile(path+"run380__run384.root","Run2011B_PromptReco_v1 ",2481e6,0,0,0));
    //filesLb.add(Datafile(path+"run434.root","Run2011A_MayReReco_v1 ",2481e6,0,0,0));
    //filesLb.add(Datafile(path+"run429.root","Run2011A_PromptReco_v4 ",2481e6,0,0,0));
    //filesLb.add(Datafile(path+"run430.root","Run2011A_PromptReco_v5 ",2481e6,0,0,0));
    //filesLb.add(Datafile(path+"run431.root","Run2011A_PromptReco_v6 ",2481e6,0,0,0));
    //filesLb.add(Datafile(path+"run435.root","Run2011B_PromptReco_v1 ",2481e6,0,0,0));

    filesLb.add(Datafile(path+"run458.root","Run2011B_PromptReco_v1 ",2481e6,0,0,0));

    Datafiles filesB0;
    //filesB0.add(Datafile(path+"run179.root","Bla",515837510.122+203016656.418,15135557+14825014,0,0));

    //filesB0.add(Datafile(path+"run198.root","Bla",515837510.122,0,0,0));
    //filesB0.add(Datafile(path+"run199.root","Bla",203016656.418,0,0,0));
    //filesB0.add(Datafile(path+"run200.root","Bla",166731765.211,0,0,0));
    //filesB0.add(Datafile(path+"run201.root","Bla",236158511.176,0,0,0));

    //filesB0.add(Datafile(path+"run259.root","Run2011A_v4",952.7e6,0,0,0));
    //filesB0.add(Datafile(path+"run260.root","Run2011A_v5",411.9e6,0,0,0));
    //filesB0.add(Datafile(path+"run261.root","Run2011A_v6",125.8e6,0,0,0));
    //filesB0.add(Datafile(path+"run262.root","Run2011A_v6",540.9e6,0,0,0));
    //filesB0.add(Datafile(path+"run263.root","Run2011B_v1",1556e6,0,0,0));

    //filesB0.add(Datafile(path+"run393__run397.root","Run2011",4913e6,0,0,0));
    //filesB0.add(Datafile(path+"run432.root","Run2011A_MayReReco",4913e6,0,0,0));
    //filesB0.add(Datafile(path+"run423.root","Run2011A_PromptReco_v4",4913e6,0,0,0));
    //filesB0.add(Datafile(path+"run424.root","Run2011A_PromptReco_v5",4913e6,0,0,0));
    //filesB0.add(Datafile(path+"run425.root","Run2011A_PromptReco_v6",4913e6,0,0,0));
    //filesB0.add(Datafile(path+"run433.root","Run2011B_PromptReco_v1",4913e6,0,0,0));
    filesB0.add(Datafile(path+"run461.root","Run2011B_PromptReco_v1",4913e6,0,0,0));

    Datafiles filesLbMC;
    //filesLbMC.add(Datafile(path+"run193.root","Bla",0,0,0,0));
    //filesLbMC.add(Datafile(path+"run385__run387.root","Bla",0,0,0,0));
    //filesLbMC.add(Datafile(path+"run407.root","Bla",0,0,0,0));
    //filesLbMC.add(Datafile(path+"run413.root","Bla",0,0,0,0));
    //filesLbMC.add(Datafile(path+"run414.root","Bla",0,0,0,0));
    //filesLbMC.add(Datafile(path+"run415.root","Bla",0,0,0,0));
    filesLbMC.add(Datafile(path+"run456.root","Bla",0,0,0,0));

    Datafiles filesB0MC;
    //filesB0MC.add(Datafile(path+"run206.root","Bla",0,0,0,0));
    //filesB0MC.add(Datafile(path+"run255.root","Bla",0,0,0,0));
    //filesB0MC.add(Datafile(path+"run398__run399.root","Bla",5e9,0,0,0));
    //filesB0MC.add(Datafile(path+"run408.root","Bla",5e9,0,0,0));
    //filesB0MC.add(Datafile(path+"run412.root","Bla",5e9,0,0,0));
    //filesB0MC.add(Datafile(path+"run417.root","Bla",5e9,0,0,0));
    //filesB0MC.add(Datafile(path+"run418.root","Bla",5e9,0,0,0));
    //filesB0MC.add(Datafile(path+"run419.root","Bla",5e9,0,0,0));
    filesB0MC.add(Datafile(path+"run460.root","Bla",5e9,0,0,0));

    // Datafiles for Lb background MC
    Datafiles filesLbBgrMC_B0;
    filesLbBgrMC_B0.add(Datafile(path+"run309.root", "B0ToPsiMuMu", 958e6,0,0,0));
    Datafiles filesLbBgrMC_Bp;
    filesLbBgrMC_Bp.add(Datafile(path+"run308.root", "BpToPsiMuMu", 883e6,0,0,0));
    Datafiles filesLbBgrMC_Bs;
    filesLbBgrMC_Bs.add(Datafile(path+"run307.root", "BsToPsiMuMu", 868e6,0,0,0));
    Datafiles filesLbBgrMC_Jp;
    filesLbBgrMC_Jp.add(Datafile(path+"run306.root", "JpToPsiMuMu", 18.6e6,0,0,0));
    Datafiles filesLbBgrMC_Xi;
    filesLbBgrMC_Xi.add(Datafile(path+"run310.root", "XiBToPsiMuMu", 5000e6,0,0,0));
    Datafiles filesLbBgrMC_Om;
    filesLbBgrMC_Om.add(Datafile(path+"run311.root", "OmegaBToPsiMuMu", 5000e6,0,0,0));

    // Datafiles for Lb background MC
    Datafiles filesB0BgrMC_Lb;
    filesB0BgrMC_Lb.add(Datafile(path+"run372.root", "LbToPsiMuMu", 4267.3,0,0,0));
    Datafiles filesB0BgrMC_Bp;
    filesB0BgrMC_Bp.add(Datafile(path+"run378.root", "BpToPsiMuMu", 882.6,0,0,0));
    Datafiles filesB0BgrMC_Bs;
    filesB0BgrMC_Bs.add(Datafile(path+"run377.root", "BsToPsiMuMu", 640.2,0,0,0));
    Datafiles filesB0BgrMC_Jp;
    filesB0BgrMC_Jp.add(Datafile(path+"run376.root", "JpToPsiMuMu", 96.5,0,0,0));

    //Datafiles filesTPjpsi;
    //filesTPjpsi.add(Datafile(path+"PromptRecov42011_00.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"PromptRecov42011_01.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"PromptRecov42011_02.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"PromptRecov42011_03.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"PromptRecov42011_04.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"May10ReReco2011_00.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"May10ReReco2011_01.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"May10ReReco2011_02.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"May10ReReco2011_03.root","T&P file",0,0,0,0));
    //filesTPjpsi.add(Datafile(path+"May10ReReco2011_04.root","T&P file",0,0,0,0));

    if (flgDoHtmlReport)
    {
	htrep->addH("Data samples used",'2');
	filesLb.mkHtmlReport(htrep,"&Lambda;<sub>b</sub> data");
	filesB0.mkHtmlReport(htrep,"B<sup>0</sup> data");
    }


    // Get tree and set up chain
    std::string treename = "events";
    TChain *chainLbData = filesLb.getTChain(treename);
    TChain *chainB0Data = filesB0.getTChain(treename);
    TChain *chainLbMC = filesLbMC.getTChain(treename);
    TChain *chainB0MC = filesB0MC.getTChain(treename);
    TChain *chainLbBgrMC_B0 = filesLbBgrMC_B0.getTChain(treename);
    TChain *chainLbBgrMC_Bp = filesLbBgrMC_Bp.getTChain(treename);
    TChain *chainLbBgrMC_Bs = filesLbBgrMC_Bs.getTChain(treename);
    TChain *chainLbBgrMC_Jp = filesLbBgrMC_Jp.getTChain(treename);
    TChain *chainLbBgrMC_Xi = filesLbBgrMC_Xi.getTChain(treename);
    TChain *chainLbBgrMC_Om = filesLbBgrMC_Om.getTChain(treename);
    TChain *chainB0BgrMC_Lb = filesB0BgrMC_Lb.getTChain(treename);
    TChain *chainB0BgrMC_Bp = filesB0BgrMC_Bp.getTChain(treename);
    TChain *chainB0BgrMC_Bs = filesB0BgrMC_Bs.getTChain(treename);
    TChain *chainB0BgrMC_Jp = filesB0BgrMC_Jp.getTChain(treename);

    // initialise the TCanvas
    setTDRStyle();
    const int maxCanvasCols = doPublicationGrade ? 1 : 3;
    const int maxCanvasRows = doPublicationGrade ? 1 : 6;
    const int nPlots = maxCanvasCols*maxCanvasRows;
    int nCanvasCols = maxCanvasCols;
    int nCanvasRows = maxCanvasRows;
    if (nPlots < maxCanvasCols*maxCanvasRows)
    {
        if(nPlots == 1)
            nCanvasCols = 1;
        else
        {
            nCanvasCols = maxCanvasCols;
            nCanvasRows = TMath::Ceil(nPlots / nCanvasCols);
        }
    }
    int nCanvasCd = nCanvasCols*nCanvasRows;
    const int nPads = nCanvasCd;
    canvas = new TCanvas("canvas","doAnalysis01_canvas",nCanvasCols*500,nCanvasRows*300);
    canvas->Divide(nCanvasCols,nCanvasRows);

    // check if we can make a preselection
    TTree* treeLbData= 0;
    if (doLb)
    {
	TFile* tmpfile = new TFile("tmp.root","RECREATE");
	if (tmpfile == 0)
	{
	    cout << "temporary file creation failed - exiting" << endl;
	    return;
	}
	cout << "tmpfile: " << tmpfile << endl;
	if(cutPreselLb.getCut().size() > 0) // we have a preselection cut, so use it
	{
	    cout << "Preselecting events..." << endl;
	    treeLbData= chainLbData->CopyTree(cutPreselLb.getCut().c_str());
	    cout << "now " << treeLbData->GetEntries() << " of " << chainLbData->GetEntries()<< endl;
	    delete chainLbData;
	}
	else
	{
	    cout << "No preselection cut given, working on full treeLbData" << endl;
	    treeLbData = chainLbData;
	}
    }

    // check if we can make a preselection
    TTree* treeB0Data= 0;
    if (doB0)
    {
	TFile* tmpfileB0 = new TFile("tmpB0.root","RECREATE");
	if (tmpfileB0 == 0)
	{
	    cout << "temporary file creation failed - exiting" << endl;
	    return;
	}
	cout << "tmpfileB0: " << tmpfileB0 << endl;
	if(cutPreselB0.getCut().size() > 0) // we have a preselection cut, so use it
	{
	    cout << "Preselecting events..." << endl;
	    treeB0Data= chainB0Data->CopyTree(cutPreselB0.getCut().c_str());
	    cout << "now " << treeB0Data->GetEntries() << " of " << chainB0Data->GetEntries()<< endl;
	    delete chainB0Data;
	    tmpfileB0->Write();
	}
	else
	{
	    cout << "No preselection cut given, working on full treeB0Data" << endl;
	    treeB0Data = chainB0Data;
	}
    }

    // check if we can make a preselection
    TTree* treeLbMC= 0;
    if (doLb && doMC)
    {
	TFile* tmpfileLbMC = new TFile("tmpLbMC.root","RECREATE");
	if (tmpfileLbMC == 0)
	{
	    cout << "temporary file creation failed - exiting" << endl;
	    return;
	}
	cout << "tmpfileLbMC: " << tmpfileLbMC << endl;
	if(cutPreselLb.getCut().size() > 0) // we have a preselection cut, so use it
	{
	    cout << "Preselecting events..." << endl;
	    treeLbMC= chainLbMC->CopyTree(cutPreselLb.getCut().c_str());
	    cout << "now " << treeLbMC->GetEntries() << " of " << chainLbMC->GetEntries()<< endl;
	    delete chainLbMC;
	}
	else
	{
	    cout << "No preselection cut given, working on full treeLbMC" << endl;
	    treeLbMC = chainLbMC;
	}
    }

    // check if we can make a preselection
    TTree* treeB0MC= 0;
    if (doB0 && doMC)
    {
	TFile* tmpfileB0MC = new TFile("tmpB0MC.root","RECREATE");
	if (tmpfileB0MC == 0)
	{
	    cout << "temporary file creation failed - exiting" << endl;
	    return;
	}
	cout << "tmpfileB0MC: " << tmpfileB0MC << endl;
	if(cutPreselB0.getCut().size() > 0) // we have a preselection cut, so use it
	{
	    cout << "Preselecting events..." << endl;
	    treeB0MC= chainB0MC->CopyTree(cutPreselB0.getCut().c_str());
	    cout << "now " << treeB0MC->GetEntries() << " of " << chainB0MC->GetEntries()<< endl;
	    delete chainB0MC;
	}
	else
	{
	    cout << "No preselection cut given, working on full treeB0MC" << endl;
	    treeB0MC = chainB0MC;
	}
    }

    // initialise page
    Canvaspager canvaspager(canvas, outpath, outpdf, nPads);

    // ===================================================================================
    // json files
    if (doJsonFiles)
    {
	if (doLb) writeJsonFile(treeLbData, cutAnalLb.getCut() + "&&mlb>5.30&&mlb<5.90", outpath + "/lbData.json");
	if (doB0) writeJsonFile(treeB0Data, cutAnalB0.getCut() + "&&mB0>5.00&&mB0<5.60", outpath + "/B0Data.json");
    }

    // ===================================================================================
    // json files
    if (doPickFiles)
    {
	if (doLb) writePickFile(treeLbData, cutAnalLb.getCut() + "&&mlb>5.30&&mlb<6.10", outpath + "/lbPicklist.txt");
	if (doB0) writePickFile(treeB0Data, cutAnalB0.getCut() + "&&mB0>5.00&&mB0<5.60", outpath + "/B0Picklist.txt");
    }

    // ===================================================================================
    // mass plot
    doMassFitLb01_fitresults resMassFitLb, resMassFitLb_HLTJpsi, resMassFitLb_HLTJpsiDispl, resMassFitLb_HLTJpsiBarrel;
    doMassFitB001_fitresults resMassFitB0, resMassFitB0_HLTJpsi, resMassFitB0_HLTJpsiDispl, resMassFitB0_HLTJpsiBarrel;
    doMassFitLb01_fitresults resMassFitLbMC;
    doMassFitB001_fitresults resMassFitB0MC;
    string curCut = (addCut=="") ? cutAnalLb.getCut() : (cutAnalLb.getCut() + "&&" + addCut);
    if (doMassPlots && doLb)
    {
	cout << "=======================================" << endl;
	cout << "Doing mass plots for Lb...." << endl;

	canvaspager.cdNext();

	// Apply cut selection
	cout << "Using cut: " << curCut << endl;
	std::auto_ptr<TTree> subtree(treeLbData->CopyTree(curCut.c_str()));
	cout << "subtree: " << subtree->GetEntries() << " entries of " << treeLbData->GetEntries()<< endl;
	const double curEntries = subtree->GetEntries();
	cout << "curEntries: " << curEntries << endl;

	if (flgDoHtmlReport) htrep->addH("Mass plots &Lambda;<sub>b</sub>",'2',true);
	resMassFitLb = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "full dataset", flgNoTitle, flgPrelim, flgOfficialPlot, true);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - full dataset (no trigger selection)");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&HLTmatch==1").c_str()));
	doMassFitLb01_fitresults resMassFitLb_both= doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT matched, #Lambda_{b}+#bar{#Lambda}_{b}", flgNoTitle, flgPrelim, flgOfficialPlot, true);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&rqha1>0&&HLTmatch==1").c_str()));
	doMassFitLb01_fitresults resMassFitLb_Lb= doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT matched, #Lambda_{b} only", flgNoTitle, flgPrelim, flgOfficialPlot, true);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched, &Lambda;<sub>b</sub> only");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&rqha1<0&&HLTmatch==1").c_str()));
	doMassFitLb01_fitresults resMassFitLb_LbBar = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT matched, #bar{#Lambda}_{b} only", flgNoTitle, flgPrelim, flgOfficialPlot, true);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched, anti-&Lambda;<sub>b</sub> only");

	if (flgDoHtmlReport) htrep->addH("Mass plots J/&psi; in &Lambda;<sub>b</sub> date",'3', true);
	canvaspager.cdNext();
	doMassFitJp01_fitresults resMassFitLbjp = doMassFitJp01(subtree.get(), filesLb.getLumiPbRounded(), "data full dataset, #Lambda_{b} #rightarrow J/#psi (#mu#mu) K_{s} (#pi#pi)", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"J/&psi; mass in &Lambda;<sub>b</sub> data full dataset (signal only)");
    }

    if (doMassPlots && doB0)
    {
	cout << "=======================================" << endl;
	cout << "Doing mass plots for B0...." << endl;

	canvaspager.cdNext();

	// Apply cut selection
	cout << "Using cut: " << cutAnalB0.getCut() << endl;
	std::auto_ptr<TTree> subtree(treeB0Data->CopyTree(cutAnalB0.getCut().c_str()));
	cout << "subtree: " << subtree->GetEntries() << " entries of " << treeB0Data->GetEntries()<< endl;
	const double curEntries = subtree->GetEntries();
	cout << "curEntries: " << curEntries << endl;

	if (flgDoHtmlReport) htrep->addH("Mass plots B<sup>0</sup>",'2', true);
	resMassFitB0 = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "full dataset", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - full dataset (no trigger selection)");

	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&HLTmatch==1").c_str()));
	doMassFitB001_fitresults resMassFitB0_hltm = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT matched, B^{0}", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLT matched");

	if (flgDoHtmlReport) htrep->addH("Mass plots J/&psi; in B<sup>0</sup> data",'3', true);
	canvaspager.cdNext();
	doMassFitJp01_fitresults resMassFitB0jp = doMassFitJp01(subtree.get(), filesB0.getLumiPbRounded(), "data full dataset, B^{0} #rightarrow J/#psi (#mu#mu) K_{s} (#pi#pi)", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"J/&psi; mass in B<sup>0</sup> data full dataset (signal only)");

	if (doMassPlotsBinned)
	{
	//------------------------
	if (flgDoHtmlReport) htrep->addP("Some binned plots in p<sub>T</sub>:", true);
	double binsPt[] = {0, 11.9198, 14.1159, 16.7824, 20.9707, 91.5254};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "ptB0>=" + toString(binsPt[i]) + "&&ptB0<" + toString(binsPt[i+1]);
	    const string binTitle = "p_{T}(B0)#in[" + toString(binsPt[i]) + "," + toString(binsPt[i+1]) + ")";
	    subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0.getLumiPbRounded(), "data, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"data signal, pt bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in &eta;:", true);
	double binsEta[] = {-2.4, -1.07, -.37, .30, 1.03, 2.4};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "etaB0>=" + toString(binsEta[i]) + "&&etaB0<" + toString(binsEta[i+1]);
	    const string binTitle = "#eta(B0)#in[" + toString(binsEta[i]) + "," + toString(binsEta[i+1]) + ")";
	    subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0.getLumiPbRounded(), "data, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"data signal, eta bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in |&eta;|:", true);
	double binsAbsEta[] = {0, .336, .660, 1.048, 1.572, 2.332};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "TMath::Abs(etaB0)>=" + toString(binsAbsEta[i]) + "&&TMath::Abs(etaB0)<" + toString(binsAbsEta[i+1]);
	    const string binTitle = "|#eta(B0)|#in[" + toString(binsAbsEta[i]) + "," + toString(binsAbsEta[i+1]) + ")";
	    subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0.getLumiPbRounded(), "data, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"data signal, eta bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in &phi;:", true);
	double binsPhi[] = {-TMath::Pi(), -.6*TMath::Pi(), -.2*TMath::Pi(), .2*TMath::Pi(), .6*TMath::Pi(), TMath::Pi()};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "phiB0>=" + toString(binsPhi[i]) + "&&phiB0<" + toString(binsPhi[i+1]);
	    const string binTitle = "#phi(B0)#in[" + toString(binsPhi[i]) + "," + toString(binsPhi[i+1]) + ")";
	    subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0.getLumiPbRounded(), "data, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"data signal, phi bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in t:", true);
	double binsT[] = {-20e-12,.1906e-12,.6676e-12,1.361e-12,2.437e-12,20e-12};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "ct3dB0>=" + toString(binsT[i]) + "&&ct3dB0<" + toString(binsT[i+1]);
	    const string binTitle = "t_{3d}(B0)#in[" + toString(binsT[i]) + "," + toString(binsT[i+1]) + ")";
	    subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0.getLumiPbRounded(), "data, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"data signal, phi bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in tE:", true);
	double binsTE[] = {-5e-12,.09707e-12,0.1477e-12,0.2371e-12,0.7715e-12,8e-12};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "ct3dB0E>=" + toString(binsTE[i]) + "&&ct3dB0E<" + toString(binsTE[i+1]);
	    const string binTitle = "tE_{3d}(B0)#in[" + toString(binsTE[i]) + "," + toString(binsTE[i+1]) + ")";
	    subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0.getLumiPbRounded(), "data, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"data signal, phi bin");
	}
	}
    }

    if (doMassPlots && doLb && doMC)
    {
	cout << "=======================================" << endl;
	cout << "Doing mass plots for Lb MC...." << endl;

	canvaspager.cdNext();

	// Apply cut selection
	cout << "Using cut: " << cutAnalLbMC.getCut() << endl;
	std::auto_ptr<TTree> subtree(treeLbMC->CopyTree(cutAnalLbMC.getCut().c_str()));
	cout << "subtree: " << subtree->GetEntries() << " entries of " << treeLbMC->GetEntries()<< endl;
	const double curEntries = subtree->GetEntries();
	cout << "curEntries: " << curEntries << endl;

	if (flgDoHtmlReport) htrep->addH("Mass plots &Lambda;<sub>b</sub> MC",'2', true);
	if (flgDoHtmlReport) htrep->addP("<strong>NOTE:</strong> mass used for MC is 5.624 GeV/c<sup>2</sup> and not the current best value 6.620 GeV/c<sup>2</sup>.");
	resMassFitLbMC = doMassFitLb01(subtree.get(), filesLbMC.getLumiPbRounded(), "full dataset", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - MC full dataset (no trigger selection)");

	canvaspager.cdNext();
	//subtree.reset(treeLbMC->CopyTree((curCut+"&&HLTmatch==1").c_str()));
	subtree.reset(treeLbMC->CopyTree((curCut).c_str()));
	doMassFitLb01_fitresults resMassFitLb_both= doMassFitLb01(subtree.get(), filesLbMC.getLumiPbRounded(), "MC, #Lambda_{b}+#bar{#Lambda}_{b}", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched");

	canvaspager.cdNext();
	subtree.reset(treeLbMC->CopyTree((curCut+"&&rqha1>0&&HLTmatch==1").c_str()));
	doMassFitLb01_fitresults resMassFitLb_Lb= doMassFitLb01(subtree.get(), filesLbMC.getLumiPbRounded(), "MC, #Lambda_{b} only", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched, &Lambda;<sub>b</sub> only");

	canvaspager.cdNext();
	subtree.reset(treeLbMC->CopyTree((curCut+"&&rqha1<0&&HLTmatch==1").c_str()));
	doMassFitLb01_fitresults resMassFitLb_LbBar = doMassFitLb01(subtree.get(), filesLbMC.getLumiPbRounded(), "MC, #bar{#Lambda}_{b} only", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched, anti-&Lambda;<sub>b</sub> only");

	if (flgDoHtmlReport) htrep->addH("Mass plots J/&psi; in &Lambda;<sub>b</sub> MC",'3', true);
	canvaspager.cdNext();
	subtree.reset(treeLbMC->CopyTree((cutAnalLbMC.getCut()+"&&isSig==1").c_str()));
	doMassFitJp01_fitresults resMassFitLbjp_isSig = doMassFitJp01(subtree.get(), filesLbMC.getLumiPbRounded(), "MC full dataset, #Lambda_{b} #rightarrow J/#psi (#mu#mu) #Lambda^{0} (#pi#pi)", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"J/&psi; mass in &Lambda;<sub>b</sub> data full dataset (signal only)");
    }

    if (doMassPlots && doB0 && doMC)
    {
	cout << "=======================================" << endl;
	cout << "Doing mass plots for B0 MC...." << endl;

	canvaspager.cdNext();

	// Apply cut selection
	cout << "Using cut: " << cutAnalB0MC.getCut() << endl;
	std::auto_ptr<TTree> subtree(treeB0MC->CopyTree(cutAnalB0MC.getCut().c_str()));
	cout << "subtree: " << subtree->GetEntries() << " entries of " << treeB0MC->GetEntries()<< endl;
	const double curEntries = subtree->GetEntries();
	cout << "curEntries: " << curEntries << endl;

	if (flgDoHtmlReport) htrep->addH("Mass plots B<sup>0</sup> MC",'2', true);
	resMassFitB0MC = doMassFitB001(subtree.get(), filesB0MC.getLumiPbRounded(), "B^{0} MC full dataset", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sub>0</sub> mass - MC full dataset (no trigger selection)");

	canvaspager.cdNext();
	subtree.reset(treeB0MC->CopyTree((cutAnalB0MC.getCut()+"&&isSig==1").c_str()));
	doMassFitB001_fitresults resMassFitB0_isSig= doMassFitB001(subtree.get(), filesB0MC.getLumiPbRounded(), "MC full dataset, B^{0} #rightarrow J/#psi (#mu#mu) K_{s} (#pi#pi)", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - MC full dataset (signal only)");

	if (flgDoHtmlReport) htrep->addH("Mass plots J/&psi; in B<sup>0</sup> MC",'3', true);
	canvaspager.cdNext();
	subtree.reset(treeB0MC->CopyTree((cutAnalB0MC.getCut()+"&&isSig==1").c_str()));
	doMassFitJp01_fitresults resMassFitB0jp_isSig = doMassFitJp01(subtree.get(), filesB0MC.getLumiPbRounded(), "MC full dataset, B^{0} #rightarrow J/#psi (#mu#mu) K_{s} (#pi#pi)", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"J/&psi; mass in B<sup>0</sup> MC full dataset (signal only)");

	if (doMassPlotsBinned)
	{
	if (flgDoHtmlReport) htrep->addP("Some binned plots in p<sub>T</sub>:", true);
	double binsPt[] = {0, 11.9198, 14.1159, 16.7824, 20.9707, 91.5254};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "ptB0>=" + toString(binsPt[i]) + "&&ptB0<" + toString(binsPt[i+1]);
	    const string binTitle = "p_{T}(B0)#in[" + toString(binsPt[i]) + "," + toString(binsPt[i+1]) + ")";
	    subtree.reset(treeB0MC->CopyTree((cutAnalB0MC.getCut()+"&&isSig&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0MC.getLumiPbRounded(), "MC, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"MC signal, pt bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in &eta;:", true);
	double binsEta[] = {-2.4, -1.07, -.37, .30, 1.03, 2.4};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "etaB0>=" + toString(binsEta[i]) + "&&etaB0<" + toString(binsEta[i+1]);
	    const string binTitle = "#eta(B0)#in[" + toString(binsEta[i]) + "," + toString(binsEta[i+1]) + ")";
	    subtree.reset(treeB0MC->CopyTree((cutAnalB0MC.getCut()+"&&isSig&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0MC.getLumiPbRounded(), "MC, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"MC signal, eta bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in |&eta;|:", true);
	double binsAbsEta[] = {0, .336, .660, 1.048, 1.572, 2.332};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "TMath::Abs(etaB0)>=" + toString(binsAbsEta[i]) + "&&TMath::Abs(etaB0)<" + toString(binsAbsEta[i+1]);
	    const string binTitle = "|#eta(B0)|#in[" + toString(binsAbsEta[i]) + "," + toString(binsAbsEta[i+1]) + ")";
	    subtree.reset(treeB0MC->CopyTree((cutAnalB0MC.getCut()+"&&isSig&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0MC.getLumiPbRounded(), "MC, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"MC signal, eta bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in &phi;:", true);
	double binsPhi[] = {-TMath::Pi(), -.6*TMath::Pi(), -.2*TMath::Pi(), .2*TMath::Pi(), .6*TMath::Pi(), TMath::Pi()};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "phiB0>=" + toString(binsPhi[i]) + "&&phiB0<" + toString(binsPhi[i+1]);
	    const string binTitle = "#phi(B0)#in[" + toString(binsPhi[i]) + "," + toString(binsPhi[i+1]) + ")";
	    subtree.reset(treeB0MC->CopyTree((cutAnalB0MC.getCut()+"&&isSig&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0MC.getLumiPbRounded(), "MC, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"MC signal, phi bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in t:", true);
	double binsT[] = {-20e-12,.1906e-12,.6676e-12,1.361e-12,2.437e-12,20e-12};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "ct3dB0>=" + toString(binsT[i]) + "&&ct3dB0<" + toString(binsT[i+1]);
	    const string binTitle = "t_{3d}(B0)#in[" + toString(binsT[i]) + "," + toString(binsT[i+1]) + ")";
	    subtree.reset(treeB0MC->CopyTree((cutAnalB0MC.getCut()+"&&isSig&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0MC.getLumiPbRounded(), "MC, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"MC signal, phi bin");
	}

	if (flgDoHtmlReport) htrep->addP("Some binned plots in tE:", true);
	double binsTE[] = {-5e-12,.09707e-12,0.1477e-12,0.2371e-12,0.7715e-12,8e-12};
	for (int i = 0; i!=5; i++)
	{
	    canvaspager.cdNext();
	    const string binCut = "ct3dB0E>=" + toString(binsTE[i]) + "&&ct3dB0E<" + toString(binsTE[i+1]);
	    const string binTitle = "tE_{3d}(B0)#in[" + toString(binsTE[i]) + "," + toString(binsTE[i+1]) + ")";
	    subtree.reset(treeB0MC->CopyTree((cutAnalB0MC.getCut()+"&&isSig&&"+binCut).c_str()));
	    doMassFitB001_CB(subtree.get(), filesB0MC.getLumiPbRounded(), "MC, signal only, "+binTitle, flgNoTitle, flgPrelim, flgOfficialPlot);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"MC signal, phi bin");
	}
	}
    }

    if (doLb && doMassPlotsHLT){
	if (flgDoHtmlReport) htrep->addH("Mass plots &Lambda;<sub>b</sub> with trigger selection",'2');
	std::auto_ptr<TTree> subtree;
	// do plots with different trigger selections
	Cuts cutHLTJpsi, cutHLTJpsiDispl, cutHLTJpsiBarrel;
	cutHLTJpsi.selectCut("HLT_matched","HLT_jpsi");
	cutHLTJpsiDispl.selectCut("HLT_matched","HLT_jpsiDispl");
	cutHLTJpsiBarrel.selectCut("HLT_matched","HLT_jpsiBarrel");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&"+cutHLTJpsi.getCut()).c_str()));
	resMassFitLb_HLTJpsi = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT jpsi, matched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched, Jpsi triggers");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&"+cutHLTJpsiDispl.getCut()).c_str()));
	resMassFitLb_HLTJpsiDispl = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT jpsi displaced, matched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched, displaced Jpsi trigger");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&"+cutHLTJpsiBarrel.getCut()).c_str()));
	resMassFitLb_HLTJpsiBarrel = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT jpsi barrel, matched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLT matched, barrel Jpsi trigger");

	doMassFitLb01_fitresults resMassDummy;
	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&"+"HLTSingleMu==1").c_str()));
	resMassDummy = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT HLTSingleMu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLTSingleMu");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&"+"HLTSingleIsoMu==1").c_str()));
	resMassDummy = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT HLTSingleIsoMu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLTSingleIsoMu");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&"+"HLTSingleL1Mu==1").c_str()));
	resMassDummy = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT HLTSingleL1Mu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLTSingleL1Mu");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&"+"HLTSingleL2Mu==1").c_str()));
	resMassDummy = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT HLTSingleL2Mu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLTSingleL2Mu");

	canvaspager.cdNext();
	subtree.reset(treeLbData->CopyTree((curCut+"&&"+"HLTSingleHLTMu==1").c_str()));
	resMassDummy = doMassFitLb01_DG(subtree.get(), filesLb.getLumiPbRounded(), "HLT HLTSingleHLTMu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - HLTSingleHLTMu");
    }

    if (doB0 && doMassPlotsHLT){
	if (flgDoHtmlReport) htrep->addH("Mass plots B<sup>0</sup> with trigger selection",'2');
	std::auto_ptr<TTree> subtree;
	// do plots with different trigger selections
	Cuts cutHLTJpsi, cutHLTJpsiDispl, cutHLTJpsiBarrel;
	cutHLTJpsi.selectCut("HLT_matched","HLT_jpsi");
	cutHLTJpsiDispl.selectCut("HLT_matched","HLT_jpsiDispl");
	cutHLTJpsiBarrel.selectCut("HLT_matched","HLT_jpsiBarrel");

	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0+cutHLTJpsi).getCut().c_str()));
	resMassFitB0_HLTJpsi = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT jpsi, matched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLT matched, Jpsi triggers");

	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0+cutHLTJpsiDispl).getCut().c_str()));
	resMassFitB0_HLTJpsiDispl = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT jpsi displaced, matched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLT matched, displaced Jpsi trigger");

	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0+cutHLTJpsiBarrel).getCut().c_str()));
	resMassFitB0_HLTJpsiBarrel = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT jpsi barrel, matched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLT matched, barrel Jpsi trigger");

	doMassFitB001_fitresults resMassDummy;
	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+"HLTSingleMu==1").c_str()));
	resMassDummy = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT HLTSingleMu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLTSingleMu");

	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+"HLTSingleIsoMu==1").c_str()));
	resMassDummy = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT HLTSingleIsoMu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLTSingleIsoMu");

	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+"HLTSingleL1Mu==1").c_str()));
	resMassDummy = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT HLTSingleL1Mu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLTSingleL1Mu");

	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+"HLTSingleL2Mu==1").c_str()));
	resMassDummy = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT HLTSingleL2Mu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLTSingleL2Mu");

	canvaspager.cdNext();
	subtree.reset(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&"+"HLTSingleHLTMu==1").c_str()));
	resMassDummy = doMassFitB001_DG(subtree.get(), filesB0.getLumiPbRounded(), "HLT HLTSingleHLTMu, unmatched", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - HLTSingleHLTMu");
    }

    if (doLb && doMassPlotsPtYbins)
    {
	if (flgDoHtmlReport) htrep->addH("&Lambda;<sub>b</sub> mass plots in kinematic bins",'2');
	// do plots with different trigger selections
	Cuts cutEmpty, cutHLTJpsi, cutHLTJpsiDispl, cutHLTJpsiBarrel;
	cutHLTJpsi.selectCut("HLT_matched","HLT_jpsi");
	cutHLTJpsiDispl.selectCut("HLT_matched","HLT_jpsiDispl");
	cutHLTJpsiBarrel.selectCut("HLT_matched","HLT_jpsiBarrel");

	if (flgDoHtmlReport) htrep->addP("no trigger selection");
	massPlotsPtYbind(treeLbData, canvaspager, cutAnalLb+cutEmpty, "", filesLb.getLumiPbRounded());
	if (flgDoHtmlReport) htrep->addP("HLT matched, Jpsi triggers");
	massPlotsPtYbind(treeLbData, canvaspager, cutAnalLb+cutHLTJpsi, "HLT jpsi matched", filesLb.getLumiPbRounded());
	if (flgDoHtmlReport) htrep->addP("HLT matched, Jpsi displaced trigger");
	massPlotsPtYbind(treeLbData, canvaspager, cutAnalLb+cutHLTJpsiDispl, "HLT jpsi displ matched", filesLb.getLumiPbRounded());
	if (flgDoHtmlReport) htrep->addP("HLT matched, Jpsi barrel trigger");
	massPlotsPtYbind(treeLbData, canvaspager, cutAnalLb+cutHLTJpsiBarrel, "HLT jpsi barrel matched", filesLb.getLumiPbRounded());
    }

    // ===================================================================================
    // plots of distributions of certain cut parameters
    if (doCutsPlots && doMassPlots && doLb) // requires mass plots result
    {
	if (flgDoHtmlReport) htrep->addH("&Lambda;<sub>b</sub> cut variables plots",'2', true);
	// do plots with different trigger selections
	Cuts cutEmpty, cutHLTJpsi, cutHLTJpsiDispl, cutHLTJpsiBarrel;
	cutHLTJpsi.selectCut("HLT_matched","HLT_jpsi");
	cutHLTJpsiDispl.selectCut("HLT_matched","HLT_jpsiDispl");
	cutHLTJpsiBarrel.selectCut("HLT_matched","HLT_jpsiBarrel");

	//cutsPlotsSidebandLb(treeLbData, canvaspager, cutAnalLb, "", "", resMassFitLb);
	if (doMassPlotsHLT)
	{
	    //cutsPlotsSidebandLb(treeLbData, canvaspager, cutAnalLb+cutHLTJpsi, "HLT Jpsi", "jpsi", resMassFitLb_HLTJpsi);
	    cutsPlotsSidebandLb(treeLbData, canvaspager, cutAnalLb+cutHLTJpsiDispl, "HLT Jpsi displ", "disp", resMassFitLb_HLTJpsiDispl);
	    cutsPlotsSidebandLb(treeLbData, canvaspager, cutAnalLb+cutHLTJpsiBarrel, "HLT Jpsi barrel", "barr", resMassFitLb_HLTJpsiBarrel);
	}
    }

    // ===================================================================================
    // Sidebandsubracted plots
    if (doSidebandPlots && doLb && doMC)
    {
	if (flgDoHtmlReport) htrep->addH("&Lambda;<sub>b</sub> cut variables plots",'2', true);
	// do plots with different trigger selections
	Cuts cutIsSig;
	cutIsSig.selectCut("isSig","isMCmatch");
	Cuts cutData = cutAnalLb;
	//cutData.cs.addCut(new cutConst("HLTDMu6p5BarJp==0")); cutData.parvec.push_back(0);
	Cuts cutMC = cutAnalLb+cutIsSig;
	if (triggerSel=="HLT_jpsiBarrel") { cutMC.cs.addCut(new cutConst("HLTokBarrelJpsiMC==1")); cutMC.parvec.push_back(0); }
	if (triggerSel=="HLT_jpsiDispl") { cutMC.cs.addCut(new cutConst("HLTokDisplJpsiMC==1")); cutMC.parvec.push_back(0); }
	// special case for d3rs requires another cut
	Cuts cutAnalLb_nod3rs = cutAnalLb;
	cutAnalLb_nod3rs.removeOneCut("d3rs");
	cutAnalLb_nod3rs.removeOneCut("d3rs/d3Ers");
	Cuts cutAnalLb_nod3rs_data = cutAnalLb_nod3rs;
	Cuts cutAnalLb_nod3rs_mc   = cutAnalLb_nod3rs + cutIsSig;

	const string addTitle(triggerSel == "HLT_jpsiBarrel" ? "Barrel" : "Displ");
	const string addHname(triggerSel == "HLT_jpsiBarrel" ? "Barrel" : "Displ");
	if (flgDoHtmlReport) htrep->addP(addTitle + ":", true);
	cout << "=======================================" << endl;
	cout << "Doing sidebandsubtracted plots for Lb.... " << addTitle << endl;
	doMassFitLb01_fitresults fitresCache; // use only for those cases where the cuts do not change, i.e. where they were not excluded
	fitresCache = cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "ptrs", false, 40, 0, 20, cutData, cutMC, "p_{T}(#Lambda^{0}) "+addTitle, "ptrs"+addHname, "p_{T}(#Lambda^{0} cand.)", "GeV/c");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "ptjp", false, 40, 0, 40, cutData, cutMC, "p_{T}(J/#psi) "+addTitle, "ptjp"+addHname, "p_{T}(J/#psi cand.)", "GeV/c", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "ptbc", false, 40, 0, 40, cutData, cutMC, "p_{T}(#Lambda_{b}) "+addTitle, "ptbc"+addHname, "p_{T}(#Lambda_{b} cand.)", "GeV/c", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "mjp", false, 60, 2.8, 3.4, cutData, cutMC, "m(J/#psi) "+addTitle, "mjp"+addHname, "m(#mu#mu)", "GeV/c^{2}", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "mrs", true, 30, 1.10, 1.13, cutData, cutMC, "m(#Lambda^{0}) "+addTitle, "mrs"+addHname, "m(p#pi)", "GeV/c^{2}");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "Kshypo", true, 40, 0.45, 0.55, cutData, cutMC, "m(K_{s}) hypothesis"+addTitle, "Kshypo"+addHname, "m(#pi#pi)", "GeV/c^{2}", fitresCache); // fit unstable, so use values from cache
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "probjp", false, 25, 0, 1, cutData, cutMC, "prob(J/#psi) "+addTitle, "probjp"+addHname, "prob(J/#psi)", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "probrs", true, 25, 0, 1, cutData, cutMC, "prob(#Lambda^{0}) "+addTitle, "probrs"+addHname, "prob(#Lambda^{0})", "");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "probmu1", false, 25, 0, 1, cutData, cutMC, "prob(#mu_{1}) "+addTitle, "probmu1"+addHname, "prob(#mu_{1})", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "probmu2", false, 25, 0, 1, cutData, cutMC, "prob(#mu_{2}) "+addTitle, "probmu2"+addHname, "prob(#mu_{2})", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "probha1", false, 25, 0, 1, cutData, cutMC, "prob(p) "+addTitle, "probha1"+addHname, "prob(p)", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "probha2", false, 25, 0, 1, cutData, cutMC, "prob(#pi}) "+addTitle, "probha2"+addHname, "prob(#pi)", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "rptmu1", true, 30, 0, 30, cutData, cutMC, "p_{T}(#mu_{1}) "+addTitle, "rptmu1"+addHname, "p_{T}(#mu_{1})", "GeV/c");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "rptmu2", true, 30, 0, 30, cutData, cutMC, "p_{T}(#mu_{2}) "+addTitle, "rptmu2"+addHname, "p_{T}(#mu_{2})", "GeV/c");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "rptha1", true, 20, 0, 20, cutData, cutMC, "p_{T}(p) "+addTitle, "rptha1"+addHname, "p_{T}(p)", "GeV/c");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "rptha2", true, 25, 0, 10, cutData, cutMC, "p_{T}(#pi) "+addTitle, "rptha2"+addHname, "p_{T}(#pi)", "GeV/c");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "rptha1-rptha2", false, 32, -1, 15, cutData, cutMC, "p_{T}(p)-p_{T}(#pi) "+addTitle, "rptprha2"+addHname, "p_{T}(p)-p_{T}(#pi)", "GeV/c", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "d3rs", true, 60, 0, 60, cutData, cutMC, "d_{3}(#Lambda^{0}) "+addTitle, "d3rs"+addHname, "d_{3}(#Lambda^{0})", "cm");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "d3rs", false, 60, 0, 60, cutAnalLb_nod3rs_data, cutAnalLb_nod3rs_mc, "d_{3}(#Lambda^{0}) no sign cut "+addTitle, "d3rsnosig"+addHname, "d_{3}(#Lambda^{0})", "cm", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "d3rs/d3Ers", true, 25, 0, 100, cutData, cutMC, "significance d_{3}(#Lambda^{0}) "+addTitle, "d3rs"+addHname, "significance d_{3}(#Lambda^{0})", "#sigma");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "d3rs/d3Ers", false, 25, 0, 100, cutAnalLb_nod3rs_data, cutAnalLb_nod3rs_mc, "significance d_{3}(#Lambda^{0}) no d_{3}(#Lambda^{0}) cut"+addTitle, "d3rs"+addHname, "significance d_{3}(#Lambda^{0})", "#sigma", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "ptgangDRrs", false, 60, 0, 0.03, cutData, cutMC, "#alpha(#Lambda^{0}) "+addTitle, "ptgangDRrs"+addHname, "#alpha(#Lambda^{0})", "#Delta R", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "ptgangDRrs", false, 60, 0, 0.1, cutData, cutMC, "#alpha(#Lambda^{0}) "+addTitle, "ptgangDRrs_2"+addHname, "#alpha(#Lambda^{0})", "#Delta R", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "alphars", true, 60, 0, 0.03, cutData, cutMC, "#alpha(#Lambda^{0}) "+addTitle, "alphars"+addHname, "#alpha(#Lambda^{0})", "rad");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "alphars", true, 25, 0, 0.005, cutData, cutMC, "#alpha(#Lambda^{0}) "+addTitle, "alphars"+addHname, "#alpha(#Lambda^{0})", "rad");
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "alphabc", false, 60, 0, 0.3, cutData, cutMC, "#alpha(#Lambda_{b}) "+addTitle, "alphabc"+addHname, "#alpha(#Lambda_{b})", "rad", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "alphabc", false, 25, 0, 0.05, cutData, cutMC, "#alpha(#Lambda_{b}) "+addTitle, "alphabc"+addHname, "#alpha(#Lambda_{b})", "rad", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "retamu1", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(#mu_{1}) "+addTitle, "retamu1"+addHname, "#eta(#mu_{1})", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "retamu2", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(#mu_{2}) "+addTitle, "retamu2"+addHname, "#eta(#mu_{2})", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "retaha1", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(p) "+addTitle, "retaha1"+addHname, "#eta(p)", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "retaha2", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(#pi) "+addTitle, "retaha2"+addHname, "#eta(#pi)", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "etars", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(#Lambda^{0}) "+addTitle, "etars"+addHname, "#eta(#Lambda^{0})", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "etajp", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(J/#psi) "+addTitle, "etajp"+addHname, "#eta(J/#psi)", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "ybc", false, 25, -2.5, 2.5, cutData, cutMC, "y(#Lambda_{b}) "+addTitle, "ybc"+addHname, "x(#Lambda_{b})", "", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "d3bc", false, 50, 0, 1, cutData, cutMC, "d_{3}(#Lambda_{b}) "+addTitle, "d3bc"+addHname, "d_{3}(#Lambda_{b})", "cm", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "pbc", false, 30, 0, 60, cutData, cutMC, "p(#Lambda_{b}) "+addTitle, "pbc"+addHname, "p(#Lambda_{b})", "GeV/c", fitresCache);
	cutPlotSidebandDataMCLb(treeLbData, treeLbMC, canvaspager, "ct3dbc", false, 50, -5e-12, 20e-12, cutData, cutMC, "t_{3}(#Lambda_{b}) "+addTitle, "ct3dbc"+addHname, "t_{3}(#Lambda_{b})", "s", fitresCache);
    }

    if (doSidebandPlots && doB0 && doMC)
    {
	if (flgDoHtmlReport) htrep->addH("B<sup>0</sup> cut variables plots",'2', true);
	// do plots with different trigger selections
	Cuts cutIsSig;
	cutIsSig.selectCut("isSig","isMCmatch");
	Cuts cutData = cutAnalB0;
	Cuts cutMC = cutAnalB0+cutIsSig;
	if (triggerSel=="HLT_jpsiBarrel") { cutMC.cs.addCut(new cutConst("HLTokBarrelJpsiMC==1")); cutMC.parvec.push_back(0); }
	if (triggerSel=="HLT_jpsiDispl") { cutMC.cs.addCut(new cutConst("HLTokDisplJpsiMC==1")); cutMC.parvec.push_back(0); }
	// special case for d3l0 requires another cut
	Cuts cutAnalB0_nod3rs = cutAnalB0;
	cutAnalB0_nod3rs.removeOneCut("d3rs");
	cutAnalB0_nod3rs.removeOneCut("d3rs/d3Ers");
	Cuts cutAnalB0_nod3rs_data = cutAnalB0_nod3rs;
	Cuts cutAnalB0_nod3rs_mc   = cutAnalB0_nod3rs + cutIsSig;

	const string addTitle(triggerSel == "HLT_jpsiBarrel" ? "Barrel" : "Displ");
	const string addHname(triggerSel == "HLT_jpsiBarrel" ? "Barrel" : "Displ");
	if (flgDoHtmlReport) htrep->addP(addTitle + ":", true);
	cout << "=======================================" << endl;
	cout << "Doing sidebandsubtracted plots for B0.... " << addTitle << endl;
	doMassFitB001_fitresults fitresCache; // use only for those cases where the cuts do not change, i.e. where they were not excluded
	fitresCache = cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "ptrs", false, 40, 0, 20, cutData, cutMC, "p_{T}(K_{s}) "+addTitle, "ptrs"+addHname, "p_{T}(K_{s} cand.)", "GeV/c");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "ptjp", false, 40, 0, 40, cutData, cutMC, "p_{T}(J/#psi) "+addTitle, "ptjp"+addHname, "p_{T}(J/#psi cand.)", "GeV/c", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "ptbc", false, 40, 0, 40, cutData, cutMC, "p_{T}(B^{0}) "+addTitle, "ptbc"+addHname, "p_{T}(B^{0} cand.)", "GeV/c", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "mjp", false, 60, 2.8, 3.4, cutData, cutMC, "m(J/#psi) "+addTitle, "mjp"+addHname, "m(#mu#mu)", "GeV/c^{2}", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "mrs", true, 30, 0.467, 0.527, cutData, cutMC, "m(K_{s}) "+addTitle, "mrs"+addHname, "m(#pi#pi)", "GeV/c^{2}");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "L0hypo", true, 40, 1.10, 1.13, cutData, cutMC, "m(#Lambda^{0}) hypothesis"+addTitle, "L0hypo"+addHname, "m(p#pi)", "GeV/c^{2}");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "probjp", false, 25, 0, 1, cutData, cutMC, "prob(J/#psi) "+addTitle, "probjp"+addHname, "prob(J/#psi)", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "probrs", true, 25, 0, 1, cutData, cutMC, "prob(K_{s}) "+addTitle, "probrs"+addHname, "prob(K_{s})", "");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "probmu1", false, 25, 0, 1, cutData, cutMC, "prob(#mu_{1}) "+addTitle, "probmu1"+addHname, "prob(#mu_{1})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "probmu2", false, 25, 0, 1, cutData, cutMC, "prob(#mu_{2}) "+addTitle, "probmu2"+addHname, "prob(#mu_{2})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "probha1", false, 25, 0, 1, cutData, cutMC, "prob(#pi_{1}) "+addTitle, "probha1"+addHname, "prob(#pi_{1})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "probha2", false, 25, 0, 1, cutData, cutMC, "prob(#pi_{2})) "+addTitle, "probha2"+addHname, "prob(#pi_{2})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "rptmu1", true, 30, 0, 30, cutData, cutMC, "p_{T}(#mu_{1}) "+addTitle, "rptmu1"+addHname, "p_{T}(#mu_{1})", "GeV/c");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "rptmu2", true, 30, 0, 30, cutData, cutMC, "p_{T}(#mu_{2}) "+addTitle, "rptmu2"+addHname, "p_{T}(#mu_{2})", "GeV/c");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "rptha1", true, 25, 0, 10, cutData, cutMC, "p_{T}(#pi_{1}) "+addTitle, "rptha1"+addHname, "p_{T}(#pi_{1})", "GeV/c");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "rptha2", true, 25, 0, 10, cutData, cutMC, "p_{T}(#pi_{2}) "+addTitle, "rptha2"+addHname, "p_{T}(#pi_{2})", "GeV/c");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "d3rs", true, 60, 0, 60, cutData, cutMC, "d_{3}(K_{s}) "+addTitle, "d3rs"+addHname, "d_{3}(K_{s})", "cm");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "d3rs", false, 60, 0, 60, cutAnalB0_nod3rs_data, cutAnalB0_nod3rs_mc, "d_{3}(K_{s}) no sign cut "+addTitle, "d3rsnosig"+addHname, "d_{3}(K_{s})", "cm", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "d3rs/d3Ers", true, 25, 0, 100, cutData, cutMC, "significance d_{3}(K_{s}) "+addTitle, "d3rs"+addHname, "significance d_{3}(K_{s})", "#sigma");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "d3rs/d3Ers", false, 25, 0, 100, cutAnalB0_nod3rs_data, cutAnalB0_nod3rs_mc, "significance d_{3}(K_{s}) no d_{3}(K_{s}) cut"+addTitle, "d3rs"+addHname, "significance d_{3}(K_{s})", "#sigma", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "ptgangDRrs", false, 60, 0, 0.03, cutData, cutMC, "#alpha(K_{s}) "+addTitle, "ptgangDRrs"+addHname, "#alpha(K_{s})", "#Delta R", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "ptgangDRrs", false, 60, 0, 0.1, cutData, cutMC, "#alpha(K_{s}) "+addTitle, "ptgangDRrs_2"+addHname, "#alpha(K_{s})", "#Delta R", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "alphars", true, 60, 0, 0.03, cutData, cutMC, "#alpha(K_{s}) "+addTitle, "alphars"+addHname, "#alpha(K_{s})", "rad");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "alphars", true, 25, 0, 0.005, cutData, cutMC, "#alpha(K_{s}) "+addTitle, "alphars"+addHname, "#alpha(K_{s})", "rad");
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "alphabc", false, 60, 0, 0.3, cutData, cutMC, "#alpha(B^{0}) "+addTitle, "alphabc"+addHname, "#alpha(B^{0})", "rad", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "alphabc", false, 25, 0, 0.05, cutData, cutMC, "#alpha(B^{0}) "+addTitle, "alphabc"+addHname, "#alpha(B^{0})", "rad", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "retamu1", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(#mu_{1}) "+addTitle, "retamu1"+addHname, "#eta(#mu_{1})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "retamu2", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(#mu_{2}) "+addTitle, "retamu2"+addHname, "#eta(#mu_{2})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "retaha1", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(#pi_{1}) "+addTitle, "retaha1"+addHname, "#eta(#pi_{1})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "retaha2", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(#pi_{2}) "+addTitle, "retaha2"+addHname, "#eta(#pi_{2})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "etars", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(K_{s}) "+addTitle, "etars"+addHname, "#eta(K_{s})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "etajp", false, 25, -2.5, 2.5, cutData, cutMC, "#eta(J/#psi) "+addTitle, "etajp"+addHname, "#eta(J/#psi)", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "ybc", false, 25, -2.5, 2.5, cutData, cutMC, "y(B^{0}) "+addTitle, "ybc"+addHname, "x(B^{0})", "", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "PvLip", false, 50, 0, .05, cutData, cutMC, "d_{lipPV}(B^{0}) "+addTitle, "PvLip"+addHname, "d_{lipPV}(B^{0})", "cm", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "PvLip/PvLipE", false, 50, 0, 20, cutData, cutMC, "significance d_{lipPV}(B^{0}) "+addTitle, "PvLipS"+addHname, "sign d_{lipPV}(B^{0})", "#sigma", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "PvLip2", false, 50, 0, 5, cutData, cutMC, "d_{lipPV2nd}(B^{0}) "+addTitle, "PvLip2"+addHname, "d_{lipPV}(B^{0})", "cm", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "PvLip2", false, 50, 0, .05, cutData, cutMC, "d_{lipPV2nd}(B^{0}) "+addTitle, "PvLip2z"+addHname, "d_{lipPV}(B^{0})", "cm", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "PvLip2/PvLipE2", false, 50, 0, 100, cutData, cutMC, "significance d_{lipPV2nd}(B^{0}) "+addTitle, "PvLip2S"+addHname, "d_{lipPV}(B^{0})", "#sigma", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "PvLip2/PvLipE2", false, 50, 0, 20, cutData, cutMC, "significance d_{lipPV2nd}(B^{0}) "+addTitle, "PvLip2Sz"+addHname, "d_{lipPV}(B^{0})", "#sigma", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "d3bc", false, 50, 0, 1, cutData, cutMC, "d_{3}(B^{0}) "+addTitle, "d3bc"+addHname, "d_{3}(B^{0})", "cm", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "pbc", false, 30, 0, 60, cutData, cutMC, "p(B^{0}) "+addTitle, "pbc"+addHname, "p(B^{0})", "GeV/c", fitresCache);
	cutPlotSidebandDataMCB0(treeB0Data, treeB0MC, canvaspager, "ct3dbc", false, 50, -5e-12, 20e-12, cutData, cutMC, "t_{3}(B^{0}) "+addTitle, "ct3dbc"+addHname, "t_{3}(B^{0})", "s", fitresCache);
    }

    // ===================================================================================
    // plots about Lb and LbBar
    if (doLbLbBarPlots && doMassPlots && doLb) // requires mass plots result
    {
	if (flgDoHtmlReport) htrep->addH("&Lambda;<sub>b</sub> <div style=\"text-decoration: overline;\">&Lambda;</div><sub>b</sub> cut variables plots",'2', true);
	// do plots with different trigger selections
	Cuts cutEmpty, cutHLTJpsi, cutHLTJpsiDispl, cutHLTJpsiBarrel;
	cutHLTJpsi.selectCut("HLT_matched","HLT_jpsi");
	cutHLTJpsiDispl.selectCut("HLT_matched","HLT_jpsiDispl");
	cutHLTJpsiBarrel.selectCut("HLT_matched","HLT_jpsiBarrel");

	cutsPlotsSidebandLbLbBar(treeLbData, canvaspager, cutAnalLb, "", "", resMassFitLb);
	if (doMassPlotsHLT)
	{
	    cutsPlotsSidebandLbLbBar(treeLbData, canvaspager, cutAnalLb+cutHLTJpsi, "HLT Jpsi", "jpsi", resMassFitLb_HLTJpsi);
	    cutsPlotsSidebandLbLbBar(treeLbData, canvaspager, cutAnalLb+cutHLTJpsiDispl, "HLT Jpsi displ", "disp", resMassFitLb_HLTJpsiDispl);
	    cutsPlotsSidebandLbLbBar(treeLbData, canvaspager, cutAnalLb+cutHLTJpsiBarrel, "HLT Jpsi barrel", "barr", resMassFitLb_HLTJpsiBarrel);
	}
    }

    // ===================================================================================
    // plots of distributions of certain cut parameters
    if (doCutsPlots && doMassPlots && doLb && doMC) // requires mass plots result
    {
	if (flgDoHtmlReport) htrep->addH("&Lambda;<sub>b</sub> cut variables plots - Monte-Carlo",'2');
	// do plots with different trigger selections
	Cuts cutEmpty, cutHLTJpsi, cutHLTJpsiDispl, cutHLTJpsiBarrel;
	cutHLTJpsi.selectCut("HLT_matched","HLT_jpsi");
	cutHLTJpsiDispl.selectCut("HLT_matched","HLT_jpsiDispl");
	cutHLTJpsiBarrel.selectCut("HLT_matched","HLT_jpsiBarrel");

	cutsPlotsSidebandLb(treeLbMC, canvaspager, cutAnalLbMC, "MC", "MC", resMassFitLbMC);
	if (doMassPlotsHLT)
	{
	    //cutsPlotsSidebandLb(treeLbData, canvaspager, cutAnalLb+cutHLTJpsi, "HLT Jpsi", "jpsi", resMassFitLb_HLTJpsi);
	    //cutsPlotsSidebandLb(treeLbData, canvaspager, cutAnalLb+cutHLTJpsiDispl, "HLT Jpsi displ", "disp", resMassFitLb_HLTJpsiDispl);
	    //cutsPlotsSidebandLb(treeLbData, canvaspager, cutAnalLb+cutHLTJpsiBarrel, "HLT Jpsi barrel", "barr", resMassFitLb_HLTJpsiBarrel);
	}
    }

    // ===================================================================================
    // plots of distributions of certain cut parameters
    if (doCutsPlots && doMassPlots && doB0) // requires mass plots result
    {
	if (flgDoHtmlReport) htrep->addH("B<sup>0</sup> cut variables plots",'2');
	// do plots with different trigger selections
	Cuts cutEmpty, cutHLTJpsi, cutHLTJpsiDispl, cutHLTJpsiBarrel;
	cutHLTJpsi.selectCut("HLT_matched","HLT_jpsi");
	cutHLTJpsiDispl.selectCut("HLT_matched","HLT_jpsiDispl");
	cutHLTJpsiBarrel.selectCut("HLT_matched","HLT_jpsiBarrel");

	cutsPlotsSidebandB0(treeB0Data, canvaspager, cutAnalB0, "B0", "B0", resMassFitB0);
	if (doMassPlotsHLT)
	{
	    cutsPlotsSidebandB0(treeB0Data, canvaspager, cutAnalB0+cutHLTJpsi, "B0 HLT Jpsi", "B0jpsi", resMassFitB0_HLTJpsi);
	    cutsPlotsSidebandB0(treeB0Data, canvaspager, cutAnalB0+cutHLTJpsiDispl, "B0 HLT Jpsi displ", "B0disp", resMassFitB0_HLTJpsiDispl);
	    cutsPlotsSidebandB0(treeB0Data, canvaspager, cutAnalB0+cutHLTJpsiBarrel, "B0 HLT Jpsi barrel", "B0barr", resMassFitB0_HLTJpsiBarrel);
	}
    }

    // ===================================================================================
    // 1d pots of MC truth comparisons
    if (doMCtruthPlots && doLb && doMC)
    {
	if (flgDoHtmlReport) htrep->addH("&Lambda;<sub>b</sub> - Monte-Carlo truth comparison plots",'2',true);

	const double ctWindow(10e-12), ctWindowZoom(1e-12);
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCt3dLb", "ct3dlb-ctlbtruth", "", 40, -ctWindow, ctWindow, "Lifetime of #Lambda_{b} 3d, all candidates found", "#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCt3dLb2", "ct3dlb-ctlbtruth", "isSig==1&&HLTokJpsi==1", 40, -ctWindow, ctWindow, "Lifetime of #Lambda_{b} 3d, require sig and HLT", "#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCt3dLb2zoom", "ct3dlb-ctlbtruth", "isSig==1&&HLTokJpsi==1", 40, -ctWindowZoom, ctWindowZoom, "Lifetime of #Lambda_{b} 3d, require sig and HLT, zoom", "#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCt3dLb2cut", "ct3dlb-ctlbtruth", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 40, -ctWindow, ctWindow, "Lifetime of #Lambda_{b} 3d, require sig and HLT, cuts", "#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCt3dLb2cutzoom", "ct3dlb-ctlbtruth", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 40, -ctWindowZoom, ctWindowZoom, "Lifetime of #Lambda_{b} 3d, require sig and HLT, cuts, zoom", "#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport)
	{
	    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("MCtruthCt3dLb2cutzoom");
	    htrep->addTableImage(canvaspager.getCurPng(),"Mean: "+toString(h1->GetMean())+" RMS: "+toString(h1->GetRMS()));
	}

	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCt3dLb2cutPull", "(ct3dlb-ctlbtruth)/ct3dlbE", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 40, -10, 10, "Lifetime pull of #Lambda_{b} 3d, require sig and HLT, cuts, zoom", "(#tau_{reco,3d}(#Lambda_{b})-#tau_{truth}(#Lambda_{b}))/#varepsilon_{3d}", "#sigma");
	if (flgDoHtmlReport)
	{
	    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("MCtruthCt3dLb2cutPull");
	    htrep->addTableImage(canvaspager.getCurPng(),"Mean: "+toString(h1->GetMean())+" RMS: "+toString(h1->GetRMS()));
	}
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthd3LbcutPull", "(d3lb-d3dlbtruth)/d3Elb", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 40, -10, 10, "Flighlength 3d pull of #Lambda_{b}, require sig and HLT, cuts, zoom", "(d_{reco,3d}(#Lambda_{b})-d_{truth}(#Lambda_{b}))/#varepsilon_{d_{3d}}", "#sigma");
	if (flgDoHtmlReport)
	{
	    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("MCtruthd3LbcutPull");
	    htrep->addTableImage(canvaspager.getCurPng(),"Mean: "+toString(h1->GetMean())+" RMS: "+toString(h1->GetRMS()));
	}

	if (flgDoHtmlReport) htrep->addP("2d");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCtxyLb", "ctxylb-ctlbtruth", "", 40, -ctWindow, ctWindow, "Lifetime of #Lambda_{b} xy, all candidates found", "#tau_{reco,xy}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCtxyLb2", "ctxylb-ctlbtruth", "isSig==1&&HLTokJpsi==1", 40, -ctWindow, ctWindow, "Lifetime of #Lambda_{b} xy, require sig and HLT", "#tau_{reco,xy}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCtxyLb2zoom", "ctxylb-ctlbtruth", "isSig==1&&HLTokJpsi==1", 40, -ctWindowZoom, ctWindowZoom, "Lifetime of #Lambda_{b} xy, require sig and HLT, zoom", "#tau_{reco,xy}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCtxyLb2cut", "ctxylb-ctlbtruth", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 40, -ctWindow, ctWindow, "Lifetime of #Lambda_{b} xy, require sig and HLT, cuts", "#tau_{reco,xy}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCtxyLb2cutzoom", "ctxylb-ctlbtruth", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 40, -ctWindowZoom, ctWindowZoom, "Lifetime of #Lambda_{b} xy, require sig and HLT, cuts, zoom", "#tau_{reco,xy}(#Lambda_{b})-#tau_{truth}(#Lambda_{b})", "s");
	if (flgDoHtmlReport)
	{
	    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("MCtruthCtxyLb2cutzoom");
	    htrep->addTableImage(canvaspager.getCurPng(),"Mean: "+toString(h1->GetMean())+" RMS: "+toString(h1->GetRMS()));
	}

	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCtxyLb2cutPull", "(ctxylb-ctlbtruth)/ctxylbE", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 40, -10, 10, "Lifetime pull of #Lambda_{b} xy, require sig and HLT, cuts, zoom", "(#tau_{reco,xy}(#Lambda_{b})-#tau_{truth}(#Lambda_{b}))/#varepsilon_{xy}", "#sigma");
	if (flgDoHtmlReport)
	{
	    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("MCtruthCtxyLb2cutPull");
	    htrep->addTableImage(canvaspager.getCurPng(),"Mean: "+toString(h1->GetMean())+" RMS: "+toString(h1->GetRMS()));
	}
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthdxyLbcutPull", "(dxylb-dxylbtruth)/dxyElb", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 40, -10, 10, "Flighlength xy pull of #Lambda_{b}, require sig and HLT, cuts, zoom", "(d_{reco,xy}(#Lambda_{b})-d_{truth}(#Lambda_{b}))/#varepsilon_{d_{xy}}", "#sigma");
	if (flgDoHtmlReport)
	{
	    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("MCtruthdxyLbcutPull");
	    htrep->addTableImage(canvaspager.getCurPng(),"Mean: "+toString(h1->GetMean())+" RMS: "+toString(h1->GetRMS()));
	}

	if (flgDoHtmlReport) htrep->addP("additional:");
	const double ctWindow3dXy(0.3e-12);
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCtxy3dLbzoom", "ct3dlb-ctxylb", "isSig==1&&HLTokJpsi==1", 50, -ctWindow3dXy, ctWindow3dXy, "Lifetime of #Lambda_{b} 3d-xy, require sig and HLT, zoom", "#tau_{reco,3d}(#Lambda_{b})-#tau_{reco,xy}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do1dPlot(treeLbMC, "MCtruthCtxy3dLb2zoom", "ct3dlb-ctxylb", cutAnalLbMC.getCut()+"&&isSig==1&&HLTokJpsi==1", 50, -ctWindow3dXy, ctWindow3dXy, "Lifetime of #Lambda_{b} 3d-xy, require sig and HLT, cuts, zoom", "#tau_{reco,3d}(#Lambda_{b})-#tau_{reco,xy}(#Lambda_{b})", "s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
    }


    // ===================================================================================
    // 2d plots of distributions of certain cut parameters
    if (doCuts2dPlots && doLb)
    {
	// Armenteros-Podolanski plots
	if (flgDoHtmlReport) htrep->addH("Armenteros-Podolanski plots for &Lambda;<sub>b</sub>",'2');
	canvaspager.cdNext();
	do2dPlot(treeLbData, "armpod", "armQt:armAl", "", 40, -1, 1, 20, 0, 1, "Armenteros-Podolanski plot for #Lambda^{0}", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbData, "armpodcut", "armQt:armAl", cutAnalLb, 40, -1, 1, 20, 0, 1, "Armenteros-Podolanski plot #Lambda^{0}, all cuts", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbData, "armpodcutPrf", "armQt:armAl", cutAnalLb, 40, -1, 1, 0, 1, "Armenteros-Podolanski profile plot #Lambda^{0}, all cuts", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dPlot(treeLbData, "armpodUnbinned", "armQt:armAl", cutAnalLb.getCut(), -1, 1, 0, 1, "Armenteros-Podolanski plot", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	// Plots to study PV sanity
	if (flgDoHtmlReport) htrep->addH("Lifetime error vs. lifetime - &Lambda;<sub>b</sub>",'2');
	canvaspager.cdNext();
	do2dPlot(treeLbData, "ct3dlbEct3dlb", "ct3dlbE:ct3dlb", cutAnalLb.getCut(), 20, 0, 9e-12, 20, 0, 9e-12, "Lifetime error vs. lifetime (log) ", "#tau_{3d}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbData, "ct3dlbEct3dlblog", "ct3dlbE:log10(ct3dlb)", cutAnalLb.getCut(), 20, -14, -11, 20, 0, 2e-12, "Lifetime error vs. lifetime (log)", "log #tau_{3d}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbData, "ct3dlbEct3dlblogprof", "ct3dlbE:log10(ct3dlb)", cutAnalLb.getCut(), 20, -13, -10, "Lifetime error vs. lifetime 3d (log profile)", "log #tau_{3d}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dProfilePlot(treeLbData, "ct3dlbSiglogprof", "ct3dlb/ct3dlbE:log10(ct3dlb)", cutAnalLb.getCut(), 20, -13, -10, "Lifetime significance vs. lifetime 3d (log profile)", "log #tau_{3d}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbData, "ctxylbEctxylblogprof", "ctxylbE:log10(ctxylb)", cutAnalLb.getCut(), 20, -13, -10, "Lifetime error vs. lifetime xy (log profile)", "log #tau_{xy}(#Lambda_{b})","s","#varepsilon_{#tau_{xy}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dProfilePlot(treeLbData, "ctxylbSiglogprof", "ctxylb/ctxylbE:log10(ctxylb)", cutAnalLb.getCut(), 20, -13, -10, "Lifetime significance vs. lifetime xy (log profile)", "log #tau_{xy}(#Lambda_{b})","s","#varepsilon_{#tau_{xy}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbData, "ctLb3dxy", "ct3dlb:ctxylb", cutAnalLb.getCut(), 20, 0, 10e-12, "Lifetime 3d vs. xy (profile)", "#tau_{xy}(#Lambda_{b})","s","#tau_{3d}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbData, "ctLb3dExyE", "ct3dlbE:ctxylbE", cutAnalLb.getCut(), 10, 0, 10e-12, "Lifetime error 3d vs. xy (profile)", "#varepsilon_{#tau_{xy}}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	// Plots to study PV-Distributions
	if (flgDoHtmlReport) htrep->addH("Primary vertex distances for &Lambda;<sub>b</sub>",'2');
	canvaspager.cdNext();
	do2dPlot(treeLbData, "PvLip2d", "PvLip2:PvLip", "", 60, 0, 0.02, 60, 0, 1, "PV distribution", "longit. I.P. to best PV","cm","longit. I.P. to 2nd best IP","cm");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbData, "PvLip2dzoom", "PvLip2:PvLip", "", 30, 0, 0.005, 30, 0, 0.25, "PV distribution, zoom", "longit. I.P. to best PV","cm","longit. I.P. to 2nd best IP","cm");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbData, "PvLip2dcut", "PvLip2:PvLip", cutAnalLb, 20, 0, 0.02, 20, 0, 1, "PV distribution, all cuts", "longit. I.P. to best PV","cm","longit. I.P. to 2nd best IP","cm");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbData, "PvLip2dcutzoom", "PvLip2:PvLip", cutAnalLb, 10, 0, 0.005, 10, 0, 0.25, "PV distribution, all cuts, zoom", "longit. I.P. to best PV","cm","longit. I.P. to 2nd best IP","cm");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
    }

    if (doCuts2dPlots && doLb && doMC)
    {
	// Armenteros-Podolanski plots
	if (flgDoHtmlReport) htrep->addH("Armenteros-Podolanski plots for &Lambda;<sub>b</sub> Monte-Carlo",'2');
	canvaspager.cdNext();
	do2dPlot(treeLbMC, "armpodMC", "armQt:armAl", "", 40, -1, 1, 20, 0, 1, "Armenteros-Podolanski plot for #Lambda^{0}", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbMC, "armpodcutMC", "armQt:armAl", cutAnalLbMC, 40, -1, 1, 20, 0, 1, "Armenteros-Podolanski plot #Lambda^{0}, all cuts", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbMC, "armpodcutPrfMC", "armQt:armAl", cutAnalLbMC, 40, -1, 1, 0, 1, "Armenteros-Podolanski profile plot #Lambda^{0}, all cuts", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dPlot(treeLbMC, "armpodUnbinnedMC", "armQt:armAl", cutAnalLbMC.getCut(), -1, 1, 0, 1, "Armenteros-Podolanski plot", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	// Plots to study PV sanity
	if (flgDoHtmlReport) htrep->addH("Lifetime error vs. lifetime - &Lambda;<sub>b</sub> Monte-Carlo",'2');
	canvaspager.cdNext();
	do2dPlot(treeLbMC, "ct3dlbEct3dlbMC", "ct3dlbE:ct3dlb", cutAnalLbMC.getCut(), 20, 0, 9e-12, 20, 0, 9e-12, "Lifetime error vs. lifetime (log) ", "#tau_{3d}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbMC, "ct3dlbEct3dlblogMC", "ct3dlbE:log10(ct3dlb)", cutAnalLbMC.getCut(), 20, -14, -11, 20, 0, 2e-12, "Lifetime error vs. lifetime (log)", "log #tau_{3d}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbMC, "ct3dlbEct3dlblogprofMC", "ct3dlbE:log10(ct3dlb)", cutAnalLbMC.getCut(), 20, -13, -10, "Lifetime error vs. lifetime (log profile)", "log #tau_{3d}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dProfilePlot(treeLbMC, "ct3dlbSiglogprofMC", "ct3dlb/ct3dlbE:log10(ct3dlb)", cutAnalLbMC.getCut(), 20, -13, -10, "Lifetime significance vs. lifetime (log profile)", "log #tau_{3d}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbMC, "ctLbMC3dxy", "ct3dlb:ctxylb", cutAnalLb.getCut(), 20, 0, 10e-12, "Lifetime 3d vs. xy (profile)", "#tau_{xy}(#Lambda_{b})","s","#tau_{3d}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeLbMC, "ctLbMC3dExyE", "ct3dlbE:ctxylbE", cutAnalLb.getCut(), 10, 0, 10e-12, "Lifetime error 3d vs. xy (profile)", "#varepsilon_{#tau_{xy}}(#Lambda_{b})","s","#varepsilon_{#tau_{3d}}(#Lambda_{b})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dPlot(treeLbMC, "ctLbMC3d", "(ct3dlb-ctlbtruth)/ct3dlbE:ct3dlb", cutAnalLb.getCut(), 10, 0, 10e-12, 10, -0.5, 0.5, "Lifetime error pull 3d vs. lifetime", "#tau_{3d}(#Lambda_{b})","s","pull(#tau_{3d}(#Lambda_{b}))","#sigma");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dPlot(treeLbMC, "ctLbMCxy", "(ctxylb-ctlbtruth)/ctxylbE:ctxylb", cutAnalLb.getCut(), 10, 0, 10e-12, 10, -5, 5, "Lifetime error pull xy vs. lifetime", "#tau_{xy}(#Lambda_{b})","s","pull(#tau_{xy}(#Lambda_{b}))","#sigma");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	// Plots to study PV-Distributions
	if (flgDoHtmlReport) htrep->addH("Primary vertex distances for &Lambda;<sub>b</sub> Monte-Carlo",'2');
	canvaspager.cdNext();
	do2dPlot(treeLbMC, "PvLip2dMC", "PvLip2:PvLip", "", 60, 0, 0.02, 60, 0, 1, "PV distribution", "longit. I.P. to best PV","cm","longit. I.P. to 2nd best IP","cm");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbMC, "PvLip2dzoomMC", "PvLip2:PvLip", "", 30, 0, 0.005, 30, 0, 0.25, "PV distribution, zoom", "longit. I.P. to best PV","cm","longit. I.P. to 2nd best IP","cm");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbMC, "PvLip2dcutMC", "PvLip2:PvLip", cutAnalLbMC, 20, 0, 0.02, 20, 0, 1, "PV distribution, all cuts", "longit. I.P. to best PV","cm","longit. I.P. to 2nd best IP","cm");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeLbMC, "PvLip2dcutzoomMC", "PvLip2:PvLip", cutAnalLbMC, 10, 0, 0.005, 10, 0, 0.25, "PV distribution, all cuts, zoom", "longit. I.P. to best PV","cm","longit. I.P. to 2nd best IP","cm");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
    }
    if (doCuts2dPlots && doB0)
    {
	// Armenteros-Podolanski plots
	if (flgDoHtmlReport) htrep->addH("Armenteros-Podolanski plots for B<sup>0</sup>",'2');
	canvaspager.cdNext();
	do2dPlot(treeB0Data, "armpodB0", "armQt:armAl", "", 40, -1, 1, 20, 0, 1, "Armenteros-Podolanski plot for K_{s}", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeB0Data, "armpodB0cut", "armQt:armAl", cutAnalB0, 40, -1, 1, 20, 0, 1, "Armenteros-Podolanski plot K_{s}, all cuts", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeB0Data, "armpodB0cutPrf", "armQt:armAl", cutAnalB0, 40, -1, 1, 0, 1, "Armenteros-Podolanski profile plot K_{s}, all cuts", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dPlot(treeB0Data, "armpodB0Unbinned", "armQt:armAl", cutAnalB0.getCut(), -1, 1, 0, 1, "Armenteros-Podolanski plot", "(|p_{+}|-|p_{-}|) / |p_{-}|","","|p_{+}#times p_{-}| / |p_{+}+p_{-}|","");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	// Plots to study PV sanity
	if (flgDoHtmlReport) htrep->addH("Lifetime error vs. lifetime - B0<sup>0</sup>",'2');
	canvaspager.cdNext();
	do2dPlot(treeB0Data, "ct3dB0Ect3dB0", "ct3dB0E:ct3dB0", cutAnalB0.getCut(), 20, 0, 9e-12, 20, 0, 9e-12, "Lifetime error vs. lifetime (log) ", "#tau_{3d}(B^{0})","s","#varepsilon_{#tau_{3d}}(B^{0})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dPlot(treeB0Data, "ct3dB0Ect3dB0log", "ct3dB0E:log10(ct3dB0)", cutAnalB0.getCut(), 20, -14, -11, 20, 0, 2e-12, "Lifetime error vs. lifetime (log)", "log #tau_{3d}(B^{0})","s","#varepsilon_{#tau_{3d}}(B^{0})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeB0Data, "ct3dB0Ect3dB0logprof", "ct3dB0E:log10(ct3dB0)", cutAnalB0.getCut(), 20, -13, -10, "Lifetime error vs. lifetime (log profile)", "log #tau_{3d}(B^{0})","s","#varepsilon_{#tau_{3d}}(B^{0})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	canvaspager.cdNext();
	do2dProfilePlot(treeB0Data, "ct3dB0Siglogprof", "ct3dB0/ct3dB0E:log10(ct3dB0)", cutAnalB0.getCut(), 20, -13, -10, "Lifetime significance vs. lifetime (log profile)", "log #tau_{3d}(B^{0})","s","significance tau_{3d}(B^{0})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeB0Data, "ctB03dxy", "ct3dB0:ctxyB0", cutAnalB0.getCut(), 20, 0, 10e-12, "Lifetime 3d vs. xy (profile)", "#tau_{xy}(B^{0})","s","#tau_{3d}(B^{0})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

	canvaspager.cdNext();
	do2dProfilePlot(treeB0Data, "ctB03dExyE", "ct3dB0E:ctxyB0E", cutAnalB0.getCut(), 10, 0, 6e-12, "Lifetime error 3d vs. xy (profile)", "#varepsilon_{#tau_{xy}}(B^{0})","s","#varepsilon_{#tau_{3d}}(B^{0})","s");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");

    }

    // ===================================================================================
    // Trigger sanity plot
    if (doTriggerPlots && doLb)
    {
	if (flgDoHtmlReport) htrep->addH("HLT trigger sanity plots",'2');
	cout << "=======================================" << endl;
	cout << "Doing HLT sanity plots for Lb...." << endl;

	TriggerPlot trgpl(treeLbData);
	for (int i = 0; i!=trgpl.getNhistos(); i++)
	{
	    canvaspager.cdNext();
	    trgpl.drawHisto(i);
	    if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"");
	}
    }

    // ===================================================================================
    // Efficiency plots for Lb trigger efficiencies
    if (doLb && doEfficiencyPlotFitterPlots)
    {
	if (flgDoHtmlReport) htrep->addH("Trigger efficiencies plots for &Lambda;<sub>b</sub>",'2');
	//TTree* subtree;
	//subtree = treeLbData->CopyTree(cutAnalLb.getCut().c_str());
	//double bins[] = {0,1e-12,2e-12,5e-12,100e-12};
	std::vector<double> bins = variableBinSizeVec(-1, 0, .5e-12, 1e-12, 1.5e-12, 2e-12, 3e-12, 1);
	doEfficiencyPlotFitterLb(treeLbData, "HLT efficiency Displ/Barrel", "ct3dlb", bins, cutAnalLb.getCut()+"&&HLTokBarrelJpsi==1", "HLTokDisplJpsi==1", "HLTokDisplJpsi==0", canvaspager, htrep);
	doEfficiencyPlotFitterLb(treeLbData, "HLT efficiency Displ/Jpsi", "ct3dlb", bins, cutAnalLb.getCut()+"&&HLTokJpsi==1", "HLTokDisplJpsi==1", "HLTokDisplJpsi==0", canvaspager, htrep);
	doEfficiencyPlotFitterLb(treeLbData, "HLT efficiency Barrel/Jpsi", "ct3dlb", bins, cutAnalLb.getCut()+"&&HLTokJpsi==1", "HLTokBarrelJpsi==1", "HLTokBarrelJpsi==0", canvaspager, htrep);
    }

    if (doB0 && doEfficiencyPlotFitterPlots)
    {
	if (flgDoHtmlReport) htrep->addH("Trigger efficiencies plots for B<sup>0</sup>",'2');
	//TTree* subtree;
	//subtree = treeLbData->CopyTree(cutAnalLb.getCut().c_str());
	//double bins[] = {0,1e-12,2e-12,5e-12,100e-12};
	std::vector<double> bins = variableBinSizeVec(-10e-12, 0, .5e-12, 1e-12, 1.5e-12, 2e-12, 3e-12, 10e-12);
	doEfficiencyPlotFitterB0(treeB0Data, "HLT efficiency Displ/Barrel", "Lifetime", "ps", "ct3dB0", bins, cutAnalB0.getCut()+"&&HLTokBarrelJpsi==1", "HLTokDisplJpsi==1", "HLTokDisplJpsi==0", canvaspager, htrep);
	doEfficiencyPlotFitterB0(treeB0Data, "HLT efficiency Displ/Jpsi", "Lifetime", "ps", "ct3dB0", bins, cutAnalB0.getCut()+"&&HLTokJpsi==1", "HLTokDisplJpsi==1", "HLTokDisplJpsi==0", canvaspager, htrep);
	doEfficiencyPlotFitterB0(treeB0Data, "HLT efficiency Barrel/Jpsi", "Lifetime", "ps", "ct3dB0", bins, cutAnalB0.getCut()+"&&HLTokJpsi==1", "HLTokBarrelJpsi==1", "HLTokBarrelJpsi==0", canvaspager, htrep);
	doEfficiencyPlotFitterB0(treeB0Data, "HLT efficiency Displ/SingleMu", "Lifetime", "ps", "ct3dB0", bins, cutAnalB0.getCut()+"&&HLTSingleMu==1", "HLTokDisplJpsi==1", "HLTokDisplJpsi==0", canvaspager, htrep);
	doEfficiencyPlotFitterB0(treeB0Data, "HLT efficiency Barrel/SingleMu", "Lifetime", "ps", "ct3dB0", bins, cutAnalB0.getCut()+"&&HLTSingleMu==1", "HLTokBarrelJpsi==1", "HLTokBarrelJpsi==0", canvaspager, htrep);
	doEfficiencyPlotFitterB0(treeB0Data, "HLT efficiency Jpsi/SingleMu", "Lifetime", "ps", "ct3dB0", bins, cutAnalB0.getCut()+"&&HLTSingleMu==1", "HLTokJpsi==1", "HLTokJpsi==0", canvaspager, htrep);
    }

    // ===================================================================================
    if (doB0 && doLifetimeFit)
    {
	cout << "Lifetime fits for B0" << endl;
	canvaspager.cdNext();

	// Apply cut selection
	cout << "Using cut: " << cutAnalB0.getCut() << endl;
	//subtree = treeB0Data->CopyTree(cutAnalB0.getCut().c_str());
	std::auto_ptr<TTree> subtree(treeB0Data->CopyTree((cutAnalB0.getCut()+"&&HLTokBarrelJpsi==1&&mB0>5.22&&mB0<5.36").c_str()));
	cout << "subtree: " << subtree->GetEntries() << " entries of " << treeB0Data->GetEntries()<< endl;
	const double curEntries = subtree->GetEntries();
	cout << "curEntries: " << curEntries << endl;

	if (flgDoHtmlReport) htrep->addH("Lifetime plots B<sub>0</sup>",'2',true);
	doCtFitB0_fitresults res_CtFull = doCtFitB0_01(subtree.get(), filesLb.getLumiPbRounded(), "full dataset", flgNoTitle, flgPrelim, flgOfficialPlot);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> lifetime - full dataset (no trigger selection)");
    }

    // ===================================================================================
    if (doTPjpsi)
    {
	cout << "=======================================" << endl;
	cout << "Doing T&P efficiencies...." << endl;

	canvaspager.cdNext();

	doMassFitJp01_fitresults res;

	// Unfortunately the centrally produced T&P trees use a directory. I did not spend too much time to figure out how to use a TChain here, so one neds to hadd them into one file
	TFile * f = TFile::Open((path+"tptrees.root").c_str());
	//TFile * f = TFile::Open((path+"May10ReReco2011_00.root").c_str());
	TDirectory *dir = (TDirectory*) f->Get("tpTree");
	TTree *t = (TTree*) dir->Get("fitter_tree");

	//if (flgDoHtmlReport) htrep->addH("Mass plots J/&psi;",'2', true);
	//res = doMassFitJp01(t, 0, "full dataset", flgNoTitle, flgPrelim, flgOfficialPlot);
	//if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"J/&psi; mass plot - full dataset (no trigger selection)");

	if (flgDoHtmlReport) htrep->addH("Trigger efficiencies plots for J/&psi; using T&amp;P",'2');
	//TTree* subtree;
	//subtree = treeLbData->CopyTree(cutAnalLb.getCut().c_str());
	//double bins[] = {0,1e-12,2e-12,5e-12,100e-12};
	std::vector<double> bins;
	unsigned int toDo(4);
	if (toDo == 1)
	{
	    //bins.push_back(-0.1);
	    bins.push_back(0);
	    bins.push_back(0.005);
	    bins.push_back(0.01);
	    bins.push_back(0.02);
	    bins.push_back(0.03);
	    bins.push_back(0.04);
	    bins.push_back(0.05);
	    bins.push_back(0.07);
	    bins.push_back(0.10);
	    bins.push_back(0.15);
	    bins.push_back(0.25);
	    bins.push_back(0.50);
	    bins.push_back(1.00);
	    //bins.push_back(10.0);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon7JpsiDisplaced_Vtx/Dimuon7JpsiDisplaced_L3", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon7JpsiDisplaced_L3==1", "Dimuon7JpsiDisplaced_Vtx==1", "Dimuon7JpsiDisplaced_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon7JpsiDisplaced_Vtx/tag_Mu5_L2Mu2_Jpsi_MU", "Flight length", "cm", "pair_VtxLxy", bins, "tag_Mu5_L2Mu2_Jpsi_MU==1", "Dimuon7JpsiDisplaced_Vtx==1", "Dimuon7JpsiDisplaced_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon7JpsiDisplaced_Vtx/tag_Mu5_Track2_Jpsi_MU", "Flight length", "cm", "pair_VtxLxy", bins, "tag_Mu5_Track2_Jpsi_MU==1", "Dimuon7JpsiDisplaced_Vtx==1", "Dimuon7JpsiDisplaced_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon7JpsiDisplaced_Vtx/tag_Mu5_L2Mu2_Jpsi_MU HLT", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon7JpsiDisplaced_L3==1&&tag_Mu5_L2Mu2_Jpsi_MU==1", "Dimuon7JpsiDisplaced_Vtx==1", "Dimuon7JpsiDisplaced_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon7JpsiDisplaced_Vtx/tag_Mu5_Track2_Jpsi_MU HLT", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon7JpsiDisplaced_L3==1&&tag_Mu5_Track2_Jpsi_MU==1", "Dimuon7JpsiDisplaced_Vtx==1", "Dimuon7JpsiDisplaced_Vtx==0", canvaspager, htrep);
	}

	if (toDo == 2)
	{
	    //bins.push_back(-0.1);
	    bins.push_back(0);
	    bins.push_back(0.005);
	    bins.push_back(0.01);
	    bins.push_back(0.02);
	    bins.push_back(0.03);
	    bins.push_back(0.04);
	    bins.push_back(0.05);
	    bins.push_back(0.07);
	    bins.push_back(0.10);
	    bins.push_back(0.15);
	    bins.push_back(0.25);
	    bins.push_back(0.50);
	    bins.push_back(1.00);
	    //bins.push_back(10.0);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/Dimuon10JpsiBarrel_L3", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon10JpsiBarrel_L3==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/tag_Mu5_L2Mu2_Jpsi_MU", "Flight length", "cm", "pair_VtxLxy", bins, "tag_Mu5_L2Mu2_Jpsi_MU==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/tag_Mu5_Track2_Jpsi_MU", "Flight length", "cm", "pair_VtxLxy", bins, "tag_Mu5_Track2_Jpsi_MU==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/tag_Mu5_L2Mu2_Jpsi_MU HLT", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon10JpsiBarrel_L3==1&&tag_Mu5_L2Mu2_Jpsi_MU==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/tag_Mu5_Track2_Jpsi_MU HLT", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon10JpsiBarrel_L3==1&&tag_Mu5_Track2_Jpsi_MU==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	}

	if (toDo == 3)
	{
	    bins.push_back(0);
	    bins.push_back(1);
	    bins.push_back(1.5);
	    bins.push_back(2);
	    bins.push_back(2.5);
	    bins.push_back(3);
	    bins.push_back(3.5);
	    bins.push_back(4);
	    bins.push_back(4.5);
	    bins.push_back(5);
	    bins.push_back(6);
	    bins.push_back(7.5);
	    bins.push_back(9);
	    bins.push_back(12);
	    bins.push_back(16);
	    bins.push_back(30);
	    //bins.push_back(90);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon7JpsiDisplaced_Vtx/Dimuon7JpsiDisplaced_L3", "Flight length significance", "#sigma", "pair_VtxLxySig", bins, "Dimuon7JpsiDisplaced_L3==1", "Dimuon7JpsiDisplaced_Vtx==1", "Dimuon7JpsiDisplaced_Vtx==0", canvaspager, htrep);
	}

	if (toDo == 4)
	{
	    //bins.push_back(-0.1);
	    bins.push_back(0);
	    bins.push_back(0.005);
	    bins.push_back(0.01);
	    bins.push_back(0.02);
	    bins.push_back(0.03);
	    bins.push_back(0.04);
	    bins.push_back(0.05);
	    bins.push_back(0.07);
	    bins.push_back(0.10);
	    bins.push_back(0.15);
	    bins.push_back(0.25);
	    bins.push_back(0.50);
	    bins.push_back(1.00);
	    //bins.push_back(10.0);
	    //doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/Dimuon10JpsiBarrel_L3", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon10JpsiBarrel_L3==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/Mu5_L2Mu2_Jpsi_L2Mu", "Flight length", "cm", "pair_VtxLxy", bins, "Mu5_L2Mu2_Jpsi_L2Mu==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/Mu5_Track2_Jpsi_TK", "Flight length", "cm", "pair_VtxLxy", bins, "Mu5_Track2_Jpsi_TK==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/Mu5_L2Mu2_Jpsi_L2Mu HLT", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon10JpsiBarrel_L3==1&&Mu5_L2Mu2_Jpsi_L2Mu==1", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	    doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon10JpsiBarrel_Vtx/Mu5_Track2_Jpsi_TK HLT", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon10JpsiBarrel_L3==1&&Mu5_Track2_Jpsi_TK", "Dimuon10JpsiBarrel_Vtx==1", "Dimuon10JpsiBarrel_Vtx==0", canvaspager, htrep);
	}

	// Nachstehende funktionieren nicht
	// doEfficiencyPlotFitterJp(t, "HLT efficiency Displ/Barrel", "Flight length", "cm", "pair_VtxLxy", bins, "tag_Mu3_Track3_Jpsi_MU==1", "Dimuon7JpsiDisplaced_L3==1", "Dimuon7JpsiDisplaced_L3==0", canvaspager, htrep);
	// doEfficiencyPlotFitterJp(t, "HLT efficiency Displ/Barrel", "Flight length", "cm", "pair_VtxLxy", bins, "Dimuon10JpsiBarrel_L3==1", "Dimuon7JpsiDisplaced_L3==1", "Dimuon7JpsiDisplaced_L3==0", canvaspager, htrep);
	// doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon7JpsiDisplaced_L3/tag_Mu7_Track5_Jpsi_MU", "Flight length", "cm", "pair_VtxLxy", bins, "tag_Mu7_Track5_Jpsi_MU==1", "Dimuon7JpsiDisplaced_L3==1", "Dimuon7JpsiDisplaced_L3==0", canvaspager, htrep);
	// doEfficiencyPlotFitterJp(t, "HLT efficiency Dimuon7JpsiDisplaced_L3/Mu3_Track3_Jpsi_TK", "Flight length", "cm", "pair_VtxLxy", bins, "Mu3_Track3_Jpsi_TK==1", "Dimuon7JpsiDisplaced_L3==1", "Dimuon7JpsiDisplaced_L3==0", canvaspager, htrep);
    }

    // ===================================================================================
    if (doMCchannelPlots && doMC && doB0)
    {
	if (flgDoHtmlReport) htrep->addH("Mass plots for B<sup>0</sup> channels in MC",'2', true);
	stackedPlotList spl;
	spl.add("B^{0} #rightarrow others","isSig==0&&nRef2G!=1007",1);
	spl.add("B^{0} #rightarrow J/#psi (#mu#mu) K_{s} (#pi#pi)","isSig==1",2);
	spl.add("B^{0} #rightarrow J/#psi (#mu#mu#gamma) K_{s} (#pi#pi)","isSig==0&&nRef2G==1007",6);
	canvaspager.cdNext();
	doStackedPlot(treeB0MC, "MCB0sig", "B^{0} decay modes in MC", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 40,4.6,6,spl);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - full dataset (no cuts except channel selection)");
	canvaspager.cdNext();
	doStackedPlot(treeB0MC, "MCB0sigCut", "B^{0} decay modes in MC", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 40,4.6,6,spl, cutAnalB0.getCut());
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - full dataset (analysis cuts)");
	canvaspager.cdNext();
	doStackedPlot(treeB0MC, "MCB0siglog", "B^{0} decay modes in MC", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 40,4.6,6,spl,"",true);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - full dataset (no cuts except channel selection)");
	canvaspager.cdNext();
	doStackedPlot(treeB0MC, "MCB0sigCutlog", "B^{0} decay modes in MC", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 40,4.6,6,spl, cutAnalB0.getCut(),true);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - full dataset (analysis cuts)");

    }

    // ===================================================================================
    if (doBgrChannelPlots && doMC && doLb)
    {
	if (flgDoHtmlReport) htrep->addH("Mass plots for &Lambda;<sub>b</sub> background channels in MC",'2', true);
	stackedPlotList spl;
	Double_t normLumi = 5000; // pb-1
	spl.add(chainLbBgrMC_Xi, "#Xi_{b}", "", normLumi/filesLbBgrMC_Xi.getLumiPbRounded(), 7);
	spl.add(chainLbBgrMC_Om, "#Omega_{b}", "", normLumi/filesLbBgrMC_Om.getLumiPbRounded(), 8);
	spl.add(chainLbBgrMC_B0, "B^{0}", "", normLumi/filesLbBgrMC_B0.getLumiPbRounded(), 2);
	spl.add(chainLbBgrMC_Bp, "B^{+}", "", normLumi/filesLbBgrMC_Bp.getLumiPbRounded(), 3);
	spl.add(chainLbBgrMC_Bs, "B_{s}", "", normLumi/filesLbBgrMC_Bs.getLumiPbRounded(), 4);
	//spl.add(chainLbBgrMC_Jp, "J/#psi", "", normLumi/filesLbBgrMC_Jp.getLumiPbRounded(), 6);
	spl.add(chainLbBgrMC_Jp, "J/#psi", "", 1, 6);
	canvaspager.cdNext();
	doStackedPlot("MCLbBgr_mlb", "#Lambda_{b} backgrounds (MC)", "mlb", "m(#mu#mup#pi)", "GeV/c^{2}", 20,5.0,6.0,spl, cutAnalLb.getCut());
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - background channels");
	canvaspager.cdNext();
	doStackedPlot("MCLbBgr_ct3dlb", "#Lambda_{b} backgrounds (MC)", "ct3dlb", "t", "s", 20,-5e-12,15e-12,spl, cutAnalLb.getCut(), true);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"&Lambda;<sub>b</sub> mass - background channels");
    }

    // ===================================================================================
    if (doBgrChannelPlots && doMC && doB0)
    {
	if (flgDoHtmlReport) htrep->addH("Mass plots for B<sup>0</sup> background channels in MC",'2', true);
	stackedPlotList spl;
	Double_t normLumi = 3000; // pb-1
	spl.add(chainB0BgrMC_Lb, "#Lambda_{b}", "", normLumi/filesB0BgrMC_Lb.getLumiPbRounded(), 2);
	//spl.add(chainB0BgrMC_Jp, "J/#psi", "", normLumi/filesB0BgrMC_Jp.getLumiPbRounded(), 6);
	spl.add(chainB0BgrMC_Jp, "J/#psi", "", .1*normLumi/filesB0BgrMC_Jp.getLumiPbRounded(), 6);
	spl.add(chainB0BgrMC_Bp, "B^{+}", "", normLumi/filesB0BgrMC_Bp.getLumiPbRounded(), 3);
	spl.add(chainB0BgrMC_Bs, "B_{s}", "", normLumi/filesB0BgrMC_Bs.getLumiPbRounded(), 4);
	spl.add(chainB0MC, "B^{0} signal", "", 1, 7);
	canvaspager.cdNext();
	doStackedPlot("MCB0Bgr_mlb", "B^{0} backgrounds (MC)", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 200,5.0,6.0,spl, cutAnalB0MC.getCut());
	//doStackedPlot("MCB0Bgr_mlb", "B^{0} backgrounds (MC)", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 20,5.0,6.0,spl);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - background channels");
	canvaspager.cdNext();

	doStackedPlot("MCB0Bgr_mlb", "B^{0} backgrounds (MC) L0veto", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 200,5.0,6.0,spl, cutAnalB0MC.getCut()+"&&(L0hypo<1.1077||L0hypo>1.1237)");
	//doStackedPlot("MCB0Bgr_mlb", "B^{0} backgrounds (MC)", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 20,5.0,6.0,spl);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - background channels");
	canvaspager.cdNext();
	//doStackedPlot("MCB0Bgr_ct3dlb", "B^{0} backgrounds (MC)", "ct3dB0", "t", "s", 20,-5e-12,15e-12,spl, cutAnalB0MC.getCut(), true);
	doStackedPlot("MCB0Bgr_ct3dlb", "B^{0} backgrounds (MC)", "ct3dB0", "t", "s", 50,-5e-12,15e-12,spl, cutAnalB0MC.getCut(), true);
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"B<sup>0</sup> mass - background channels");
    }

    // ===================================================================================
    // Plots for lifetime reconstruction sanity
    if (doB0tTruthPlots && doB0 && doMC)
    {
	const string obj2plot("B0");
	const string obj4html("B<sup>0</sup>");
	const string obj4root("B^{0}");

	const bool noTitle(true);

	if (flgDoHtmlReport) htrep->addH(obj4html+" - Monte-Carlo truth comparison plots for lifetime",'2',true);

	const double ctWindow(1), ctWindowZoom(2);
	const double d3dWindow(2000), d3dWindowZoom(4000);
	const double pWindow(1), pWindowZoom(2);
	const double massB0(5.2794), massLo(massB0-.1), massHi(massB0+.1);

	canvaspager.cdNext("MCtruthCt3d"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthCt3d"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth)", "", 200, -ctWindow, ctWindow, noTitle ? "" : "Lifetime of "+obj4root+" 3d, all candidates found", "t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime of "+obj4html+" 3d, all candidates found");
	canvaspager.cdNext("MCtruthCt3dTruth"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthCt3dTruth"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth)", "isMCmatch==1", 200, -ctWindow, ctWindow, noTitle ? "" : "Lifetime of "+obj4root+" 3d, truth matched candidates", "t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime of "+obj4html+" 3d, truth matched candidates");
	canvaspager.cdNext("MCtruthCt3dCutTruth"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthCt3dCutTruth"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth)", cutAnalB0MC.getCut()+"&&isMCmatch==1", 100, -ctWindow, ctWindow, noTitle ? "" : "Lifetime of "+obj4root+" 3d, truth matched candidates, cuts", "t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime of "+obj4html+" 3d, truth matched candidates, cuts");

	if (flgDoHtmlReport) htrep->addH("Flight length",'2');
	canvaspager.cdNext("MCtruthd3d"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthd3d"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth)", "", 200, -d3dWindow, d3dWindow, noTitle ? "" : "Flight length of "+obj4root+" 3d, all candidates found", "d_{reco,3d}("+obj4root+")-d_{truth}("+obj4root+")", "#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Flight length of "+obj4html+" 3d, all candidates found");
	canvaspager.cdNext("MCtruthd3dTruth"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthd3dTruth"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth)", "isMCmatch==1", 200, -d3dWindow, d3dWindow, noTitle ? "" : "Flight length of "+obj4root+" 3d, truth matched candidates", "d_{reco,3d}("+obj4root+")-d_{truth}("+obj4root+")", "#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Flight length of "+obj4html+" 3d, truth matched candidates");
	canvaspager.cdNext("MCtruthd3dCutTruth"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthd3dCutTruth"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth)", cutAnalB0MC.getCut()+"&&isMCmatch==1", 100, -d3dWindow, d3dWindow, noTitle ? "" : "Flight length of "+obj4root+" 3d, truth matched candidates, cuts", "d_{reco,3d}("+obj4root+")-d_{truth}("+obj4root+")", "#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Flight length of "+obj4html+" 3d, truth matched candidates, cuts");

	if (flgDoHtmlReport) htrep->addH("Momentum",'2');
	canvaspager.cdNext("MCtruthp"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthp"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth", "", 200, -pWindow, pWindow, noTitle ? "" : "Momentum of "+obj4root+" all candidates found", "p_{reco}("+obj4root+")-p_{truth}("+obj4root+")", "GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Momentum of "+obj4html+" all candidates found");
	canvaspager.cdNext("MCtruthpTruth"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthpTruth"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth", "isMCmatch==1", 200, -pWindow, pWindow, noTitle ? "" : "Momentum of "+obj4root+" truth matched candidates", "p_{reco}("+obj4root+")-p_{truth}("+obj4root+")", "GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Momentum of "+obj4html+" truth matched candidates");
	canvaspager.cdNext("MCtruthpCutTruth"+obj2plot);
	do1dPlotGaus(treeB0MC, "MCtruthpCutTruth"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth", cutAnalB0MC.getCut()+"&&isMCmatch==1", 100, -pWindow, pWindow, noTitle ? "" : "Momentum of "+obj4root+" truth matched candidates, cuts", "p_{reco}("+obj4root+")-p_{truth}("+obj4root+")", "GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Momentum of "+obj4html+" truth matched candidates, cuts");


	// ---------------------------------------------------
	if (flgDoHtmlReport) htrep->addH("Profile plots for 3d measurement",'2');

	canvaspager.cdNext("TtruthVsT3dprof"+obj2plot);
	do2dProfilePlot(treeB0MC, "TtruthVsT3dprof"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth):1e12*(ct"+obj2plot+"truth)", "isMCmatch==1", 32, 0, 16, -ctWindow, +ctWindow, noTitle ? "" : "t_{3d}-t_truth vs. t","t_{truth}","ps", "#Deltat","ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"t_<sub>3d</sub>-t_truth vs. t, truth matched");

	canvaspager.cdNext("TtruthVsT3dprof_cuts"+obj2plot);
	do2dProfilePlot(treeB0MC, "TtruthVsT3dprof_cuts"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth):1e12*(ct"+obj2plot+"truth)", cutAnalB0MC.getCut()+"&&isMCmatch==1", 32, 0, 16, -ctWindow, +ctWindow, noTitle ? "" : "t_{3d}-t_truth vs. t","t_{truth}","ps", "#Deltat","ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"t_<sub>3d</sub>-t_truth vs. t, cuts");

	if (flgDoHtmlReport) htrep->addH("Flight length",'2');
	canvaspager.cdNext("DtruthVsD3dprof"+obj2plot);
	do2dProfilePlot(treeB0MC, "DtruthVsD3dprof"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth):10000*(d3d"+obj2plot+"truth)", "isMCmatch==1", 32, 0, 2*d3dWindow, -d3dWindow, +d3dWindow, noTitle ? "" : "d_{3d}-d_truth vs. d","d_{truth}","#mum", "#Deltad","#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"d_<sub>3d</sub>-d<sub>truth</sub> vs. d, truth matched");

	canvaspager.cdNext("DtruthVsD3dprof_cuts"+obj2plot);
	do2dProfilePlot(treeB0MC, "DtruthVsD3dprof_cuts"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth):10000*(d3d"+obj2plot+"truth)", cutAnalB0MC.getCut()+"&&isMCmatch==1", 32, 0, 2*d3dWindow, -d3dWindow, +d3dWindow, noTitle ? "" : "d_{3d}-d_truth vs. d","d_{truth}","#mum", "#Deltat","#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"d_<sub>3d</sub>-d<sub>truth</sub> vs. d, cuts");

	if (flgDoHtmlReport) htrep->addH("Momentum",'2');
	canvaspager.cdNext("PtruthVsP3dprof"+obj2plot);
	do2dProfilePlot(treeB0MC, "PtruthVsP3dprof"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth:p"+obj2plot+"truth", "isMCmatch==1", 25, 10, 100, -4*pWindow, +4*pWindow, noTitle ? "" : "p-p_truth vs. p","p_{truth}","GeV/c", "#Deltap","GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"p_-p<sub>truth</sub> vs. p, truth matched");

	canvaspager.cdNext("PtruthVsP3dprof_cuts"+obj2plot);
	do2dProfilePlot(treeB0MC, "PtruthVsP3dprof_cuts"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth:p"+obj2plot+"truth", cutAnalB0MC.getCut()+"&&isMCmatch==1", 25, 10, 100, -4*pWindow, +4*pWindow, noTitle ? "" : "p-p_truth vs. p","p_{truth}","GeV/c", "#Deltap","GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"p-p<sub>truth</sub> vs. p, cuts");

    }

    if (doLbtTruthPlots && doLb && doMC)
    {
	const string obj2plot("lb");
	const string obj4html("&Lambda;<sub>b</sub>");
	const string obj4root("#Lambda_{b}");

	const bool noTitle(true);

	if (flgDoHtmlReport) htrep->addH(obj4html+" - Monte-Carlo truth comparison plots for lifetime",'2',true);

	const double ctWindow(1), ctWindowZoom(2);
	const double d3dWindow(2000), d3dWindowZoom(4000);
	const double pWindow(1), pWindowZoom(2);
	const double massLb(5.2794), massLo(massLb-.1), massHi(massLb+.1);

	canvaspager.cdNext("MCtruthCt3d"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthCt3d"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth)", "", 200, -ctWindow, ctWindow, noTitle ? "" : "Lifetime of "+obj4root+" 3d, all candidates found", "t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime of "+obj4html+" 3d, all candidates found");
	canvaspager.cdNext("MCtruthCt3dTruth"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthCt3dTruth"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth)", "isMCmatch==1", 200, -ctWindow, ctWindow, noTitle ? "" : "Lifetime of "+obj4root+" 3d, truth matched candidates", "t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime of "+obj4html+" 3d, truth matched candidates");
	canvaspager.cdNext("MCtruthCt3dCutTruth"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthCt3dCutTruth"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth)", cutAnalLbMC.getCut()+"&&isMCmatch==1", 100, -ctWindow, ctWindow, noTitle ? "" : "Lifetime of "+obj4root+" 3d, truth matched candidates, cuts", "t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime of "+obj4html+" 3d, truth matched candidates, cuts");

	canvaspager.cdNext("MCtruthCt3dPull"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthCt3dPull"+obj2plot, "(ct3d"+obj2plot+"-ct"+obj2plot+"truth)/ct3d"+obj2plot+"E", "", 200, -5, 5, noTitle ? "" : "Lifetime pull of "+obj4root+" 3d, all candidates found", "pull t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "#sigma");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime pulls of "+obj4html+" 3d, all candidates found");
	canvaspager.cdNext("MCtruthCt3dPullTruth"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthCt3dPullTruth"+obj2plot, "(ct3d"+obj2plot+"-ct"+obj2plot+"truth)/ct3d"+obj2plot+"E", "isMCmatch==1", 200, -5, 5, noTitle ? "" : "Lifetime pulls of "+obj4root+" 3d, truth matched candidates", "t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime pulls of "+obj4html+" 3d, truth matched candidates");
	canvaspager.cdNext("MCtruthCt3dPullCutTruth"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthCt3dPullCutTruth"+obj2plot, "(ct3d"+obj2plot+"-ct"+obj2plot+"truth)/ct3d"+obj2plot+"E", cutAnalLbMC.getCut()+"&&isMCmatch==1", 100, -5, 5, noTitle ? "" : "Lifetime pulls of "+obj4root+" 3d, truth matched candidates, cuts", "t_{reco,3d}("+obj4root+")-t_{truth}("+obj4root+")", "ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Lifetime pulls of "+obj4html+" 3d, truth matched candidates, cuts");

	//-----
	if (flgDoHtmlReport) htrep->addH("Flight length",'2');
	canvaspager.cdNext("MCtruthd3d"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthd3d"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth)", "", 200, -d3dWindow, d3dWindow, noTitle ? "" : "Flight length of "+obj4root+" 3d, all candidates found", "d_{reco,3d}("+obj4root+")-d_{truth}("+obj4root+")", "#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Flight length of "+obj4html+" 3d, all candidates found");
	canvaspager.cdNext("MCtruthd3dTruth"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthd3dTruth"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth)", "isMCmatch==1", 200, -d3dWindow, d3dWindow, noTitle ? "" : "Flight length of "+obj4root+" 3d, truth matched candidates", "d_{reco,3d}("+obj4root+")-d_{truth}("+obj4root+")", "#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Flight length of "+obj4html+" 3d, truth matched candidates");
	canvaspager.cdNext("MCtruthd3dCutTruth"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthd3dCutTruth"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth)", cutAnalLbMC.getCut()+"&&isMCmatch==1", 100, -d3dWindow, d3dWindow, noTitle ? "" : "Flight length of "+obj4root+" 3d, truth matched candidates, cuts", "d_{reco,3d}("+obj4root+")-d_{truth}("+obj4root+")", "#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Flight length of "+obj4html+" 3d, truth matched candidates, cuts");

	//-----
	if (flgDoHtmlReport) htrep->addH("Momentum",'2');
	canvaspager.cdNext("MCtruthp"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthp"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth", "", 200, -pWindow, pWindow, noTitle ? "" : "Momentum of "+obj4root+" all candidates found", "p_{reco}("+obj4root+")-p_{truth}("+obj4root+")", "GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Momentum of "+obj4html+" all candidates found");
	canvaspager.cdNext("MCtruthpTruth"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthpTruth"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth", "isMCmatch==1", 200, -pWindow, pWindow, noTitle ? "" : "Momentum of "+obj4root+" truth matched candidates", "p_{reco}("+obj4root+")-p_{truth}("+obj4root+")", "GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Momentum of "+obj4html+" truth matched candidates");
	canvaspager.cdNext("MCtruthpCutTruth"+obj2plot);
	do1dPlotGaus(treeLbMC, "MCtruthpCutTruth"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth", cutAnalLbMC.getCut()+"&&isMCmatch==1", 100, -pWindow, pWindow, noTitle ? "" : "Momentum of "+obj4root+" truth matched candidates, cuts", "p_{reco}("+obj4root+")-p_{truth}("+obj4root+")", "GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(), "Momentum of "+obj4html+" truth matched candidates, cuts");


	// ---------------------------------------------------
	if (flgDoHtmlReport) htrep->addH("Profile plots for 3d measurement",'2');

	canvaspager.cdNext("TtruthVsT3dprof"+obj2plot);
	do2dProfilePlot(treeLbMC, "TtruthVsT3dprof"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth):1e12*(ct"+obj2plot+"truth)", "isMCmatch==1", 32, 0, 16, -ctWindow, +ctWindow, noTitle ? "" : "t_{3d}-t_truth vs. t","t_{truth}","ps", "#Deltat","ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"t_<sub>3d</sub>-t_truth vs. t, truth matched");

	canvaspager.cdNext("TtruthVsT3dprof_cuts"+obj2plot);
	do2dProfilePlot(treeLbMC, "TtruthVsT3dprof_cuts"+obj2plot, "1e12*(ct3d"+obj2plot+"-ct"+obj2plot+"truth):1e12*(ct"+obj2plot+"truth)", cutAnalLbMC.getCut()+"&&isMCmatch==1", 32, 0, 16, -ctWindow, +ctWindow, noTitle ? "" : "t_{3d}-t_truth vs. t","t_{truth}","ps", "#Deltat","ps");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"t_<sub>3d</sub>-t_truth vs. t, cuts");

	if (flgDoHtmlReport) htrep->addH("Flight length",'2');
	canvaspager.cdNext("DtruthVsD3dprof"+obj2plot);
	do2dProfilePlot(treeLbMC, "DtruthVsD3dprof"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth):10000*(d3d"+obj2plot+"truth)", "isMCmatch==1", 32, 0, 2*d3dWindow, -d3dWindow, +d3dWindow, noTitle ? "" : "d_{3d}-d_truth vs. d","d_{truth}","#mum", "#Deltad","#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"d_<sub>3d</sub>-d<sub>truth</sub> vs. d, truth matched");

	canvaspager.cdNext("DtruthVsD3dprof_cuts"+obj2plot);
	do2dProfilePlot(treeLbMC, "DtruthVsD3dprof_cuts"+obj2plot, "10000*(d3"+obj2plot+"-d3d"+obj2plot+"truth):10000*(d3d"+obj2plot+"truth)", cutAnalLbMC.getCut()+"&&isMCmatch==1", 32, 0, 2*d3dWindow, -d3dWindow, +d3dWindow, noTitle ? "" : "d_{3d}-d_truth vs. d","d_{truth}","#mum", "#Deltat","#mum");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"d_<sub>3d</sub>-d<sub>truth</sub> vs. d, cuts");

	if (flgDoHtmlReport) htrep->addH("Momentum",'2');
	canvaspager.cdNext("PtruthVsP3dprof"+obj2plot);
	do2dProfilePlot(treeLbMC, "PtruthVsP3dprof"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth:p"+obj2plot+"truth", "isMCmatch==1", 25, 10, 100, -4*pWindow, +4*pWindow, noTitle ? "" : "p-p_truth vs. p","p_{truth}","GeV/c", "#Deltap","GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"p_-p<sub>truth</sub> vs. p, truth matched");

	canvaspager.cdNext("PtruthVsP3dprof_cuts"+obj2plot);
	do2dProfilePlot(treeLbMC, "PtruthVsP3dprof_cuts"+obj2plot, "p"+obj2plot+"-p"+obj2plot+"truth:p"+obj2plot+"truth", cutAnalLbMC.getCut()+"&&isMCmatch==1", 25, 10, 100, -4*pWindow, +4*pWindow, noTitle ? "" : "p-p_truth vs. p","p_{truth}","GeV/c", "#Deltap","GeV/c");
	if (flgDoHtmlReport) htrep->addTableImage(canvaspager.getCurPng(),"p-p<sub>truth</sub> vs. p, cuts");

    }

    // ===================================================================================
    // finalize page
    canvaspager.forceSave();
    if (flgDoHtmlReport) delete htrep;

    //delete canvas;
}

int main()
{
    doAnalysis01("Test",true);
    return 0;
}


