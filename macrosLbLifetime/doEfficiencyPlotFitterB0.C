#include <string>
#include <vector>
#include "TTree.h"
#include "TFile.h"
//#include "doMassFitB001.C"
#include "HtmlReport.h"
#include "Canvaspager.h"
#include "doGraphError.C"
#include <iostream>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void doEfficiencyPlotFitterB0(TTree *tree, const string &title, const string &titleObs, const string &titleUnitObs,
	const string &binVar, const std::vector<double> &bins, const string &cut, const string &cutPass, const string &cutFail, Canvaspager &cp, HtmlReport *htrep)
{
    //TFile* f = new TFile("doEffPlot_tmp.root","RECREATE");

    bool flgDoHtmlReport = (htrep != 0);
    if (flgDoHtmlReport) htrep->addH(title,3);

    vector<doMassFitB001_fitresults> resFull, resPass, resFail;

    const int nBins = bins.size()-1;

    // fits for all, i.e. no additional requirement except bin
    if (flgDoHtmlReport) htrep->addP("All (=100%)");
    if (flgDoHtmlReport) htrep->addP("Cut: " + cut);
    for (int i = 0; i!=nBins; i++)
    {
	const string curTitle = title + " bin " + toString(i) + " (" + toString(bins[i]) + "," + toString(bins[i+1]) + "]";
	const string curCut = cut + "&&" + binVar + ">" + toString(bins[i]) + "&&" + binVar + "<=" + toString(bins[i+1]);
	TTree *subtree = tree->CopyTree(curCut.c_str());
	cp.cdNext();
	doMassFitB001_fitresults res = doMassFitB001(subtree, 0, curTitle, false, false, false);
	resFull.push_back(res);
	if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),curTitle);
    }
    // fits for pass
    if (flgDoHtmlReport) htrep->addP("Pass");
    if (flgDoHtmlReport) htrep->addP("Cut: " + cut + "&&" + cutPass);
    for (int i = 0; i!=nBins; i++)
    {
	const string curTitle = title + " bin " + toString(i) + " (" + toString(bins[i]) + "," + toString(bins[i+1]) + "]";
	const string curCut = cut + "&&" + cutPass + "&&" + binVar + ">" + toString(bins[i]) + "&&" + binVar + "<=" + toString(bins[i+1]);
	TTree *subtree = tree->CopyTree(curCut.c_str());
	cp.cdNext();
	doMassFitB001_fitresults res = doMassFitB001(subtree, 0, curTitle, false, false, false);
	resPass.push_back(res);
	if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),curTitle);
    }
    // fits for fail
    if (flgDoHtmlReport) htrep->addP("Fail");
    if (flgDoHtmlReport) htrep->addP("Cut: " + cut + "&&" + cutFail);
    for (int i = 0; i!=nBins; i++)
    {
	const string curTitle = title + " bin " + toString(i) + " (" + toString(bins[i]) + "," + toString(bins[i+1]) + "]";
	const string curCut = cut + "&&" + cutFail + "&&" + binVar + ">" + toString(bins[i]) + "&&" + binVar + "<=" + toString(bins[i+1]);
	TTree *subtree = tree->CopyTree(curCut.c_str());
	cp.cdNext();
	doMassFitB001_fitresults res = doMassFitB001(subtree, 0, curTitle, false, false, false);
	resFail.push_back(res);
	if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),curTitle);
    }
    //delete f;

    cout << "no. full pass fail" << endl;
    if (flgDoHtmlReport) {
	htrep->beginTable();
	htrep->beginTableRow();
	htrep->addTableCell("i");
	htrep->addTableCell("");
	htrep->addTableCell("bin");
	htrep->addTableCell("tot");
	htrep->addTableCell("pass");
	htrep->addTableCell("efficiency");
	htrep->addTableCell("fail");
	htrep->endTableRow();
    }
    double sumFull(0), sumPass(0), sumFail(0);
    std::vector<double> xVec, effVec, eff_errVec;
    for (int i = 0; i!=nBins; i++)
    {
	sumFull += resFull[i].sig;
	sumPass += resPass[i].sig;
	sumFail += resFail[i].sig;
	cout << i << ": " << resFull[i].sig << " " << resPass[i].sig << " " << resFail[i].sig << endl;
	//if (flgDoHtmlReport) htrep->addTableRow(toString(i)+": ("+toString(bins[i])+ "," + toString(bins[i+1]) + "]",roundToString(resFull[i].sig,1),roundToString(resPass[i].sig,1),roundToString(resFail[i].sig,1));
	const double Ntot = resFull[i].sig;
	const double Ntot_err = resFull[i].sig_err;
	const double Npass = resPass[i].sig;
	const double Npass_err = resPass[i].sig_err;
	const double eff = Npass/Ntot;
	double eff_err = sqrt(1/(Ntot*Ntot)*(Npass_err*Npass_err+eff*eff*Ntot_err*Ntot_err));
	xVec.push_back((bins[i]+bins[i+1])/2);
	effVec.push_back(eff);
	eff_errVec.push_back(eff_err);
	if (flgDoHtmlReport)
	{
	    htrep->beginTableRow();
	    htrep->addTableCell(toString(i)+": (");
	    htrep->addTableCell(toString(bins[i]));
	    htrep->addTableCell(" , "+toString(bins[i+1])+"]");
	    htrep->addTableCell(roundToString(resFull[i].sig,1)+" &plusmn; "+roundToString(resFull[i].sig_err,1));
	    htrep->addTableCell(roundToString(resPass[i].sig,1)+" &plusmn; "+roundToString(resPass[i].sig_err,1));
	    htrep->addTableCell(roundToString(eff,3)+" &plusmn; "+roundToString(eff_err,3));
	    htrep->addTableCell(roundToString(resFail[i].sig,1)+" &plusmn; "+roundToString(resFail[i].sig_err,1));
	    htrep->endTableRow();
	}
    }
    if (flgDoHtmlReport) htrep->addTableRow("Sum:", roundToString(sumFull,1), roundToString(sumPass,1), roundToString(sumFail,1));
    if (flgDoHtmlReport) htrep->endTable();
    // Make plot
    cp.cdNext();
    doGraphError(xVec, effVec, eff_errVec, title, titleObs, titleUnitObs, "Efficiency", "", 0);
    if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),"Efficiency plot");
}

