#include <string>
#include <vector>
#include "TTree.h"
#include "TFile.h"
//#include "doMassFitLb01.C"
#include "HtmlReport.h"
#include "Canvaspager.h"
#include <iostream>

using std::string;
using std::vector;
using std::cout;
using std::endl;

void doEfficiencyPlotFitterLb(TTree *tree, const string &title, const string &binVar, const std::vector<double> &bins, const string &cut, const string &cutPass, const string &cutFail, Canvaspager &cp, HtmlReport *htrep)
{
    //TFile* f = new TFile("doEffPlot_tmp.root","RECREATE");

    bool flgDoHtmlReport = (htrep != 0);

    vector<doMassFitLb01_fitresults> resFull, resPass, resFail;

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
	doMassFitLb01_fitresults res = doMassFitLb01(subtree, 0, curTitle, false, false, false);
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
	doMassFitLb01_fitresults res = doMassFitLb01(subtree, 0, curTitle, false, false, false);
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
	doMassFitLb01_fitresults res = doMassFitLb01(subtree, 0, curTitle, false, false, false);
	resFail.push_back(res);
	if (flgDoHtmlReport) htrep->addTableImage(cp.getCurPng(),curTitle);
    }
    //delete f;

    cout << "no. full pass fail" << endl;
    if (flgDoHtmlReport) htrep->beginTable();
    double sumFull(0), sumPass(0), sumFail(0);
    for (int i = 0; i!=nBins; i++)
    {
	sumFull += resFull[i].sig;
	sumPass += resPass[i].sig;
	sumFail += resFail[i].sig;
	cout << i << ": " << resFull[i].sig << " " << resPass[i].sig << " " << resFail[i].sig << endl;
	//if (flgDoHtmlReport) htrep->addTableRow(toString(i)+": ("+toString(bins[i])+ "," + toString(bins[i+1]) + "]",roundToString(resFull[i].sig,1),roundToString(resPass[i].sig,1),roundToString(resFail[i].sig,1));
	if (flgDoHtmlReport)
	{
	    htrep->beginTableRow();
	    htrep->addTableCell(toString(i)+": (");
	    htrep->addTableCell(toString(bins[i]));
	    htrep->addTableCell(" , "+toString(bins[i+1])+"]");
	    htrep->addTableCell(roundToString(resFull[i].sig,1));
	    htrep->addTableCell(roundToString(resPass[i].sig,1));
	    htrep->addTableCell(roundToString(resPass[i].sig/resFull[i].sig*100,1)+"%");
	    htrep->addTableCell(roundToString(resFail[i].sig,1));
	    htrep->endTableRow();
	}
    }
    if (flgDoHtmlReport) htrep->addTableRow("Sum:", roundToString(sumFull,1), roundToString(sumPass,1), roundToString(sumFail,1));
    if (flgDoHtmlReport) htrep->endTable();
}

