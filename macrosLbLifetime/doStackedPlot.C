#include <string>
#include <vector>
#include <iostream>

#include "THStack.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TPad.h"

#include "utils.h"

using std::string;
using std::cout;
using std::endl;

struct stackedPlotDef
{
    stackedPlotDef() : title(""), cut(""), scale(1.0), color(1) {};
    stackedPlotDef(string t, string c, int col) : title(t), cut(c), scale(1.0), color(col) {};
    stackedPlotDef(string t, string c, Double_t sc, int col) : title(t), cut(c), scale(sc), color(col) {};
    stackedPlotDef(TTree *tr, string t, string c, int col) : tree(tr), title(t), cut(c), scale(1.0), color(col) {};
    stackedPlotDef(TTree *tr, string t, string c, Double_t sc, int col) : tree(tr), title(t), cut(c), scale(sc), color(col) {};
    TTree *tree;
    string title;
    string cut;
    Double_t scale;
    int color;
};

struct stackedPlotList
{
    typedef std::vector<stackedPlotDef> stackedPlotListVec_t;
    void add(const string &title, const string &cut, const int &color)
    {
	vec.push_back(stackedPlotDef(title, cut, color));
    }
    void add(TTree *tr, const string &title, const string &cut, const Double_t sc, const int &color)
    {
	vec.push_back(stackedPlotDef(tr, title, cut, sc, color));
    }
    stackedPlotListVec_t vec;

    string getAllCutsORed()
    {
	string ret("");
	int i(0);
	for(stackedPlotListVec_t::const_iterator it = vec.begin(); it != vec.end(); it++, i++)
	{
	    if (i!=0) ret += "||";
	    ret += it->cut;
	}
	return ret;
    }
};

void doStackedPlot(TTree *tree, string name, string title, string toDraw, string xtitle, string xunit, int nBins, double xmin, double xmax, stackedPlotList spl, string cut = "", bool setLog = false)
{
    THStack *hs = new THStack((name + "_stack").c_str(), title.c_str());

    const float legX(0.55), legXwidth(0.40);
    const float legY(0.9), legYwidth1(0.05);
    TLegend *legend = new TLegend(legX,legY-legYwidth1*spl.vec.size(),legX+legXwidth,legY);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);

    const string curYtitle = entriesPerBin(nBins, xmax, xmin, xunit);
    const string curXtitle = valueWithUnit(xtitle, xunit);

    int i(0);
    for(stackedPlotList::stackedPlotListVec_t::const_iterator it = spl.vec.begin(); it != spl.vec.end(); it++, i++)
    {
	const string curName = "h" + name + toString(i);
	const string drawstring = toDraw + ">>" + curName + "(" + toString(nBins) + "," + toString(xmin) + "," + toString(xmax) + ")";
	const string curcut = ((cut.size() > 0) ? cut + "&&" : "") + "(" + it->cut + ")";
	cout << drawstring.c_str() << " - " << curcut << " scale: " << it->scale << endl;
	tree->Draw(drawstring.c_str(),curcut.c_str());
	TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(curName.c_str());
	h->SetFillColor(it->color);
	h->GetXaxis()->SetTitle(curXtitle.c_str());
	h->GetYaxis()->SetTitle(curYtitle.c_str());
	h->GetYaxis()->SetTitle(curYtitle.c_str());
	h->Scale(it->scale);
	legend->AddEntry(curName.c_str(),it->title.c_str(),"f");
	hs->Add(h);
    }
    hs->Draw();
    hs->GetXaxis()->SetTitle(curXtitle.c_str());
    hs->GetYaxis()->SetTitle(curYtitle.c_str());
    hs->GetYaxis()->SetTitleOffset(1.25);
    hs->SetTitle(title.c_str());
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.1);
    gPad->SetLogy(setLog);
    legend->Draw();
}

// case where source data comes from different trees
void doStackedPlot(string name, string title, string toDraw, string xtitle, string xunit, int nBins, double xmin, double xmax, stackedPlotList spl, string cut = "", bool setLog = false)
{
    THStack *hs = new THStack((name + "_stack").c_str(), title.c_str());

    const float legX(0.55), legXwidth(0.40);
    const float legY(0.9), legYwidth1(0.05);
    TLegend *legend = new TLegend(legX,legY-legYwidth1*spl.vec.size(),legX+legXwidth,legY);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);

    const string curYtitle = entriesPerBin(nBins, xmax, xmin, xunit);
    const string curXtitle = valueWithUnit(xtitle, xunit);

    int i(0);
    for(stackedPlotList::stackedPlotListVec_t::const_iterator it = spl.vec.begin(); it != spl.vec.end(); it++, i++)
    {
	const string curName = "h" + name + toString(i);
	const string drawstring = toDraw + ">>" + curName + "(" + toString(nBins) + "," + toString(xmin) + "," + toString(xmax) + ")";
	string curcut("");
	if (cut.size() > 0) curcut = cut;
	if (cut.size() > 0 && it->cut.size() > 0) curcut += "&&";
	if (it->cut.size() > 0) curcut += it->cut;

	cout << drawstring.c_str() << " - " << curcut << " scale: " << it->scale << endl;
	it->tree->Draw(drawstring.c_str(),curcut.c_str());
	TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(curName.c_str());
	h->SetFillColor(it->color);
	h->GetXaxis()->SetTitle(curXtitle.c_str());
	h->GetYaxis()->SetTitle(curYtitle.c_str());
	h->GetYaxis()->SetTitle(curYtitle.c_str());
	h->Scale(it->scale);
	legend->AddEntry(curName.c_str(),it->title.c_str(),"f");
	hs->Add(h);
    }
    hs->Draw();
    hs->GetXaxis()->SetTitle(curXtitle.c_str());
    hs->GetYaxis()->SetTitle(curYtitle.c_str());
    hs->GetYaxis()->SetTitleOffset(1.25);
    hs->SetTitle(title.c_str());
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.1);
    gPad->SetLogy(setLog);
    legend->Draw();
}

void dummy(TTree *t)
{
    stackedPlotList spl;
    spl.add("B^{0} #rightarrow J/#psi K_{s} (#pi#pi)","nRef2G==942",2); // our signal
    spl.add("B^{0} #rightarrow J/#psi K^{*0} (K^{0}#pi^{0})","nRef2G==946",3);
    spl.add("B^{0} #rightarrow J/#psi K^{*0} (K^{+}#pi^{-})","nRef2G==950",4);
    spl.add("B^{0} #rightarrow J/#psi K^{*0} (K^{+}#pi^{-}) #pi #pi","nRef2G==951",5);
    spl.add("B^{0} #rightarrow J/#psi K^{*+} (K^{0}#pi^{+})","nRef2G==969",6);
    spl.add("B^{0} #rightarrow J/#psi K_{1}(1279)^{0} (#rho^{0} K^{0})","nRef2G==889",7);
    spl.add("B^{0} #rightarrow J/#psi K_{1}(1279)^{0} (#rho(779)^{+} K^{+})","nRef2G==890",8);
    spl.add("B^{0} #rightarrow J/#psi K^{0} (K_{s})","nRef2G==939",9); 
    spl.add("B^{0} #rightarrow others","!(nRef2G==942||nRef2G==946||nRef2G==950||nRef2G==951||nRef2G==969||nRef2G==889||nRef2G==890||nRef2G==939)",11); 
    doStackedPlot(t, "name", "B^{0} decay modes in MC", "mB0", "m(#mu#mu#pi#pi)", "GeV/c^{2}", 40,4.6,6,spl);
}

