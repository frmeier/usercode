#include <vector>
#include <string>
#include <memory>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TPad.h"
#include "TLine.h"
#include "TF1.h"
#include "TLegend.h"
#include "TEventList.h"

#include "utils.h"
#include "Cuts.C"
#include "setTDRStyle_modified.C"

using std::string;
using std::auto_ptr;

void drawBitPlot(TTree *t, string hname, int N, string leafname, string cut)
{
    TH1F *h = new TH1F(hname.c_str(), hname.c_str(), N, 0, N);
    t->Draw(">>lst", cut.c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");

    int value;
    t->SetBranchAddress(leafname.c_str(), &value);

    for(Long64_t i = 0; i!=lst->GetN(); i++)
    {
	t->GetEntry(lst->GetEntry(i));
	for (int j=0; j!=N; j++)
	{
	    const int v = 0x1<<j;
	    if ((value & v) == v) h->Fill(j);
	}
    }
    h->Draw();

    delete lst;
}

// makes an TEfficiency plot. cutPass is just what you need in addition to cutAll
void doComparisonPlot(TCanvas *c, TTree *t, string hname, string title, string toDraw, int nBins, double lo, double hi,
	string atitle, string unit, const string &cut1, const string &cut2, const string &title1, const string &title2, bool setLog = false,  bool doScale = false)
{
    //cout << "cut2: " << cut2 << endl;
    //cout << "cut1: " << cut1 << endl;

    const string plotstring = makePlotsString(toDraw, hname, nBins, lo, hi);

    t->Draw(makePlotsString(toDraw, hname+"_1", nBins, lo, hi).c_str(), cut1.c_str());
    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject((hname+"_1").c_str());

    t->Draw(makePlotsString(toDraw, hname+"_2", nBins, lo, hi).c_str(), cut2.c_str());
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject((hname+"_2").c_str());
    if (doScale) h2->Scale(h1->GetEntries()/h2->GetEntries());

    double max = 1.1*(h1->GetMaximum() > h2->GetMaximum() ? h1->GetMaximum() : h2->GetMaximum());
    cout << max << endl;

    h1->SetMaximum(max);
    if (!setLog) h1->SetMinimum(0);
    h1->SetLineColor(2);
    h1->SetFillColor(2);
    h1->SetFillStyle(3004);
    h1->GetXaxis()->SetTitle(valueWithUnit(atitle, unit).c_str());
    h1->SetTitle(title.c_str());

    h2->SetMaximum(max);
    if (!setLog) h2->SetMinimum(0);
    h2->SetLineColor(4);
    h2->SetFillColor(4);
    h2->SetFillStyle(3005);
    h2->GetXaxis()->SetTitle(valueWithUnit(atitle, unit).c_str());
    h2->SetTitle(title.c_str());

    h1->Draw();
    h2->Draw("same");

    TLegend *leg = new TLegend(0.70,0.75, 0.90,0.85);
    leg->AddEntry(h1, title1.c_str(), "f");
    leg->AddEntry(h2, title2.c_str(), "f");
    leg->Draw();
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    gPad->SetTopMargin(0.10);
    gPad->SetLogy(setLog);

    c->SaveAs((hname+".pdf").c_str());
}

// makes an TEfficiency plot. cutPass is just what you need in addition to cutAll
void doComparisonPlotBit(TCanvas *c, TTree *t, string hname, string title, string toDraw, int nBins,
	string atitle, const string &cut1, const string &cut2, const string &title1, const string &title2, bool setLog = false, bool doScale = false)
{
    //cout << "cut2: " << cut2 << endl;
    //cout << "cut1: " << cut1 << endl;

    drawBitPlot(t, hname+"_1", nBins, toDraw, cut1);
    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject((hname+"_1").c_str());

    drawBitPlot(t, hname+"_2", nBins, toDraw, cut2);
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject((hname+"_2").c_str());
    if (doScale) h2->Scale(h1->GetEntries()/h2->GetEntries());

    double max = 1.1*(h1->GetMaximum() > h2->GetMaximum() ? h1->GetMaximum() : h2->GetMaximum());
    cout << max << endl;

    h1->SetMaximum(max);
    if (!setLog) h1->SetMinimum(0);
    h1->SetLineColor(2);
    h1->SetFillColor(2);
    h1->SetFillStyle(3004);
    h1->GetXaxis()->SetTitle(atitle.c_str());
    h1->SetTitle(title.c_str());

    h2->SetMaximum(max);
    if (!setLog) h2->SetMinimum(0);
    h2->SetLineColor(4);
    h2->SetFillColor(4);
    h2->SetFillStyle(3005);
    h2->GetXaxis()->SetTitle(atitle.c_str());
    h2->SetTitle(title.c_str());

    h1->Draw();
    h2->Draw("same");

    TLegend *leg = new TLegend(0.70,0.75, 0.90,0.85);
    leg->AddEntry(h1, title1.c_str(), "f");
    leg->AddEntry(h2, title2.c_str(), "f");
    leg->Draw();
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    gPad->SetTopMargin(0.10);
    gPad->SetLogy(setLog);

    c->SaveAs((hname+".pdf").c_str());
}

string getParticleLatex(string p)
{
    if (p == "lb") return "#Lambda_{b}";
    if (p == "Lb") return "#Lambda_{b}";
    if (p == "B0") return "B^{0}";
    if (p == "b0") return "B^{0}";
    return "dummy";
}

// =====================================================================================================
void doSomeComparisonPlots()
{
    setTDRStyle();
    const string p = "lb";
    TCanvas *c = new TCanvas("c", "c", 400, 400);
    //TFile *f = TFile::Open("../data/run456.root"); // MC
    //TFile *f = TFile::Open("../data/run458.root"); // Data
    //TFile *f = TFile::Open("../data/run477_478.root"); // Data
    //TFile *f = TFile::Open("../data/run480.root"); // Data
    TFile *f = TFile::Open("../data/run548.root"); // Data
    TTree *events = (TTree*)f->Get("events");

    Cuts cutLb;
    cutLb.selectCut("lb13exp", "acc06Lb", "muSoft", "HLT_jpsiBarrel", "HLT_matched");
    //cutLb.selectCut("lb12exp", "acc06Lb", "muSoft", "HLT_jpsiBarrel", "HLT_matched");
    //cutLb.selectCut("acc06Lb", "muSoft", "HLT_jpsiBarrel", "HLT_matched");
    //const string ana = cutLb.getCut();
    //const string ana = cutLb.getCut() + "&&npixha1>0&&npixha2>0&&d3rs/d3Ers>3";
    const string ana = cutLb.getCut();
    const string cutSig = "mbc>5.588&&mbc<5.650";
    const string cutSide = "((mbc>5.35&&mbc<5.56)||(mbc>5.67&&mbc<6))";
    const string anaSig = ana+"&&"+cutSig;
    const string anaSide = ana+"&&"+cutSide;

    const int nBins(40);

    /*
    doComparisonPlot(c, events, "mbc", "All after cuts", "mbc", nBins*2, 5.55, 5.70, "m(#Lambda_{b})", "GeV/c^{2}", ana+"&&rqha1>0", ana+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");

    doComparisonPlot(c, events, "ptrs_all", "All after cuts", "ptrs", nBins, 0, 40, "p_{T}(#Lambda^{0})", "GeV/c", ana+"&&rqha1>0", ana+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "ptrs_sig", "Signal region after cuts", "ptrs", nBins, 0, 40, "p_{T}(#Lambda^{0})", "GeV/c", anaSig+"&&rqha1>0", anaSig+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "ptrs_side", "Sideband region after cuts", "ptrs", nBins, 0, 40, "p_{T}(#Lambda^{0})", "GeV/c", anaSide+"&&rqha1>0", anaSide+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");

    doComparisonPlot(c, events, "mrs_all", "All after cuts", "mrs", nBins, 1.100, 1.130, "m(#Lambda^{0})", "GeV/c^{2}", ana+"&&rqha1>0", ana+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "mrs_sig", "Signal region after cuts", "mrs", nBins, 1.100, 1.130, "m(#Lambda^{0})", "GeV/c^{2}", anaSig+"&&rqha1>0", anaSig+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "mrs_side", "Sideband region after cuts", "mrs", nBins, 1.100, 1.130, "m(#Lambda^{0})", "GeV/c^{2}", anaSide+"&&rqha1>0", anaSide+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");

    doComparisonPlot(c, events, "probrs_all", "All after cuts", "probrs", .5*nBins, 0.0, 1.0, "prob(#Lambda^{0})", "GeV/c^{2}", ana+"&&rqha1>0", ana+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "probrs_sig", "Signal region after cuts", "probrs", .5*nBins, 0.0, 1.0, "prob(#Lambda^{0})", "GeV/c^{2}", anaSig+"&&rqha1>0", anaSig+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "probrs_side", "Sideband region after cuts", "probrs", .5*nBins, 0.0, 1.0, "prob(#Lambda^{0})", "GeV/c^{2}", anaSide+"&&rqha1>0", anaSide+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");

    doComparisonPlot(c, events, "probha1_all", "All after cuts", "probha1", .5*nBins, 0.0, 1.0, "prob(p)", "GeV/c^{2}", ana+"&&rqha1>0", ana+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "probha1_sig", "Signal region after cuts", "probha1", .5*nBins, 0.0, 1.0, "prob(p)", "GeV/c^{2}", anaSig+"&&rqha1>0", anaSig+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "probha1_side", "Sideband region after cuts", "probha1", .5*nBins, 0.0, 1.0, "prob(p)", "GeV/c^{2}", anaSide+"&&rqha1>0", anaSide+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");

    doComparisonPlot(c, events, "ndofha1_all", "All after cuts", "ndofha1", 50, 0.0, 50, "ndof(p)", "GeV/c^{2}", ana+"&&rqha1>0", ana+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "ndofha1_sig", "Signal region after cuts", "ndofha1", 50, 0.0, 50, "ndof(p)", "GeV/c^{2}", anaSig+"&&rqha1>0", anaSig+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "ndofha1_side", "Sideband region after cuts", "ndofha1", 50, 0.0, 50, "ndof(p)", "GeV/c^{2}", anaSide+"&&rqha1>0", anaSide+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    */

    /*
    doComparisonPlot(c, events, "ndofha1_1_all", "All after cuts", "ndofha1", 50, 0.0, 50, "ndof(p)", "GeV/c^{2}", ana+"&&rqha1>0&&ndofha1!=9", ana+"&&rqha1<0&&ndofha1!=9", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "ndofha1_1_sig", "Signal region after cuts", "ndofha1", 50, 0.0, 50, "ndof(p)", "GeV/c^{2}", anaSig+"&&rqha1>0&&ndofha1!=9", anaSig+"&&rqha1<0&&ndofha1!=9", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    doComparisonPlot(c, events, "ndofha1_1_side", "Sideband region after cuts", "ndofha1", 50, 0.0, 50, "ndof(p)", "GeV/c^{2}", anaSide+"&&rqha1>0&&ndofha1!=9", anaSide+"&&rqha1<0&&ndofha1!=9", "#Lambda_{b}", "#bar{#Lambda}_{b}");
    */

    /*
    doComparisonPlot(c, events, "ct3dbc_all", "All after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", ana+"&&rqha1>0", ana+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);
    doComparisonPlot(c, events, "ct3dbc_sig", "Signal region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", anaSig+"&&rqha1>0", anaSig+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);
    doComparisonPlot(c, events, "ct3dbc_side", "Sideband region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", anaSide+"&&rqha1>0", anaSide+"&&rqha1<0", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);
    */

    // Plots for investigating ndof=9 and bump region
    //const string sel = "&&ct3dbc>2.5e-12&&ct3dbc<4e-12";
    //const string nosel = "&&((ct3dbc<2.5e-12)||(ct3dbc>4e-12))";
    //const string textsel = "bump region";
    //const string textnosel = "outside";
    const string sel = "&&ndofha1==9";
    const string nosel = "&&ndofha1!=9";
    const string textsel = "ndof==9";
    const string textnosel = "ndof!=9";
    const string lb = "&&rqha1>0";
    const string lbbar = "&&rqha1<0";

    doComparisonPlot(c, events, "ct3dbc_sigbump_lb_lbbar", "Signal region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", anaSig+lb, anaSig+lbbar, "#Lambda_{b}", "#bar{#Lambda}_{b}", true, true);

    /*
    doComparisonPlot(c, events, "ct3dbc_sigbump_lb", "Signal region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, true, true);
    doComparisonPlot(c, events, "ct3dbc_sigbump_lbb", "Signal region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#bar{#Lambda}_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, true, true);

    doComparisonPlot(c, events, "etabc_sigbump_lb", "Signal region after cuts", "etabc", nBins, -2.4, 2.4, "#eta(#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "etabc_sigbump_lbb", "Signal region after cuts", "etabc", nBins, -2.4, 2.4, "#eta(#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phibc_sigbump_lb", "Signal region after cuts", "phibc", nBins, -3.14, 3.14, "#phi(#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phibc_sigbump_lbb", "Signal region after cuts", "phibc", nBins, -3.14, 3.14, "#phi(#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptbc_sigbump_lb", "Signal region after cuts", "ptbc", nBins, 0, 40, "p_{t}(#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptbc_sigbump_lbb", "Signal region after cuts", "ptbc", nBins, 0, 40, "p_{t}(#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);

    doComparisonPlot(c, events, "etajp_sigbump_lb", "Signal region after cuts", "etajp", nBins, -2.4, 2.4, "#eta(J/#psi) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "etajp_sigbump_lbb", "Signal region after cuts", "etajp", nBins, -2.4, 2.4, "#eta(J/#psi) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phijp_sigbump_lb", "Signal region after cuts", "phijp", nBins, -3.14, 3.14, "#phi(J/#psi) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phijp_sigbump_lbb", "Signal region after cuts", "phijp", nBins, -3.14, 3.14, "#phi(J/#psi) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptjp_sigbump_lb", "Signal region after cuts", "ptjp", nBins, 0, 40, "p_{t}(J/#psi) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptjp_sigbump_lbb", "Signal region after cuts", "ptjp", nBins, 0, 40, "p_{t}(J/#psi) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);

    doComparisonPlot(c, events, "etars_sigbump_lb", "Signal region after cuts", "etars", nBins, -2.4, 2.4, "#eta(#Lambda^{0}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "etars_sigbump_lbb", "Signal region after cuts", "etars", nBins, -2.4, 2.4, "#eta(#Lambda^{0}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phirs_sigbump_lb", "Signal region after cuts", "phirs", nBins, -3.14, 3.14, "#phi(#Lambda^{0}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phirs_sigbump_lbb", "Signal region after cuts", "phirs", nBins, -3.14, 3.14, "#phi(#Lambda^{0}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptrs_sigbump_lb", "Signal region after cuts", "ptrs", nBins, 0, 40, "p_{t}(#Lambda^{0}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptrs_sigbump_lbb", "Signal region after cuts", "ptrs", nBins, 0, 40, "p_{t}(#Lambda^{0}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);

    doComparisonPlot(c, events, "etamu1_sigbump_lb", "Signal region after cuts", "retamu1", nBins, -2.4, 2.4, "#eta(#mu_{1}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "etamu1_sigbump_lbb", "Signal region after cuts", "retamu1", nBins, -2.4, 2.4, "#eta(#mu_{1}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phimu1_sigbump_lb", "Signal region after cuts", "rphimu1", nBins, -3.14, 3.14, "#phi(#mu_{1}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phimu1_sigbump_lbb", "Signal region after cuts", "rphimu1", nBins, -3.14, 3.14, "#phi(#mu_{1}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptmu1_sigbump_lb", "Signal region after cuts", "rptmu1", nBins, 0, 40, "p_{t}(#mu_{1}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptmu1_sigbump_lbb", "Signal region after cuts", "rptmu1", nBins, 0, 40, "p_{t}(#mu_{1}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ndofmu1_sigbump_lb", "Signal region after cuts", "ndofmu1", 20, 0, 20, "ndof(#mu_{1}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ndofmu1_sigbump_lbb", "Signal region after cuts", "ndofmu1", 20, 0, 20, "ndof(#mu_{1}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "npixmu1_sigbump_lb", "Signal region after cuts", "npixmu1", 10, 0, 10, "npix(#mu_{1}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "npixmu1_sigbump_lbb", "Signal region after cuts", "npixmu1", 10, 0, 10, "npix(#mu_{1}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlotBit(c, events, "qualmu1_sigbump_lb", "Signal region after cuts", "qualmu1", 10, "qual(#mu_{1}) (#Lambda_{b})", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlotBit(c, events, "qualmu1_sigbump_lbb", "Signal region after cuts", "qualmu1", 10, "qual(#mu_{1}) (#bar#Lambda_{b})", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);

    doComparisonPlot(c, events, "etamu2_sigbump_lb", "Signal region after cuts", "retamu2", nBins, -2.4, 2.4, "#eta(#mu_{2}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "etamu2_sigbump_lbb", "Signal region after cuts", "retamu2", nBins, -2.4, 2.4, "#eta(#mu_{2}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phimu2_sigbump_lb", "Signal region after cuts", "rphimu2", nBins, -3.14, 3.14, "#phi(#mu_{2}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phimu2_sigbump_lbb", "Signal region after cuts", "rphimu2", nBins, -3.14, 3.14, "#phi(#mu_{2}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptmu2_sigbump_lb", "Signal region after cuts", "rptmu2", nBins, 0, 40, "p_{t}(#mu_{2}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptmu2_sigbump_lbb", "Signal region after cuts", "rptmu2", nBins, 0, 40, "p_{t}(#mu_{2}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ndofmu2_sigbump_lb", "Signal region after cuts", "ndofmu2", 20, 0, 20, "ndof(#mu_{2}) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ndofmu2_sigbump_lbb", "Signal region after cuts", "ndofmu2", 20, 0, 20, "ndof(#mu_{2}) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlotBit(c, events, "qualmu2_sigbump_lb", "Signal region after cuts", "qualmu2", 10, "qual(#mu_{2}) (#Lambda_{b})", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlotBit(c, events, "qualmu2_sigbump_lbb", "Signal region after cuts", "qualmu2", 10, "qual(#mu_{2}) (#bar#Lambda_{b})", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);

    doComparisonPlot(c, events, "etaha1_sigbump_lb", "Signal region after cuts", "retaha1", nBins, -2.4, 2.4, "#eta(p) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "etaha1_sigbump_lbb", "Signal region after cuts", "retaha1", nBins, -2.4, 2.4, "#eta(p) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phiha1_sigbump_lb", "Signal region after cuts", "rphiha1", nBins, -3.14, 3.14, "#phi(p) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phiha1_sigbump_lbb", "Signal region after cuts", "rphiha1", nBins, -3.14, 3.14, "#phi(p) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptha1_sigbump_lb", "Signal region after cuts", "rptha1", nBins, 0, 15, "p_{t}(p) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptha1_sigbump_lbb", "Signal region after cuts", "rptha1", nBins, 0, 15, "p_{t}(p) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ndofha1_sigbump_lb", "Signal region after cuts", "ndofha1", 20, 0, 20, "ndof(p) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ndofha1_sigbump_lbb", "Signal region after cuts", "ndofha1", 20, 0, 20, "ndof(p) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "npixha1_sigbump_lb", "Signal region after cuts", "npixha1", 10, 0, 10, "npix(p) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "npixha1_sigbump_lbb", "Signal region after cuts", "npixha1", 10, 0, 10, "npix(p) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlotBit(c, events, "qualha1_sigbump_lb", "Signal region after cuts", "qualha1", 10, "qual(p) (#Lambda_{b})", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlotBit(c, events, "qualha1_sigbump_lbb", "Signal region after cuts", "qualha1", 10, "qual(p) (#bar#Lambda_{b})", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);

    doComparisonPlot(c, events, "etaha2_sigbump_lb", "Signal region after cuts", "retaha2", nBins, -2.4, 2.4, "#eta(#pi) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "etaha2_sigbump_lbb", "Signal region after cuts", "retaha2", nBins, -2.4, 2.4, "#eta(#pi) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phiha2_sigbump_lb", "Signal region after cuts", "rphiha2", nBins, -3.14, 3.14, "#phi(#pi) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "phiha2_sigbump_lbb", "Signal region after cuts", "rphiha2", nBins, -3.14, 3.14, "#phi(#pi) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptha2_sigbump_lb", "Signal region after cuts", "rptha2", nBins, 0, 5, "p_{t}(#pi) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ptha2_sigbump_lbb", "Signal region after cuts", "rptha2", nBins, 0, 5, "p_{t}(#pi) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ndofha2_sigbump_lb", "Signal region after cuts", "ndofha2", 20, 0, 20, "ndof(#pi) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ndofha2_sigbump_lbb", "Signal region after cuts", "ndofha2", 20, 0, 20, "ndof(#pi) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "npixha2_sigbump_lb", "Signal region after cuts", "npixha2", 10, 0, 10, "npix(#pi) (#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "npixha2_sigbump_lbb", "Signal region after cuts", "npixha2", 10, 0, 10, "npix(#pi) (#bar#Lambda_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlotBit(c, events, "qualha2_sigbump_lb", "Signal region after cuts", "qualha2", 10, "qual(#pi) (#Lambda_{b})", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlotBit(c, events, "qualha2_sigbump_lbb", "Signal region after cuts", "qualha2", 10, "qual(#pi) (#bar#Lambda_{b})", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);

    doComparisonPlot(c, events, "d3Ebc_sigbump_lb", "Signal region after cuts", "d3Ebc", nBins, 0, .04, "d3E(#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "d3Ebc_sigbump_lbb", "Signal region after cuts", "d3Ebc", nBins, 0, .04, "d3E(#bar(#Lambda)_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ct3dbcE_sigbump_lb", "Signal region after cuts", "ct3dbcE", nBins, 0, .2e-12, "ct3E(#Lambda_{b})", "s", anaSig+sel+lb, anaSig+nosel+lb, textsel, textnosel, false, true);
    doComparisonPlot(c, events, "ct3dbcE_sigbump_lbb", "Signal region after cuts", "ct3dbcE", nBins, 0, .2e-12, "ct3E(#bar(#Lambda)_{b})", "s", anaSig+sel+lbbar, anaSig+nosel+lbbar, textsel, textnosel, false, true);
    */

    /*
    string anaSigLb = anaSig+"&&rqha1>0";
    string anaSigLbbar = anaSig+"&&rqha1<0";
    doComparisonPlot(c, events, "phiha1_Lb", "Signal region after cuts", "rphiha1", nBins, -3.141, 3.141, "#phi(p)", "rad", anaSigLb+"&&ndofha1!=9", anaSigLb+"&&ndofha1==9", "#Lambda_{b}, ndof!=9", "#Lambda_{b}, ndof==9", false);
    doComparisonPlot(c, events, "phiha1_Lbbar", "Signal region after cuts", "rphiha1", nBins, -3.141, 3.141, "#phi(p)", "rad", anaSigLbbar+"&&ndofha1!=9", anaSigLbbar+"&&ndofha1==9", "#bar{#Lambda}_{b}, ndof!=9", "#bar{#Lambda}_{b}, ndof==9", false);
    doComparisonPlot(c, events, "phiha1_sig", "Signal region after cuts", "rphiha1", nBins, -3.141, 3.141, "#phi(p)", "rad", anaSigLb, anaSigLbbar, "#Lambda_{b}", "#bar{#Lambda}_{b}", false);

    doComparisonPlot(c, events, "etaha1_Lb", "Signal region after cuts", "retaha1", nBins, -2.5, 2.5, "#eta(p)", "", anaSigLb+"&&ndofha1!=9", anaSigLb+"&&ndofha1==9", "#Lambda_{b}, ndof!=9", "#Lambda_{b}, ndof==9", false);
    doComparisonPlot(c, events, "etaha1_Lbbar", "Signal region after cuts", "retaha1", nBins, -2.5, 2.5, "#eta(p)", "", anaSigLbbar+"&&ndofha1!=9", anaSigLbbar+"&&ndofha1==9", "#bar{#Lambda}_{b}, ndof!=9", "#bar{#Lambda}_{b}, ndof==9", false);
    doComparisonPlot(c, events, "etaha1_sig", "Signal region after cuts", "retaha1", nBins, -2.5, 2.5, "#eta(p)", "", anaSigLb, anaSigLbbar, "#Lambda_{b}", "#bar{#Lambda}_{b}", false);

    doComparisonPlot(c, events, "ptha1_Lb", "Signal region after cuts", "rptha1", nBins, 0, 10, "#pt(p)", "GeV/c", anaSigLb+"&&ndofha1!=9", anaSigLb+"&&ndofha1==9", "#Lambda_{b}, ndof!=9", "#Lambda_{b}, ndof==9", false);
    doComparisonPlot(c, events, "ptha1_Lbbar", "Signal region after cuts", "rptha1", nBins, 0, 10, "#pt(p)", "GeV/c", anaSigLbbar+"&&ndofha1!=9", anaSigLbbar+"&&ndofha1==9", "#bar{#Lambda}_{b}, ndof!=9", "#bar{#Lambda}_{b}, ndof==9", false);
    doComparisonPlot(c, events, "ptha1_sig", "Signal region after cuts", "rptha1", nBins, 0, 10, "#pt(p)", "GeV/c", anaSigLb, anaSigLbbar, "#Lambda_{b}", "#bar{#Lambda}_{b}", false);

    doComparisonPlot(c, events, "ndofha2_Lb", "Signal region after cuts", "ndofha2", 50, 0,50, "ndof(#pi)", "", anaSigLb+"&&ndofha1!=9", anaSigLb+"&&ndofha1==9", "#Lambda_{b}, ndof!=9", "#Lambda_{b}, ndof==9", false);
    doComparisonPlot(c, events, "ndofha2_Lbbar", "Signal region after cuts", "ndofha2", 50, 0,50, "ndof(#pi)", "", anaSigLbbar+"&&ndofha1!=9", anaSigLbbar+"&&ndofha1==9", "#bar{#Lambda}_{b}, ndof!=9", "#bar{#Lambda}_{b}, ndof==9", false);
    doComparisonPlot(c, events, "ndofha2_sig", "Signal region after cuts", "ndofha2", 50, 0,50, "ndof(#pi)", "", anaSigLb, anaSigLbbar, "#Lambda_{b}", "#bar{#Lambda}_{b}", false);
    */

    /*
    doComparisonPlot(c, events, "ct3dbc_2_all", "All after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", ana+"&&rqha1>0&&ndofha1!=9", ana+"&&rqha1<0&&ndofha1!=9", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);
    doComparisonPlot(c, events, "ct3dbc_2_sig", "Signal region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", anaSig+"&&rqha1>0&&ndofha1!=9", anaSig+"&&rqha1<0&&ndofha1!=9", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);
    doComparisonPlot(c, events, "ct3dbc_2_side", "Sideband region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", anaSide+"&&rqha1>0&&ndofha1!=9", anaSide+"&&rqha1<0&&ndofha1!=9", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);

    doComparisonPlot(c, events, "ct3dbc_3_all", "All after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", ana+"&&rqha1>0&&(probha1<.25||probha1>.45)", ana+"&&rqha1<0&&(probha1<.25||probha1>.45)", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);
    doComparisonPlot(c, events, "ct3dbc_3_sig", "Signal region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", anaSig+"&&rqha1>0&&(probha1<.25||probha1>.45)", anaSig+"&&rqha1<0&&(probha1<.25||probha1>.45)", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);
    doComparisonPlot(c, events, "ct3dbc_3_side", "Sideband region after cuts", "ct3dbc", nBins, -5e-12, 15e-12, "t(#Lambda_{b})", "s", anaSide+"&&rqha1>0&&(probha1<.25||probha1>.45)", anaSide+"&&rqha1<0&&(probha1<.25||probha1>.45)", "#Lambda_{b}", "#bar{#Lambda}_{b}", true);
    */
}


