#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <stdexcept>
#include "setTDRStyle_modified.C"

TCanvas *c;
TFile *f;
const std::string suffix("_16");
//const std::string suffix("");
const std::string prefix("rooFitCutScan");
//const std::string run1("run129p"), title1("MC #Lambda_{b}");
const std::string run1("run129p16"), title1("MC #Lambda_{b}");
const std::string run2("run127"), title2("data");

enum legposenum { posUL, posUR, posLL, posLR };

TH1F* extractHisto(std::string filename, std::string histoname, std::string newname)
{
    cout << "Request to open " << filename << endl;
    TFile *curf = TFile::Open(filename.c_str());
    if (curf == 0) throw std::runtime_error("File " + filename + " not found");
    TH1F* h = (TH1F*)curf->Get(histoname.c_str());
    //cout << h << endl;
    if (h == 0) throw std::runtime_error("Histo " + histoname + " not present in file");
    f->cd();
    TH1F* hret = (TH1F*)h->Clone(newname.c_str());
    delete h;
    delete curf;
    //cout << "ok" << endl;
    return hret;
}

void rooFitCutScanPlot(std::string cutname, std::string cutTexname = "", legposenum legpos = posLR)
{
    setTDRStyle();

    if(0==c) c = new TCanvas("c","c",400,400);
    f = new TFile("tmp.root","recreate");

    TH1F *h1 = extractHisto(prefix+"_"+cutname+"_"+run1+suffix+".root", "hsig", "h1");
    h1->SetMarkerColor(4);
    h1->SetLineColor(4);
    h1->SetFillColor(33);
    h1->GetYaxis()->SetTitle(h1->GetTitle());
    h1->GetXaxis()->SetTitle(("Variation of cut on " + cutTexname.size()==0 ? cutname : cutTexname).c_str());
    h1->SetTitle("");
    h1->Draw("P0E3");
    h1->Draw("P0E2same");
    TH1F *h2 = extractHisto(prefix+"_"+cutname+"_"+run2+suffix+".root", "hsig", "h1");
    h2->SetMarkerColor(1);
    h2->SetLineColor(1);
    h2->SetFillColor(1);
    h2->SetTitle(cutname.c_str());
    h2->Draw("P0E2same");

    TLegend *leg;
    const double legwidth = 0.17;
    const double legheight = 0.14;
    if (legpos == posUL) { double x(0.22), y(0.70); leg = new TLegend(x,y,x+legwidth,y+legheight); }
    if (legpos == posUR) { double x(0.72), y(0.70); leg = new TLegend(x,y,x+legwidth,y+legheight); }
    if (legpos == posLL) { double x(0.22), y(0.20); leg = new TLegend(x,y,x+legwidth,y+legheight); }
    if (legpos == posLR) { double x(0.72), y(0.20); leg = new TLegend(x,y,x+legwidth,y+legheight); }
    leg->AddEntry(h1,title1.c_str(),"pf");
    leg->AddEntry(h2,title2.c_str(),"lp");
    leg->SetBorderSize(0);
    leg->Draw();
}

struct Cutplot
{
    Cutplot(std::string cn, std::string ctn, legposenum pos = posLR) : cutname(cn), cutTexname(ctn), legpos(pos) {};
    std::string cutname, cutTexname;
    legposenum legpos;
};

void rooFitCutScanPlotAll()
{
    std::vector<Cutplot> cutVec;
    //cutnameVec.push_back("mlb"); cutTexnameVec.push_back("m(#Lambda)");
    // list of plots to make
    cutVec.push_back(Cutplot("ml0","m(#Lambda) [GeV/c^{2}]"));
    cutVec.push_back(Cutplot("mjp","m(J/#psi) [GeV/c^{2}]"));
    cutVec.push_back(Cutplot("alphal0","#alpha(#Lambda) [rad]"));
    cutVec.push_back(Cutplot("alphalb","#alpha(#Lambda_{b}) [rad]"));
    cutVec.push_back(Cutplot("d3l0d3El0","d3(#Lambda)/d3_{err}(#Lambda)"));
    cutVec.push_back(Cutplot("d3l0","d3(#Lambda) [cm]", posLL));
    cutVec.push_back(Cutplot("Kshypo","veto m(K_{s}) [GeV/c^{2}]", posLL));
    cutVec.push_back(Cutplot("problb","Prob(#chi^{2}(#Lambda_{b}))"));
    cutVec.push_back(Cutplot("ptjp","p_{T}(J/#psi) [GeV/c]", posLL));
    cutVec.push_back(Cutplot("ptl0","p_{T}(#Lambda) [GeV/c]", posUR));
    cutVec.push_back(Cutplot("rptpi","p_{T}(#pi) [GeV/c]", posUR));
    cutVec.push_back(Cutplot("rptpr","p_{T}(p) [GeV/c]", posUR));

    // now do the plots
    for(std::vector<Cutplot>::const_iterator it = cutVec.begin(); it!=cutVec.end(); it++)
    {
	rooFitCutScanPlot(it->cutname,it->cutTexname,it->legpos);
	c->SaveAs(("rooFitCutScanPlot_" + it->cutname + suffix + ".pdf").c_str());
    }
}


