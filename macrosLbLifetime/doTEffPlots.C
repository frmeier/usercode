#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <list>

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TEventList.h"
#include "TPad.h"
#include "TLine.h"
#include "TF1.h"

#include "utils.h"
#include "Cuts.C"

using std::string;
using std::auto_ptr;

std::vector<double> makeDynamicBins(TTree *t, string branchname, string cut, int nBins, double lo, double hi, const double factor = 1)
{
    typedef std::list<double> listtype;
    listtype valuelist;

    t->Draw(">>lst", cut.c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");

    double value;
    t->SetBranchAddress(branchname.c_str(), &value);
    const int N = lst->GetN();

    for (int i=0; i!=N; i++)
    {
	t->GetEntry(lst->GetEntry(i));
	value*=factor;
	if (value>=lo && value <=hi)
	{
	    valuelist.push_back(value);
	}
    }
    valuelist.sort();
    const int entriesPerBin = valuelist.size() / nBins;

    std::vector<double> binvector;
    binvector.push_back(lo);

    int i = 0;
    for (listtype::const_iterator it = valuelist.begin(); it!=valuelist.end(); it++, i++)
    {
	if (i==0) continue;
	if (i % entriesPerBin == 0) binvector.push_back(*it);
    }
    binvector.pop_back(); // remove last one (ugly, but easier than awful logic in loop)
    binvector.push_back(hi); // and replace it with given boundary
    for (std::vector<double>::const_iterator it = binvector.begin(); it!=binvector.end(); it++)
	cout << *it << " ";
    cout << endl;

    return binvector;
}


// makes an TEfficiency plot. cutPass is just what you need in addition to cutAll
void doTEffPlot(TCanvas *c, TTree *t, string hname, string title, string toDraw, int nBins, double lo, double hi,
	string unit, const string &cutAll, const string &cutPass)
{
    const string plotstring = makePlotsString(toDraw, hname, nBins, lo, hi);

    t->Draw(makePlotsString(toDraw, hname+"_all", nBins, lo, hi).c_str(), cutAll.c_str());
    TH1F *hall = (TH1F*)gDirectory->GetList()->FindObject((hname+"_all").c_str());

    t->Draw(makePlotsString(toDraw, hname+"_pass", nBins, lo, hi).c_str(), (cutAll.size() > 0 ? cutAll+"&&"+cutPass : cutPass).c_str());
    TH1F *hpass = (TH1F*)gDirectory->GetList()->FindObject((hname+"_pass").c_str());

    auto_ptr<TEfficiency> heff (new TEfficiency(*hpass, *hall));
    heff->SetTitle((title+";"+unit+";#varepsilon").c_str());
    heff->Draw();

    const double averageEff = (double)hpass->GetEntries()/(double)hall->GetEntries();
    double x1, x2, y1, y2;
    gPad->Update();
    gPad->GetRangeAxis(x1, y1, x2, y2);
    TLine *tl = new TLine(x1, averageEff, x2, averageEff);
    tl->SetLineColor(3);
    tl->Draw();
    gPad->Update();
    c->SaveAs((hname+".pdf").c_str());

    cout << hname << ": " << hpass->GetEntries() << "/" << hall->GetEntries() << ": " << averageEff << endl;
}

void doTEffPlotFit(TCanvas *c, TTree *t, string hname, string title, string toDraw, int nBins, double lo, double hi,
	string unit, const string &cutAll, const string &cutPass)
{
    const string plotstring = makePlotsString(toDraw, hname, nBins, lo, hi);

    t->Draw(makePlotsString(toDraw, hname+"_all", nBins, lo, hi).c_str(), cutAll.c_str());
    TH1F *hall = (TH1F*)gDirectory->GetList()->FindObject((hname+"_all").c_str());

    t->Draw(makePlotsString(toDraw, hname+"_pass", nBins, lo, hi).c_str(), (cutAll.size() > 0 ? cutAll+"&&"+cutPass : cutPass).c_str());
    TH1F *hpass = (TH1F*)gDirectory->GetList()->FindObject((hname+"_pass").c_str());

    auto_ptr<TEfficiency> heff (new TEfficiency(*hpass, *hall));
    heff->SetTitle((title+";"+unit+";#varepsilon").c_str());
    heff->Draw();

    const double averageEff = (double)hpass->GetEntries()/(double)hall->GetEntries();
    double x1, x2, y1, y2;
    gPad->Update();
    gPad->GetRangeAxis(x1, y1, x2, y2);
    TLine *tl = new TLine(x1, averageEff, x2, averageEff);
    tl->SetLineColor(3);
    tl->Draw();
    gPad->Update();

    // now we fit a pol1
    TF1* f1 = new TF1("f1","[0]*(1+[1]*x)",lo>0 ? lo : 0,hi);
    f1->SetParameters(averageEff, -0.00);
    heff->Fit(f1, "I");

    c->SaveAs((hname+".pdf").c_str());

    cout << hname << ": " << hpass->GetEntries() << "/" << hall->GetEntries() << ": " << averageEff << endl;
}

void doTEffPlotFit(TCanvas *c, TTree *t, string hname, string title, string toDraw, const std::vector<double> &bins,
	string unit, const string &cutAll, const string &cutPass, double factor = 1)
{
    cout << bins.size()-1 << endl;
    TH1F *hall = new TH1F((hname+"_all").c_str(), (hname+"_all").c_str(), bins.size()-1, &bins[0]);
    t->Draw(makePlotsString(toDraw+(factor!=1?"*"+toString(factor):""), hname+"_all").c_str(), cutAll.c_str());
    c->SaveAs("test_all.pdf");

    TH1F *hpass = new TH1F((hname+"_pass").c_str(), (hname+"_pass").c_str(), bins.size()-1, &bins[0]);
    t->Draw(makePlotsString(toDraw+(factor!=1?"*"+toString(factor):""), hname+"_pass").c_str(), (cutAll.size() > 0 ? cutAll+"&&"+cutPass : cutPass).c_str());
    c->SaveAs("test_pass.pdf");

    auto_ptr<TEfficiency> heff (new TEfficiency(*hpass, *hall));
    heff->SetTitle((title+";"+unit+";#varepsilon").c_str());
    heff->Draw();
    gPad->Update();

    heff->GetPaintedGraph()->GetXaxis()->SetLabelSize(0.05);
    heff->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.05);

    heff->GetPaintedGraph()->GetYaxis()->SetLabelSize(0.05);
    heff->GetPaintedGraph()->GetYaxis()->SetTitleOffset(.8);
    heff->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.05);

    const double averageEff = (double)hpass->GetEntries()/(double)hall->GetEntries();
    double x1, x2, y1, y2;
    gPad->Update();
    gPad->GetRangeAxis(x1, y1, x2, y2);
    TLine *tl = new TLine(x1, averageEff, x2, averageEff);
    tl->SetLineColor(3);
    tl->Draw();
    gPad->SetLeftMargin(.08);
    gPad->SetRightMargin(.02);
    gPad->SetBottomMargin(.12);
    gPad->Update();

    // now we fit a pol1
    const double lo = bins[0];
    const double hi = bins[bins.size()-1];
    TF1* f1 = new TF1("f1","[0]*(1+[1]*x)",lo>0 ? lo : 0,hi);
    f1->SetParameters(averageEff, -0.00);
    heff->Fit(f1, "I");

    c->SaveAs((hname+".pdf").c_str());

    cout << hname << ": " << hpass->GetEntries() << "/" << hall->GetEntries() << ": " << averageEff << endl;
}

void doTEffPlotFit(TCanvas *c, string hnamepass, string hnameall, string hname, string title, string unit)
{
    TH1F *hall = (TH1F*)gDirectory->GetList()->FindObject(hnameall.c_str());
    TH1F *hpass = (TH1F*)gDirectory->GetList()->FindObject(hnamepass.c_str());

    auto_ptr<TEfficiency> heff (new TEfficiency(*hpass, *hall));
    heff->SetTitle((title+";"+unit+";#varepsilon").c_str());
    heff->Draw();

    gPad->Update();

    heff->GetPaintedGraph()->GetXaxis()->SetLabelSize(0.05);
    heff->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.05);

    heff->GetPaintedGraph()->GetYaxis()->SetLabelSize(0.05);
    heff->GetPaintedGraph()->GetYaxis()->SetTitleOffset(.8);
    heff->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.05);

    const double averageEff = (double)hpass->GetEntries()/(double)hall->GetEntries();
    double x1, x2, y1, y2;
    gPad->Update();
    gPad->GetRangeAxis(x1, y1, x2, y2);
    TLine *tl = new TLine(x1, averageEff, x2, averageEff);
    tl->SetLineColor(3);
    tl->Draw();
    gPad->SetLeftMargin(.08);
    gPad->SetRightMargin(.02);
    gPad->SetBottomMargin(.12);
    gPad->Update();

    const double lo = hpass->GetBinLowEdge(1);
    const double hi = hpass->GetBinLowEdge(hpass->GetNbinsX()+1);

    // now we fit a pol1
    TF1* f1 = new TF1("f1","[0]*(1+[1]*x)",lo,hi);
    //f1->SetParameters(averageEff, -0.03e12);
    f1->SetParameters(averageEff, 0);
    heff->Fit(f1, "I");

    c->SaveAs((hname+".pdf").c_str());

    cout << hname << ": " << hpass->GetEntries() << "/" << hall->GetEntries() << ": " << averageEff << endl;
}

// makes an TEfficiency plot. cutPass is just what you need in addition to cutAll
void doTEffPlot(TCanvas *c, TTree *t, string hname, string title, string toDraw, const std::vector<double> &bins,
	string unit, const string &cutAll, const string &cutPass)
{
    TH1F *hall = new TH1F((hname+"_all").c_str(), title.c_str(), bins.size()-1, &bins[0]);
    t->Draw((toDraw+">>"+hname+"_all").c_str(), cutAll.c_str());

    TH1F *hpass = new TH1F((hname+"_pass").c_str(), title.c_str(), bins.size()-1, &bins[0]);
    t->Draw((toDraw+">>"+hname+"_pass").c_str(), (cutAll.size() > 0 ? cutAll+"&&"+cutPass : cutPass).c_str());

    auto_ptr<TEfficiency> heff (new TEfficiency(*hpass, *hall));
    heff->SetTitle((title+";"+unit+";#varepsilon").c_str());
    heff->Draw();

    const double averageEff = (double)hpass->GetEntries()/(double)hall->GetEntries();
    double x1, x2, y1, y2;
    gPad->Update();
    gPad->GetRangeAxis(x1, y1, x2, y2);
    TLine *tl = new TLine(x1, averageEff, x2, averageEff);
    tl->SetLineColor(3);
    tl->Draw();
    gPad->Update();
    c->SaveAs((hname+".pdf").c_str());

    cout << hname << ": " << hpass->GetEntries() << "/" << hall->GetEntries() << ": " << averageEff << endl;
}

// helper function for getting particles in root-LaTex style
string getParticleLatex(string p)
{
    if (p == "lb") return "#Lambda_{b}";
    if (p == "Lb") return "#Lambda_{b}";
    if (p == "B0") return "B^{0}";
    if (p == "b0") return "B^{0}";
    if (p == "L0") return "#Lambda^{0}";
    if (p == "l0") return "#Lambda^{0}";
    if (p == "Ks") return "K_{s}";
    return "dummy";
}

// =====================================================================================================
void doTEffPlotSetGen(TCanvas *c, TTree *t, string hname, string title, string particle, int nBins, const string &cutAll, const string &cutPass)
{
    //doTEffPlot(c, t, particle+"_"+hname+"_p", title, "pbc", nBins, 0, 100, valueWithUnit("p("+getParticleLatex(particle)+")","GeV/c"), cutAll, cutPass);
    //doTEffPlot(c, t, particle+"_"+hname+"_d3", title, "d3dbc", nBins, 0, 1, valueWithUnit("d_{3}("+getParticleLatex(particle)+")","cm"), cutAll, cutPass);
    doTEffPlotFit(c, t, particle+"_"+hname+"_ct", title, "ctbc*1e12", nBins, 0, 10, valueWithUnit("t("+getParticleLatex(particle)+")","ps"), cutAll, cutPass);
    doTEffPlotFit(c, t, particle+"_"+hname+"_ctz", title+" (zoom)", "ctbc*1e12", .4*nBins, 0, 2, valueWithUnit("t("+getParticleLatex(particle)+")","ps"), cutAll, cutPass);
}

void doTEffPlotSetGen_t(TCanvas *c, TTree *t, string hname, string title, string particle, const string &cutAll, const string &cutPass, std::vector<double> bins)
{
    doTEffPlotFit(c, t, particle+"_"+hname+"_ct", title, "ctbc", bins, valueWithUnit("t("+getParticleLatex(particle)+")","ps"), cutAll, cutPass, 1e12);
}

void doTEffPlotSetReco(TCanvas *c, TTree *t, string hname, string title, string particle, int nBins, const string &cutAll, const string &cutPass)
{
    //doTEffPlot(c, t, particle+"_"+hname+"_p", title, "pbctruth", nBins, 0, 100, valueWithUnit("p("+getParticleLatex(particle)+")","GeV/c"), cutAll, cutPass);
    //doTEffPlot(c, t, particle+"_"+hname+"_d3", title, "d3dbctruth", nBins, 0, 1, valueWithUnit("d_{3}("+getParticleLatex(particle)+")","cm"), cutAll, cutPass);
    doTEffPlotFit(c, t, particle+"_"+hname+"_ct", title, "ctbctruth*1e12", nBins, 0, 10, valueWithUnit("t("+getParticleLatex(particle)+")","ps"), cutAll, cutPass);
    doTEffPlotFit(c, t, particle+"_"+hname+"_ctz", title+" (zoom)", "ctbctruth*1e12", .4*nBins, 0, 2, valueWithUnit("t("+getParticleLatex(particle)+")","ps"), cutAll, cutPass);
}

void doTEffPlotSetReco_t(TCanvas *c, TTree *t, string hname, string title, string particle, const string &cutAll, const string &cutPass, std::vector<double> bins)
{
    doTEffPlotFit(c, t, particle+"_"+hname+"_ct", title, "ctbctruth", bins, valueWithUnit("t("+getParticleLatex(particle)+")","ps"), cutAll, cutPass,  1e12);
}

void doTEffPlotSetRecoRecoval(TCanvas *c, TTree *t, string hname, string title, string particle, int nBins, const string &cutAll, const string &cutPass)
{
    //doTEffPlot(c, t, particle+"_"+hname+"_reco_p", title, "pbc", nBins, 0, 100, valueWithUnit("p("+getParticleLatex(particle)+") reco","GeV/c"), cutAll, cutPass);
    //doTEffPlot(c, t, particle+"_"+hname+"_reco_d3", title, "d3bc", nBins, -0.2, 1, valueWithUnit("d_{3}("+getParticleLatex(particle)+") reco","cm"), cutAll, cutPass);
    doTEffPlotFit(c, t, particle+"_"+hname+"_reco_ct", title, "ct3dbc*1e12", nBins, -2, 10, valueWithUnit("t("+getParticleLatex(particle)+") reco","ps"), cutAll, cutPass);
    doTEffPlotFit(c, t, particle+"_"+hname+"_reco_ctz", title+" (zoom)", "ct3dbc*1e12", nBins, -.5, 2, valueWithUnit("t("+getParticleLatex(particle)+") reco","ps"), cutAll, cutPass);
}

void doTEffPlotR(TCanvas *c, TTree *t, string hname, string title, string particle, int nBins, const string &cutAll, const string &cutPass)
{
    string rsname = (particle == "lb") ? "L0" : "Ks";
    doTEffPlot(c, t, particle+"_"+hname+"_rrs", title, "vrrsPV", nBins, 0, 100, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_rrs_zoom", title, "vrrsPV", nBins, 0, 30, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_rrs_zoom2", title, "vrrsPV", nBins, 0, 5, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAll, cutPass);
    if (particle == "lb")
    {
	const string cutAllL0  = ((cutAll.size()!=0)  ? (cutAll + "&&")  : "");
	const string cutPassL0 = ((cutPass.size()!=0) ? (cutPass + "&&") : "");
	doTEffPlot(c, t, particle+"_"+hname+"_l0_rrs", title + " #Lambda^{0}", "vrrs", nBins, 0, 100, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAllL0+"qha1>0", cutPassL0+"qha1>0");
	doTEffPlot(c, t, particle+"_"+hname+"_l0_rrs_zoom", title + " #Lambda^{0}", "vrrs", nBins, 0, 30, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAllL0+"qha1>0", cutPassL0+"qha1>0");
	doTEffPlot(c, t, particle+"_"+hname+"_l0bar_rrs", title + " #bar{#Lambda}^{0}", "vrrs", nBins, 0, 100, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAllL0+"qha1<0", cutPassL0+"qha1<0");
	doTEffPlot(c, t, particle+"_"+hname+"_l0bar_rrs_zoom", title + " #bar{#Lambda}^{0}", "vrrs", nBins, 0, 30, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAllL0+"qha1<0", cutPassL0+"qha1<0");
    }

    /*
    std::vector<double> bins, bins2;
    bins.push_back(0.0);
    bins2.push_back(0.0);
    for (int i=0; i<=25; i++) bins.push_back(pow(10,i*.03-.3)+.5);
    for (int i=0; i<=40; i++) bins2.push_back(pow(10,i*.03-.3)+.5);
    doTEffPlot(c, t, particle+"_"+hname+"_rrs", title, "vrrs", bins, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_rrs2", title, "vrrs", bins2, valueWithUnit("r("+getParticleLatex(rsname)+")","cm"), cutAll, cutPass);
    */
}

// =====================================================================================================
void doTEffPlotL0()
{
    const string p = "lb";
    //std::vector<double> bins = variableBinSizeVec(0, 12, 14, 16, 18, 20, 24, 28, 36, 60, 100);

    TCanvas *c = new TCanvas("c", "c", 400, 400);
    TFile *f = TFile::Open("../data/run533.root");
    TTree *genevents = (TTree*)f->Get("genevents");
    //TTree *events = (TTree*)f->Get("events");

    doTEffPlotR(c, genevents, "rrs", "Reco efficiency", "lb", 100, "", "hasCand==1");
}

// =====================================================================================================
void doTEffPlotKs()
{
    const string p = "B0";
    //std::vector<double> bins = variableBinSizeVec(0, 12, 14, 16, 18, 20, 24, 28, 36, 60, 100);
    TCanvas *c = new TCanvas("c", "c", 400, 400);
    //TFile *f = TFile::Open("../data/run532.root");
    TFile *f = TFile::Open("../data/run546.root");
    TTree *genevents = (TTree*)f->Get("genevents");
    //TTree *events = (TTree*)f->Get("events");

    doTEffPlotR(c, genevents, "rrs", "Reco efficiency", "B0", 100, "", "hasCand==1");
    //doTEffPlotR(c, genevents, "rrs", "Reco efficiency", "B0", 100, "TMath::Abs(etamu1)<1.5&&TMath::Abs(etamu2)<1.5", "hasCand==1");
}

// =====================================================================================================
void doTEffPlotLb()
{
    const string p = "lb";
    //std::vector<double> bins = variableBinSizeVec(0, 12, 14, 16, 18, 20, 24, 28, 36, 60, 100);
    //TH1F *h1 = new TH1F("h1","h1", bins.size()-1, &bins[0]);
    //TH1F *h2 = new TH1F("h2","h2", bins.size()-1, &bins[0]);
    TCanvas *c = new TCanvas("c", "c", 400, 400);
    //TFile *f = TFile::Open("../data/run385__run387.root");
    //TFile *f = TFile::Open("../data/run413_415.root");
    //TFile *f = TFile::Open("../data/run456.root");
    TFile *f = TFile::Open("../data/run547.root");
    TTree *genevents = (TTree*)f->Get("genevents");
    TTree *events = (TTree*)f->Get("events");

    // ------------------------------------------------------------------------- step 1: acceptance
    // Muon acceptance
    const string acc_mu1_1 = "TMath::Abs(etamu1)<1.2&&ptmu1>3.5";
    const string acc_mu1_2 = "TMath::Abs(etamu1)>=1.2&&TMath::Abs(etamu1)<1.6&&ptmu1>8-3.75*TMath::Abs(etamu1)";
    const string acc_mu1_3 = "TMath::Abs(etamu1)>=1.6&&TMath::Abs(etamu1)<2.4&&ptmu1>2.0";
    const string acc_mu2_1 = "TMath::Abs(etamu2)<1.2&&ptmu2>3.5";
    const string acc_mu2_2 = "TMath::Abs(etamu2)>=1.2&&TMath::Abs(etamu2)<1.6&&ptmu2>8-3.75*TMath::Abs(etamu2)";
    const string acc_mu2_3 = "TMath::Abs(etamu2)>=1.6&&TMath::Abs(etamu2)<2.4&&ptmu2>2.0";
    const string acc_mu = "(("+acc_mu1_1+")||("+acc_mu1_2+")||("+acc_mu1_3+"))&&(("+acc_mu2_1+")||("+acc_mu2_2+")||("+acc_mu2_3+"))";

    //doTEffPlotSetGen(c, genevents, "eff_muAcc", "Step 1: #mu acceptance", p, 50, "", acc_mu);

    // ------------------------------------------------------------------------- step 2: Jpsi Kandidat
    const string jpCandAll = acc_mu;
    const string jpCandPass = "hasJpCand==1";

    ////doTEffPlotSetGen(c, genevents, "eff_JpCand", "Step 2: J/#psi candidate creation efficiency", p, 25, jpCandAll, jpCandPass);

    // ------------------------------------------------------------------------- step 3: L0 Kandidat
    const string l0CandAll = jpCandAll + "&&" + jpCandPass;
    const string l0CandPass = "hasrsCand==1";

    //doTEffPlotSetGen(c, genevents, "eff_rsCand", "Step 3: #Lambda^{0} candidate creation efficiency", p, 25, l0CandAll, l0CandPass);

    // ------------------------------------------------------------------------- step 4: Lb Kandidat
    const string lbCandAll = l0CandAll + "&&" + l0CandPass;
    const string lbCandPass = "hasCand==1";

    //doTEffPlotSetGen(c, genevents, "eff_LbCand", "Step 4: #Lambda_{b} candidate creation efficiency", p, 25, lbCandAll, lbCandPass);
    //doTEffPlotSetGen(c, genevents, "eff_LbCandOnly", "Step 4: #Lambda_{b} candidate only", p, 25, "", lbCandPass);

    // ------------------------------------------------------------------------- step 4a: Lb Kandidat mit Truthmatcher
    const string lbCandMatchAll = lbCandAll + "&&" + lbCandPass;
    const string lbCandMatchPass = "isMCmatch==1";

    doTEffPlotSetGen(c, genevents, "eff_LbCandMatch", "#Lambda_{b} candidate, truth matched", p, 25, lbCandAll, lbCandPass);
    doTEffPlotSetGen(c, genevents, "eff_LbCandMatchCum", "#Lambda_{b} candidate, truth matched, cum", p, 25, "", lbCandPass);

    // ------------------------------------------------------------------------- step 5: Trigger
    const string trigAll = lbCandAll + "&&" + lbCandPass;
    const string trigPassBarrel = "HLTmatch==1&&HLTokBarrelJpsi==1";
    const string trigPassDispl = "HLTmatch==1&&HLTokDisplJpsi==1";

    //doTEffPlotSetGen(c, genevents, "eff_HLTbarrEff", "Step 5: HLTBarrelJpsi efficiency", p, 50, trigAll, trigPassBarrel);
    //doTEffPlotSetGen(c, genevents, "eff_HLTdispEff", "Step 5: HLTDisplJpsi efficiency", p, 50, trigAll, trigPassDispl);

    // do the full from start for comparison
    //doTEffPlotSetGen(c, genevents, "eff_muAcc_HLTbarrEff", "Steps 1-5: Cumulative HLTBarrelJpsi efficiency", p, 50, "", trigAll+"&&"+trigPassBarrel);
    //doTEffPlotSetGen(c, genevents, "eff_muAcc_HLTdispEff", "Steps 1-5: Cumulative HLTDisplJpsi efficiency", p, 50, "", trigAll+"&&"+trigPassDispl);

    // and trigger only
    //doTEffPlotSetGen(c, genevents, "eff_only_HLTbarrEff", "Step 5 alone: HLTBarrelJpsi efficiency", p, 50, "", trigPassBarrel);
    //doTEffPlotSetGen(c, genevents, "eff_only_HLTdispEff", "Step 5 alone: HLTDisplJpsi efficiency", p, 50, "", trigPassDispl);

    // =====================================================================================================
    // ------------------------------------------------------------------------- step 6: analysis cuts
    Cuts cutLb;
    cutLb.selectCut("lb14", "acc06Lb", "muSoft");
    const string anaAll = "";
    const string anaPass = cutLb.getCut();
    const string anaPassBarrel = anaPass+"&&HLTokBarrelJpsi==1&&HLTmatch";
    const string anaPassDispl = anaPass+"&&HLTokDisplJpsi==1&&HLTmatch";

    const int nBinsReco(40);
    doTEffPlotSetReco(c, events, "eff_anal_notrig", "analysis cuts efficiency (no trig)", p, nBinsReco, anaAll, anaPass);
    doTEffPlotSetReco(c, events, "eff_anal_barr", "analysis cuts efficiency (barrel)", p, nBinsReco, anaAll, anaPassBarrel);
    doTEffPlotSetReco(c, events, "eff_anal_disp", "analysis cuts efficiency (displ)", p, nBinsReco, anaAll, anaPassDispl);

    doTEffPlotSetRecoRecoval(c, events, "eff_anal_barr", "analysis cuts efficiency (barrel)", p, nBinsReco, anaAll, anaPassBarrel);
    doTEffPlotSetRecoRecoval(c, events, "eff_anal_disp", "analysis cuts efficiency (displ)", p, nBinsReco, anaAll, anaPassDispl);

    // ------------------------------------------------------------------------- some additional plots
    Cuts cutMuSoft;
    cutMuSoft.selectCut("muSoft");
    const string anaMuSoft = cutMuSoft.getCut();
    //doTEffPlotSetReco(c, events, "eff_analMuSoft", "isol. soft #mu selection efficiency", p, nBinsReco, anaAll, anaMuSoft);

    Cuts cutAccLb;
    cutAccLb.selectCut("acc06Lb");
    const string anaAccLb = cutAccLb.getCut();
    //doTEffPlotSetReco(c, events, "eff_analAcc", "isol. reco acceptance efficiency", p, nBinsReco, anaAll, anaAccLb);
}

void doTEffPlotLbfinal()
{
    const string p = "lb";
    //TCanvas *c = new TCanvas("c", "c", 400, 400);
    TCanvas *c = new TCanvas("c", "c", 1000, 300);
    TFile *f = TFile::Open("../data/run547.root");
    TTree *genevents = (TTree*)f->Get("genevents");
    TTree *events = (TTree*)f->Get("events");

    // define binning
    const int nBins(50);
    const double lo(0), hi(10);
    std::vector<double> bins_t = makeDynamicBins(genevents, "ctbc", "hasCand==1", nBins, lo, hi, 1e12);

    // -------------------------------------------------------------------------

    //doTEffPlotSetGen(c, genevents, "eff_LbCandCum", "#Lambda_{b}+#bar{#Lambda}_{b} candidate creation from gen", p, nBins, "", "hasCand==1");
    //doTEffPlotSetGen(c, genevents, "eff_LbCandCum_lb", "#Lambda_{b} candidate creation from gen", p, nBins, "qha1>0", "qha1>0&&hasCand==1");
    //doTEffPlotSetGen(c, genevents, "eff_LbCandCum_lbbar", "#bar{#Lambda}_{b} candidate creation from gen", p, nBins, "qha1<0", "qha1<0&&hasCand==1");
    doTEffPlotSetGen_t(c, genevents, "eff_LbCandCum", "#Lambda_{b}+#bar{#Lambda}_{b} candidate creation from gen", p, "", "hasCand==1", bins_t);
    doTEffPlotSetGen_t(c, genevents, "eff_LbCandCum_lb", "#Lambda_{b} candidate creation from gen", p, "qha1>0", "qha1>0&&hasCand==1", bins_t);
    doTEffPlotSetGen_t(c, genevents, "eff_LbCandCum_lbbar", "#bar{#Lambda}_{b} candidate creation from gen", p, "qha1<0", "qha1<0&&hasCand==1", bins_t);

    // =====================================================================================================
    // ------------------------------------------------------------------------- step 6: analysis cuts
    Cuts cutLb;
    cutLb.selectCut("lb14", "acc06Lb", "muSoft");
    //const string anaAll = "";
    const string anaAll = "isSig==1";
    const string anaPass = cutLb.getCut();
    const string anaPassBarrel = anaPass+"&&HLTokBarrelJpsi==1&&HLTmatch";
    const string anaPassDispl = anaPass+"&&HLTokDisplJpsi==1&&HLTmatch";

    doTEffPlotSetReco_t(c, events, "eff_anal_notrig", "#Lambda_{b}+#bar{#Lambda}_{b} acc, muon, cuts from cand", p, anaAll, anaPass, bins_t);
    doTEffPlotSetReco_t(c, events, "eff_anal_barr", "#Lambda_{b}+#bar{#Lambda}_{b} barrel trg, acc, muon, cuts from cand", p, anaAll, anaPassBarrel, bins_t);
    doTEffPlotSetReco_t(c, events, "eff_anal_disp", "#Lambda_{b}+#bar{#Lambda}_{b} displ trg, acc, muon, cuts from cand", p, anaAll, anaPassDispl, bins_t);

    doTEffPlotSetReco_t(c, events, "eff_anal_barr_lb", "#Lambda_{b} barrel trg, acc, muon, cuts from cand", p, anaAll+"&&rqha1>0", anaPassBarrel+"&&rqha1>0", bins_t);
    doTEffPlotSetReco_t(c, events, "eff_anal_barr_lbbar", "#bar{#Lambda}_{b} barrel trg, acc, muon, cuts from cand", p, anaAll+"&&rqha1<0", anaPassBarrel+"&&rqha1<0", bins_t);

    //doTEffPlotSetRecoRecoval(c, events, "eff_anal_barr", "analysis cuts efficiency (barrel)", p, nBinsReco, anaAll, anaPassBarrel);
    //doTEffPlotSetRecoRecoval(c, events, "eff_anal_disp", "analysis cuts efficiency (displ)", p, nBinsReco, anaAll, anaPassDispl);

    doTEffPlotFit(c, "lb_eff_anal_barr_ct_pass", "lb_eff_LbCandCum_ct_all", "lb_eff_overall_ct", "#Lambda_{b}+#bar{#Lambda}_{b} overall efficiency", valueWithUnit("t("+getParticleLatex(p)+")","ps"));
    doTEffPlotFit(c, "lb_eff_anal_barr_lb_ct_pass", "lb_eff_LbCandCum_lb_ct_all", "lb_eff_overall_lb_ct", "#Lambda_{b} overall efficiency", valueWithUnit("t("+getParticleLatex(p)+")","ps"));
    doTEffPlotFit(c, "lb_eff_anal_barr_lbbar_ct_pass", "lb_eff_LbCandCum_lbbar_ct_all", "lb_eff_overall_lbbar_ct", "#bar{#Lambda}_{b} overall efficiency", valueWithUnit("t("+getParticleLatex(p)+")","ps"));
}


void doTEffPlotB0()
{
    const string p = "B0";
    //std::vector<double> bins = variableBinSizeVec(0, 12, 14, 16, 18, 20, 24, 28, 36, 60, 100);
    //TH1F *h1 = new TH1F("h1","h1", bins.size()-1, &bins[0]);
    //TH1F *h2 = new TH1F("h2","h2", bins.size()-1, &bins[0]);
    TCanvas *c = new TCanvas("c", "c", 400, 400);
    //TFile *f = TFile::Open("../data/run417_420.root");
    //TFile *f = TFile::Open("../data/run460_472.root");
    TFile *f = TFile::Open("../data/run546.root");
    TTree *genevents = (TTree*)f->Get("genevents");
    TTree *events = (TTree*)f->Get("events");

    // ------------------------------------------------------------------------- step 1: acceptance
    // Muon acceptance
    const string acc_mu1_1 = "TMath::Abs(etamu1)<1.2&&ptmu1>3.5";
    const string acc_mu1_2 = "TMath::Abs(etamu1)>=1.2&&TMath::Abs(etamu1)<1.6&&ptmu1>8-3.75*TMath::Abs(etamu1)";
    const string acc_mu1_3 = "TMath::Abs(etamu1)>=1.6&&TMath::Abs(etamu1)<2.4&&ptmu1>2.0";
    const string acc_mu2_1 = "TMath::Abs(etamu2)<1.2&&ptmu2>3.5";
    const string acc_mu2_2 = "TMath::Abs(etamu2)>=1.2&&TMath::Abs(etamu2)<1.6&&ptmu2>8-3.75*TMath::Abs(etamu2)";
    const string acc_mu2_3 = "TMath::Abs(etamu2)>=1.6&&TMath::Abs(etamu2)<2.4&&ptmu2>2.0";
    const string acc_mu = "(("+acc_mu1_1+")||("+acc_mu1_2+")||("+acc_mu1_3+"))&&(("+acc_mu2_1+")||("+acc_mu2_2+")||("+acc_mu2_3+"))";

    //doTEffPlotSetGen(c, genevents, "eff_muAcc", "Step 1: #mu acceptance", p, 50, "", acc_mu);

    // ------------------------------------------------------------------------- step 2: Jpsi Kandidat
    const string jpCandAll = acc_mu;
    const string jpCandPass = "hasJpCand==1";

    //doTEffPlotSetGen(c, genevents, "eff_JpCand", "Step 2:J/#psi candidate creation efficiency", p, 25, jpCandAll, jpCandPass);

    // ------------------------------------------------------------------------- step 3: Ks Kandidat
    const string KsCandAll = jpCandAll + "&&" + jpCandPass;
    const string KsCandPass = "hasrsCand==1";

    //doTEffPlotSetGen(c, genevents, "eff_rsCand", "Step 3: K_{s} candidate creation efficiency", p, 25, KsCandAll, KsCandPass);

    // ------------------------------------------------------------------------- step 4: B0 Kandidat
    const string b0CandAll = KsCandAll + "&&" + KsCandPass;
    const string b0CandPass = "hasCand==1";

    //doTEffPlotSetGen(c, genevents, "eff_B0Cand", "Step 4: B^{0} candidate creation efficiency", p, 25, b0CandAll, b0CandPass);

    // ------------------------------------------------------------------------- step 4a: B0 Kandidat mit Truthmatcher
    const string b0CandMatchAll = b0CandAll + "&&" + b0CandPass;
    const string b0CandMatchPass = "isMCmatch==1";

    doTEffPlotSetGen(c, genevents, "eff_B0CandMatch", "B^{0} candidate, truth matched", p, 25, b0CandAll, b0CandPass);
    doTEffPlotSetGen(c, genevents, "eff_B0CandMatchCum", "B^{0} candidate, truth matched, cum", p, 25, "", b0CandPass);

    // ------------------------------------------------------------------------- step 5: Trigger
    const string trigAll = b0CandAll + "&&" + b0CandPass;;
    const string trigPassBarrel = "HLTmatch==1&&HLTokBarrelJpsi==1";
    const string trigPassDispl = "HLTmatch==1&&HLTokDisplJpsi==1";

    //doTEffPlotSetGen(c, genevents, "eff_HLTbarrEff", "Step 5: HLTBarrelJpsi efficiency (barrel)", p, 50, trigAll, trigPassBarrel);
    //doTEffPlotSetGen(c, genevents, "eff_HLTdispEff", "Step 5: HLTDisplJpsi efficiency (displ)", p, 50, trigAll, trigPassDispl);

    // do the full from start for comparison
    //doTEffPlotSetGen(c, genevents, "eff_muAcc_HLTbarrEff", "Steps 1-5: Cumulative HLTBarrelJpsi efficiency", p, 50, "", trigPassBarrel);
    //doTEffPlotSetGen(c, genevents, "eff_muAcc_HLTdispEff", "Steps 1-5: Cumulative HLTDisplJpsi efficiency", p, 50, "", trigPassDispl);


    // =====================================================================================================
    // ------------------------------------------------------------------------- step 7: analysis cuts
    Cuts cutB0;
    cutB0.selectCut("B008", "acc06B0", "muSoft");
    const string anaAll = "isMCmatch==1&&isSig==1";
    const string anaPass = cutB0.getCut();
    const string anaPassBarrel = anaPass+"&&HLTokBarrelJpsi==1&&HLTmatch";
    const string anaPassDispl = anaPass+"&&HLTokDisplJpsi==1&&HLTmatch";

    doTEffPlotSetReco(c, events, "eff_anal_notrig", "analysis cuts efficiency (no trig)", p, 25, anaAll, anaPass);
    doTEffPlotSetReco(c, events, "eff_anal_barr", "analysis cuts efficiency (barrel)", p, 25, anaAll, anaPassBarrel);
    doTEffPlotSetReco(c, events, "eff_anal_disp", "analysis cuts efficiency (displ)", p, 25, anaAll, anaPassDispl);

    // ------------------------------------------------------------------------- some additional plots
    Cuts cutMuSoft;
    cutMuSoft.selectCut("muSoft");
    const string anaMuSoft = cutMuSoft.getCut();
    //doTEffPlotSetReco(c, events, "eff_analMuSoft", "isol. soft #mu selection efficiency", p, 25, anaAll, anaMuSoft);

    Cuts cutAccB0;
    cutAccB0.selectCut("acc06B0");
    const string anaAccB0 = cutAccB0.getCut();
    //doTEffPlotSetReco(c, events, "eff_analAcc", "isol. reco acceptance efficiency", p, 25, anaAll, anaAccB0);
}

void doTEffPlotB0final()
{
    const string p = "B0";
    //TCanvas *c = new TCanvas("c", "c", 400, 400);
    TCanvas *c = new TCanvas("c", "c", 1000, 300);
    TFile *f = TFile::Open("../data/run546.root");
    TTree *genevents = (TTree*)f->Get("genevents");
    TTree *events = (TTree*)f->Get("events");

    const int nBins(50);
    const double lo(0), hi(10);
    std::vector<double> bins_t = makeDynamicBins(genevents, "ctbc", "hasCand==1", nBins, lo, hi, 1e12);

    // -------------------------------------------------------------------------

    //doTEffPlotSetGen(c, genevents, "eff_B0CandMatch", "B^{0} candidate, truth matched", p, nBins, b0CandAll, b0CandPass);
    //doTEffPlotSetGen(c, genevents, "eff_B0CandCum", "B^{0} candidate creation from gen", p, nBins, "", "hasCand==1");
    //doTEffPlotSetGen(c, genevents, "eff_B0CandMatchCum", "B^{0} candidate, truth matched, cum", p, nBins, "", b0CandMatchPass+"&&"+b0CandPass);
    //doTEffPlotSetGen(c, genevents, "eff_B0CandMatchCumRad", "B^{0} candidate, truth matched, radcut, cum", p, nBins, "vrrs>3", b0CandMatchPass+"&&"+b0CandPass);
    doTEffPlotSetGen_t(c, genevents, "eff_B0CandCum", "B^{0} candidate creation from gen", p, "", "hasCand==1", bins_t);


    // =====================================================================================================
    // ------------------------------------------------------------------------- step 7: analysis cuts
    Cuts cutB0;
    cutB0.selectCut("B008", "acc06B0", "muSoft");
    //const string anaAll = "isMCmatch==1&&isSig==1";
    const string anaAll = "isSig==1";
    const string anaPass = cutB0.getCut();
    const string anaPassBarrel = anaPass+"&&HLTokBarrelJpsi==1&&HLTmatch";
    const string anaPassDispl = anaPass+"&&HLTokDisplJpsi==1&&HLTmatch";

    doTEffPlotSetReco_t(c, events, "eff_anal_notrig", "B^{0} acc, muon, cuts from cand", p, anaAll, anaPass, bins_t);
    doTEffPlotSetReco_t(c, events, "eff_anal_barr", "B^{0} barrel trg, acc, muon, cuts from cand", p, anaAll, anaPassBarrel, bins_t);
    doTEffPlotSetReco_t(c, events, "eff_anal_disp", "B^{0} displ trg, acc, muon, cuts from cand", p, anaAll, anaPassDispl, bins_t);

    doTEffPlotFit(c, "B0_eff_anal_barr_ct_pass", "B0_eff_B0CandCum_ct_all", "B0_eff_overall_ct", "B^{0} overall efficiency", valueWithUnit("t("+getParticleLatex(p)+")","ps"));
}


