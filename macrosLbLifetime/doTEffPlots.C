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

#include "utils.h"
#include "Cuts.C"

using std::string;
using std::auto_ptr;

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
    TF1* f1 = new TF1("f1","[0]*(1+[1]*x)",lo,hi);
    f1->SetParameters(averageEff, -0.03e12);
    heff->Fit(f1);

    c->SaveAs((hname+".pdf").c_str());

    cout << hname << ": " << hpass->GetEntries() << "/" << hall->GetEntries() << ": " << averageEff << endl;
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
void doTEffPlotSetGen(TCanvas *c, TTree *t, string hname, string title, string particle, int nBins, const string &cutAll, const string &cutPass)
{
    doTEffPlot(c, t, particle+"_"+hname+"_p", title, "pbc", nBins, 0, 100, valueWithUnit("p("+getParticleLatex(particle)+")","GeV/c"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_d3", title, "d3dbc", nBins, 0, 1, valueWithUnit("d_{3}("+getParticleLatex(particle)+")","cm"), cutAll, cutPass);
    doTEffPlotFit(c, t, particle+"_"+hname+"_ct", title, "ctbc", nBins, 0, 10e-12, valueWithUnit("t("+getParticleLatex(particle)+")","s"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_ctz", title+" (zoom)", "ctbc", .4*nBins, 0, 2e-12, valueWithUnit("t("+getParticleLatex(particle)+")","s"), cutAll, cutPass);
}

void doTEffPlotSetReco(TCanvas *c, TTree *t, string hname, string title, string particle, int nBins, const string &cutAll, const string &cutPass)
{
    doTEffPlot(c, t, particle+"_"+hname+"_p", title, "pbctruth", nBins, 0, 100, valueWithUnit("p("+getParticleLatex(particle)+")","GeV/c"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_d3", title, "d3dbctruth", nBins, 0, 1, valueWithUnit("d_{3}("+getParticleLatex(particle)+")","cm"), cutAll, cutPass);
    doTEffPlotFit(c, t, particle+"_"+hname+"_ct", title, "ctbctruth", nBins, 0, 10e-12, valueWithUnit("t("+getParticleLatex(particle)+")","s"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_ctz", title+" (zoom)", "ctbctruth", .4*nBins, 0, 2e-12, valueWithUnit("t("+getParticleLatex(particle)+")","s"), cutAll, cutPass);
}

void doTEffPlotSetRecoRecoval(TCanvas *c, TTree *t, string hname, string title, string particle, int nBins, const string &cutAll, const string &cutPass)
{
    doTEffPlot(c, t, particle+"_"+hname+"_reco_p", title, "pbc", nBins, 0, 100, valueWithUnit("p("+getParticleLatex(particle)+") reco","GeV/c"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_reco_d3", title, "d3bc", nBins, -0.2, 1, valueWithUnit("d_{3}("+getParticleLatex(particle)+") reco","cm"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_reco_ct", title, "ct3dbc", nBins, -2e-12, 10e-12, valueWithUnit("t("+getParticleLatex(particle)+") reco","s"), cutAll, cutPass);
    doTEffPlot(c, t, particle+"_"+hname+"_reco_ctz", title+" (zoom)", "ct3dbc", nBins, -.5e-12, 2e-12, valueWithUnit("t("+getParticleLatex(particle)+") reco","s"), cutAll, cutPass);
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
    TFile *f = TFile::Open("../data/run456.root");
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

    doTEffPlotSetGen(c, genevents, "eff_muAcc", "Step 1: #mu acceptance", p, 50, "", acc_mu);

    // ------------------------------------------------------------------------- step 2: Jpsi Kandidat
    const string jpCandAll = acc_mu;
    const string jpCandPass = "hasJpCand==1";

    doTEffPlotSetGen(c, genevents, "eff_JpCand", "Step 2: J/#psi candidate creation efficiency", p, 25, jpCandAll, jpCandPass);

    // ------------------------------------------------------------------------- step 3: L0 Kandidat
    const string l0CandAll = jpCandAll + "&&" + jpCandPass;
    const string l0CandPass = "hasrsCand==1";

    doTEffPlotSetGen(c, genevents, "eff_rsCand", "Step 3: #Lambda^{0} candidate creation efficiency", p, 25, l0CandAll, l0CandPass);

    // ------------------------------------------------------------------------- step 4: Lb Kandidat
    const string lbCandAll = l0CandAll + "&&" + l0CandPass;
    const string lbCandPass = "hasCand==1";

    doTEffPlotSetGen(c, genevents, "eff_LbCand", "Step 4: #Lambda_{b} candidate creation efficiency", p, 25, lbCandAll, lbCandPass);
    doTEffPlotSetGen(c, genevents, "eff_LbCandOnly", "Step 4: #Lambda_{b} candidate only", p, 25, "", lbCandPass);

    // ------------------------------------------------------------------------- step 4a: Lb Kandidat mit Truthmatcher
    const string lbCandMatchAll = lbCandAll + "&&" + lbCandPass;
    const string lbCandMatchPass = "isMCmatch==1";

    doTEffPlotSetGen(c, genevents, "eff_LbCandMatch", "Step 4b: #Lambda_{b} candidate, truth matched", p, 25, lbCandAll, lbCandPass);
    doTEffPlotSetGen(c, genevents, "eff_LbCandMatchCum", "Step 4b: #Lambda_{b} candidate, truth matched, cum", p, 25, "", lbCandPass);

    // ------------------------------------------------------------------------- step 5: Trigger
    const string trigAll = lbCandAll + "&&" + lbCandPass;
    const string trigPassBarrel = "HLTmatch==1&&HLTokBarrelJpsi==1";
    const string trigPassDispl = "HLTmatch==1&&HLTokDisplJpsi==1";

    doTEffPlotSetGen(c, genevents, "eff_HLTbarrEff", "Step 5: HLTBarrelJpsi efficiency", p, 50, trigAll, trigPassBarrel);
    doTEffPlotSetGen(c, genevents, "eff_HLTdispEff", "Step 5: HLTDisplJpsi efficiency", p, 50, trigAll, trigPassDispl);

    // do the full from start for comparison
    doTEffPlotSetGen(c, genevents, "eff_muAcc_HLTbarrEff", "Steps 1-5: Cumulative HLTBarrelJpsi efficiency", p, 50, "", trigAll+"&&"+trigPassBarrel);
    doTEffPlotSetGen(c, genevents, "eff_muAcc_HLTdispEff", "Steps 1-5: Cumulative HLTDisplJpsi efficiency", p, 50, "", trigAll+"&&"+trigPassDispl);

    // and trigger only
    doTEffPlotSetGen(c, genevents, "eff_only_HLTbarrEff", "Step 5 alone: HLTBarrelJpsi efficiency", p, 50, "", trigPassBarrel);
    doTEffPlotSetGen(c, genevents, "eff_only_HLTdispEff", "Step 5 alone: HLTDisplJpsi efficiency", p, 50, "", trigPassDispl);

    // =====================================================================================================
    // ------------------------------------------------------------------------- step 6: analysis cuts
    Cuts cutLb;
    cutLb.selectCut("lb12", "acc06Lb", "muSoft");
    const string anaAll = "";
    const string anaPass = cutLb.getCut();
    const string anaPassBarrel = anaPass+"&&HLTokBarrelJpsi==1&&HLTmatch";
    const string anaPassDispl = anaPass+"&&HLTokDisplJpsi==1&&HLTmatch";

    const int nBinsReco(40);
    doTEffPlotSetReco(c, events, "eff_anal_notrig", "Step 6: analysis cuts efficiency (no trig)", p, nBinsReco, anaAll, anaPass);
    doTEffPlotSetReco(c, events, "eff_anal_barr", "Step 6: analysis cuts efficiency (barrel)", p, nBinsReco, anaAll, anaPassBarrel);
    doTEffPlotSetReco(c, events, "eff_anal_disp", "Step 6: analysis cuts efficiency (displ)", p, nBinsReco, anaAll, anaPassDispl);

    doTEffPlotSetRecoRecoval(c, events, "eff_anal_barr", "Step 6: analysis cuts efficiency (barrel)", p, nBinsReco, anaAll, anaPassBarrel);
    doTEffPlotSetRecoRecoval(c, events, "eff_anal_disp", "Step 6: analysis cuts efficiency (displ)", p, nBinsReco, anaAll, anaPassDispl);

    // ------------------------------------------------------------------------- some additional plots
    Cuts cutMuSoft;
    cutMuSoft.selectCut("muSoft");
    const string anaMuSoft = cutMuSoft.getCut();
    doTEffPlotSetReco(c, events, "eff_analMuSoft", "isol. soft #mu selection efficiency", p, nBinsReco, anaAll, anaMuSoft);

    Cuts cutAccLb;
    cutAccLb.selectCut("acc06Lb");
    const string anaAccLb = cutAccLb.getCut();
    doTEffPlotSetReco(c, events, "eff_analAcc", "isol. reco acceptance efficiency", p, nBinsReco, anaAll, anaAccLb);
}

void doTEffPlotB0()
{
    const string p = "B0";
    //std::vector<double> bins = variableBinSizeVec(0, 12, 14, 16, 18, 20, 24, 28, 36, 60, 100);
    //TH1F *h1 = new TH1F("h1","h1", bins.size()-1, &bins[0]);
    //TH1F *h2 = new TH1F("h2","h2", bins.size()-1, &bins[0]);
    TCanvas *c = new TCanvas("c", "c", 400, 400);
    //TFile *f = TFile::Open("../data/run417_420.root");
    TFile *f = TFile::Open("../data/run460.root");
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

    doTEffPlotSetGen(c, genevents, "eff_muAcc", "Step 1: #mu acceptance", p, 50, "", acc_mu);

    // ------------------------------------------------------------------------- step 2: Jpsi Kandidat
    const string jpCandAll = acc_mu;
    const string jpCandPass = "hasJpCand==1";

    doTEffPlotSetGen(c, genevents, "eff_JpCand", "Step 2:J/#psi candidate creation efficiency", p, 25, jpCandAll, jpCandPass);

    // ------------------------------------------------------------------------- step 3: Ks Kandidat
    const string KsCandAll = jpCandAll + "&&" + jpCandPass;
    const string KsCandPass = "hasrsCand==1";

    doTEffPlotSetGen(c, genevents, "eff_rsCand", "Step 3: K_{s} candidate creation efficiency", p, 25, KsCandAll, KsCandPass);

    // ------------------------------------------------------------------------- step 4: B0 Kandidat
    const string b0CandAll = KsCandAll + "&&" + KsCandPass;
    const string b0CandPass = "hasCand==1";

    doTEffPlotSetGen(c, genevents, "eff_B0Cand", "Step 4: B^{0} candidate creation efficiency", p, 25, b0CandAll, b0CandPass);

    // ------------------------------------------------------------------------- step 4a: B0 Kandidat mit Truthmatcher
    const string b0CandMatchAll = b0CandAll + "&&" + b0CandPass;
    const string b0CandMatchPass = "isMCmatch==1";

    doTEffPlotSetGen(c, genevents, "eff_B0CandMatch", "Step 4b: B^{0} candidate, truth matched", p, 25, b0CandAll, b0CandPass);
    doTEffPlotSetGen(c, genevents, "eff_B0CandMatchCum", "Step 4b: B^{0} candidate, truth matched, cum", p, 25, "", b0CandPass);

    // ------------------------------------------------------------------------- step 5: Trigger
    const string trigAll = b0CandAll + "&&" + b0CandPass;;
    const string trigPassBarrel = "HLTmatch==1&&HLTokBarrelJpsi==1";
    const string trigPassDispl = "HLTmatch==1&&HLTokDisplJpsi==1";

    doTEffPlotSetGen(c, genevents, "eff_HLTbarrEff", "Step 5: HLTBarrelJpsi efficiency (barrel)", p, 50, trigAll, trigPassBarrel);
    doTEffPlotSetGen(c, genevents, "eff_HLTdispEff", "Step 5: HLTDisplJpsi efficiency (displ)", p, 50, trigAll, trigPassDispl);

    // do the full from start for comparison
    doTEffPlotSetGen(c, genevents, "eff_muAcc_HLTbarrEff", "Steps 1-5: Cumulative HLTBarrelJpsi efficiency", p, 50, "", trigPassBarrel);
    doTEffPlotSetGen(c, genevents, "eff_muAcc_HLTdispEff", "Steps 1-5: Cumulative HLTDisplJpsi efficiency", p, 50, "", trigPassDispl);


    // =====================================================================================================
    // ------------------------------------------------------------------------- step 7: analysis cuts
    Cuts cutB0;
    cutB0.selectCut("B006", "acc06B0", "muSoft");
    const string anaAll = "isMCmatch==1&&isSig==1";
    const string anaPass = cutB0.getCut();
    const string anaPassBarrel = anaPass+"&&HLTokBarrelJpsi==1&&HLTmatch";
    const string anaPassDispl = anaPass+"&&HLTokDisplJpsi==1&&HLTmatch";

    doTEffPlotSetReco(c, events, "eff_anal_notrig", "Step 6: analysis cuts efficiency (no trig)", p, 25, anaAll, anaPass);
    doTEffPlotSetReco(c, events, "eff_anal_barr", "Step 6: analysis cuts efficiency (barrel)", p, 25, anaAll, anaPassBarrel);
    doTEffPlotSetReco(c, events, "eff_anal_disp", "Step 6: analysis cuts efficiency (displ)", p, 25, anaAll, anaPassDispl);

    // ------------------------------------------------------------------------- some additional plots
    Cuts cutMuSoft;
    cutMuSoft.selectCut("muSoft");
    const string anaMuSoft = cutMuSoft.getCut();
    doTEffPlotSetReco(c, events, "eff_analMuSoft", "isol. soft #mu selection efficiency", p, 25, anaAll, anaMuSoft);

    Cuts cutAccB0;
    cutAccB0.selectCut("acc06B0");
    const string anaAccB0 = cutAccB0.getCut();
    doTEffPlotSetReco(c, events, "eff_analAcc", "isol. reco acceptance efficiency", p, 25, anaAll, anaAccB0);
}

