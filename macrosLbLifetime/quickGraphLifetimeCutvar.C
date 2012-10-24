#include <string>

#include "TROOT.h"
#include "TPad.h"
#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLine.h"
#include "TMultiGraph.h"

#include "utils.h"

TGraph* quickGraphLifetimeCutvar(string filename, string title, double tau, string pdfFile)
{
    const double xLo(-.5), xHi(4.5);
    TGraphErrors *tgr = new TGraphErrors(filename.c_str());
    tgr->SetMarkerStyle(21);
    tgr->SetTitle(title.c_str());
    tgr->Draw("AP");
    tgr->GetXaxis()->SetTitle("Cut on r(V^{0}) (cm)");
    tgr->GetXaxis()->SetLimits(xLo, xHi);
    tgr->GetYaxis()->SetTitle("Lifetime (ps)");
    TLine *line = new TLine(xLo, tau, xHi, tau);
    line->Draw();
    line->SetLineColor(3);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    gPad->SaveAs(pdfFile.c_str());
    return tgr;
}

void doSomeGraphs()
{
    const double tauB0(1.536);
    const double tauLb(1.424);
    quickGraphLifetimeCutvar("b0_journal_p102_all.dat", "B^{0} truth required a cand", tauB0, "b0_rrvscuts_all.pdf");
    quickGraphLifetimeCutvar("b0_journal_p102_cut.dat", "B^{0} truth, cand+cuts", tauB0, "b0_rrvscuts_cut.pdf");

    TMultiGraph *mgr_lb = new TMultiGraph();
    mgr_lb->Add(quickGraphLifetimeCutvar("lbboth_journal_p103_all.dat", "#Lambda_{b}+#bar{#Lambda_{b}} truth required a cand", tauLb, "lball_rrvscuts_all.pdf"));
    mgr_lb->Add(quickGraphLifetimeCutvar("lb_journal_p103_all.dat", "#Lambda_{b} truth required a cand", tauLb, "lb_rrvscuts_all.pdf"));
    mgr_lb->Add(quickGraphLifetimeCutvar("lbbar_journal_p103_all.dat", "#bar{#Lambda_{b}} truth required a cand", tauLb, "lbbar_rrvscuts_all.pdf"));
    mgr_lb->Draw("AP");
}

