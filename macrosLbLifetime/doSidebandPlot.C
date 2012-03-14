#include "TROOT.h"
#include "TH1F.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "utils.h"
#include <string>
#include <algorithm>
#include <iostream>

using std::endl;
using std::cout;

using std::string;

void printH(TH1F *h)
{
    cout << h->GetName() << ": ";
    for(int i=1; i<=h->GetNbinsX(); i++)
	cout << h->GetBinContent(i) << " ";
    cout << endl;
}

void doSidebandPlot(TTree *tree, string name, string toPlot, string title, string xtitle, string xunit, int nBins, double min, double max, string cutSignal, string cutSideband, double nSig, double nBgr, double drawLine1 = -9999, double drawLine2 = -9999)
{
    cout << name << endl;
    gPad->SetLeftMargin(.12);
    gPad->SetRightMargin(.08);
    gPad->SetTopMargin(.08);
    gPad->SetBottomMargin(.12);

    string histSigName = name + "_sig";
    string histSideName = name + "_sideband";
    string histSideCorrName = name + "_sidebandcorr";
    string histformat = "(" + toString(nBins) + "," + toString(min) + "," + toString(max) + ")";

    // draw signal region
    tree->Draw((toPlot + ">>" + histSigName + histformat).c_str(), cutSignal.c_str());
    TH1F *hsig = (TH1F*)gDirectory->Get(histSigName.c_str());
    hsig->SetTitle(title.c_str());
    hsig->GetXaxis()->SetTitle(valueWithUnit(xtitle, xunit).c_str());
    const double binsize = (max-min)/nBins;
    hsig->GetYaxis()->SetTitle(entriesPerBin(binsize, xunit).c_str());
    hsig->GetYaxis()->SetTitleOffset(0.90);
    hsig->SetLineStyle(2);
    hsig->SetLineColor(9);
    hsig->SetFillStyle(0);

    // draw sideband region
    tree->Draw((toPlot + ">>" + histSideName + histformat).c_str(), cutSideband.c_str(),"same");
    TH1F *hside = (TH1F*)gDirectory->Get(histSideName.c_str());
    hside->SetTitle("");
    hside->SetLineStyle(2);
    hside->SetLineColor(50);
    hside->SetFillStyle(0);

    // draw sidebandcorrected
    TH1F *hsidecorr = new TH1F(histSideCorrName.c_str(), title.c_str(), nBins, min, max);
    const int Nhsig = hsig->GetEntries();
    const int Nhside = hside->GetEntries();
    const double weightBgr = Nhsig/Nhside*nBgr/(nSig+nBgr); // see p. 126 in Journal #13
    hsidecorr->Add(hsig, hside, 1.0, -weightBgr);
    hsidecorr->SetTitle("");
    hsidecorr->SetLineColor(4);
    hsidecorr->SetFillStyle(0);

    // scale sideband histo to signal histo
    hside->Scale(weightBgr);

    // determine max value of y axis
    double curMax = 1.2 * std::max(std::max(hsig->GetBinContent(hsig->GetMaximumBin()), hside->GetBinContent(hside->GetMaximumBin())),
	    hsidecorr->GetBinContent(hsidecorr->GetMaximumBin()));
    hsig->SetMaximum(curMax);
    hsig->SetMinimum(0.0);
    hside->SetMaximum(curMax);
    hside->SetMinimum(0.0);
    hsidecorr->SetMaximum(curMax);
    hsidecorr->SetMinimum(0.0);

    // final draw of histos
    hsig->Draw();
    hside->Draw("same");
    hsidecorr->Draw("same");

    // Add legend
    TLegend *legend = new TLegend(0.55,0.75,0.90,0.90);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    legend->AddEntry(hsidecorr,"Sidebandsubtracted data","l");
    legend->AddEntry(hsig,"Signal region, unscaled","l");
    legend->AddEntry(hside,"Sideband, scaled to bgr","l");
    legend->Draw();

    // draw lines if necessary
    if (drawLine1 != -9999)
    {
	const double length = hsig->GetMaximum() * .7;
	TLine *l1 = new TLine(drawLine1, hsig->GetMinimum(), drawLine1, hsig->GetMinimum()+length);
	l1->SetLineColor(2);
	l1->SetLineWidth(2);
	l1->Draw();
    }
    if (drawLine2 != -9999)
    {
	const double length = hsig->GetMaximum() * .7;
	TLine *l2 = new TLine(drawLine2, hsig->GetMinimum(), drawLine2, hsig->GetMinimum()+length);
	l2->SetLineColor(2);
	l2->SetLineWidth(2);
	l2->Draw();
    }

    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->Update();
}

void doSidebandPlot(TTree *datatree, TTree *mctree, string name, string toPlot, string title, string xtitle, string xunit, int nBins, double min, double max, 
	string cutSignal, string cutSideband, string cutSignalMC, double nSig, double nBgr, double drawLine1 = -9999, double drawLine2 = -9999)
{
    gPad->SetLeftMargin(.12);
    gPad->SetRightMargin(.08);
    gPad->SetTopMargin(.08);
    gPad->SetBottomMargin(.12);

    string histSigName = name + "_sig";
    string histSideName = name + "_sideband";
    string histSideCorrName = name + "_sidebandcorr";
    string histMCName = name + "_mc";

    string histformat = "(" + toString(nBins) + "," + toString(min) + "," + toString(max) + ")";

    // draw signal region
    datatree->Draw((toPlot + ">>" + histSigName + histformat).c_str(), cutSignal.c_str());
    TH1F *hsig = (TH1F*)gDirectory->Get(histSigName.c_str());
    hsig->SetTitle(title.c_str());
    hsig->GetXaxis()->SetTitle(valueWithUnit(xtitle, xunit).c_str());
    const double binsize = (max-min)/nBins;
    hsig->GetYaxis()->SetTitle(entriesPerBin(binsize, xunit).c_str());
    hsig->GetYaxis()->SetTitleOffset(0.90);
    hsig->SetLineStyle(2);
    hsig->SetLineColor(9);
    hsig->SetFillStyle(0);
    cout << "hsig has " << hsig->GetEntries() << " entries" << endl;

    // draw sideband region
    datatree->Draw((toPlot + ">>" + histSideName + histformat).c_str(), cutSideband.c_str(),"same");
    TH1F *hside = (TH1F*)gDirectory->Get(histSideName.c_str());
    hside->SetTitle("");
    hside->SetLineStyle(2);
    hside->SetLineColor(50);
    hside->SetFillStyle(0);
    cout << "hside has " << hside->GetEntries() << " entries" << endl;

    // draw sidebandcorrected
    TH1F *hsidecorr = new TH1F(histSideCorrName.c_str(), title.c_str(), nBins, min, max);
    const double Nhsig = hsig->GetEntries();
    const double Nhside = hside->GetEntries();
    const double weightBgr = Nhsig/Nhside*nBgr/(nSig+nBgr); // see p. 126 in Journal #13
    hsidecorr->Add(hsig, hside, 1.0, -weightBgr);
    hsidecorr->SetTitle("");
    hsidecorr->SetLineColor(4);
    hsidecorr->SetLineWidth(2);
    hsidecorr->SetFillStyle(0);

    // scale sideband histo to signal histo
    hside->Scale(weightBgr);

    // draw MC
    cout << "mctree has " << mctree->GetEntries() << " entries" << endl;
    mctree->Draw((toPlot + ">>" + histMCName + histformat).c_str(), cutSignalMC.c_str(), "same");
    TH1F *hmc = (TH1F*)gDirectory->Get(histMCName.c_str());
    hmc->SetTitle("");
    hmc->SetLineStyle(0);
    hmc->SetLineWidth(2);
    hmc->SetLineColor(8);
    hmc->SetFillStyle(0);

    const double Nhmc = hmc->GetEntries();
    if (Nhmc != 0) hmc->Scale(Nhsig/Nhmc*nSig/(nSig+nBgr));
    else cout << "WARNING: histo for MC is empty!" << endl;

    // determine max value of y axis
    double curMax = 1.2 * std::max(std::max(hsig->GetBinContent(hsig->GetMaximumBin()), hside->GetBinContent(hside->GetMaximumBin())),
	    hsidecorr->GetBinContent(hsidecorr->GetMaximumBin()));
    hsig->SetMaximum(curMax);
    hsig->SetMinimum(0.0);
    hside->SetMaximum(curMax);
    hside->SetMinimum(0.0);
    hsidecorr->SetMaximum(curMax);
    hsidecorr->SetMinimum(0.0);
    hmc->SetMaximum(curMax);
    hmc->SetMinimum(0.0);

    // set error bar
    // following https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
    // but open issues on error propagation : https://hypernews.cern.ch/HyperNews/CMS/get/statistics/225.html
    //hsig->SetBinErrorOption(TH1::kPoisson);
    //hside->SetBinErrorOption(TH1::kPoisson);
    //hsidecorr->SetBinErrorOption(TH1::kPoisson);
    //hsig->SetBinErrorOption(TH1::kPoisson);
    //hsig->SetMarkerStyle(20);

    // final draw of histos
    hsig->Draw();
    hside->Draw("same");
    hsidecorr->Draw("same");
    hmc->Draw("same");

    // Add legend
    TLegend *legend = new TLegend(0.55,0.75,0.90,0.90);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    legend->AddEntry(hsidecorr,"Sidebandsubtracted data","l");
    legend->AddEntry(hsig,"Signal region, unscaled","l");
    legend->AddEntry(hside,"Sideband, scaled to bgr","l");
    legend->AddEntry(hmc,"MC, scaled to signal","l");
    legend->Draw();

    // draw lines if necessary
    {
	const int linecol(2);
	const double length = hsig->GetMaximum() * .7;
	const double width = 0.1*(hsig->GetBinCenter(hsig->GetNbinsX()) - hsig->GetBinCenter(0));
	const double hminX = 0.5*(hsig->GetBinCenter(0)+hsig->GetBinCenter(1));
	const double hmaxX = 0.5*(hsig->GetBinCenter(hsig->GetNbinsX())+hsig->GetBinCenter(hsig->GetNbinsX()+1));
	const double hminY = hsig->GetMinimum();
	if (drawLine1 != -9999)
	{
	    // draw a line
	    TLine *l1 = new TLine(drawLine1, hminY, drawLine1, hminY+length);
	    l1->SetLineColor(linecol);
	    l1->SetLineWidth(linecol);
	    l1->Draw();
	    // draw shaded box left to it
	    double left = drawLine1-width;
	    if (drawLine1 < hminX) left = hminX;
	    if (drawLine2 != -9999 && drawLine2<drawLine1) left = 0.5*(drawLine2+drawLine1);
	    TBox *b1 = new TBox(left, hminY, drawLine1, hminY+length);
	    b1->SetLineColor(linecol);
	    b1->SetFillColor(linecol);
	    b1->SetFillStyle(3004);
	    b1->Draw();
	}
	if (drawLine2 != -9999)
	{
	    // draw a line
	    TLine *l2 = new TLine(drawLine2, hminY, drawLine2, hminY+length);
	    l2->SetLineColor(2);
	    l2->SetLineWidth(2);
	    l2->Draw();
	    // draw shaded box right to it
	    double right = drawLine2+width;
	    if (drawLine2 > hmaxX) right = hmaxX;
	    if (drawLine1 != -9999 && drawLine2<drawLine1) right = 0.5*(drawLine2+drawLine1);
	    TBox *b1 = new TBox(drawLine2, hminY, right, hminY+length);
	    b1->SetLineColor(linecol);
	    b1->SetFillColor(linecol);
	    b1->SetFillStyle(3004);
	    b1->Draw();
	}
    }

    gPad->SetLogy(0);
    gPad->SetLogx(0);
    gPad->Update();

}

void doSidebandPlot(TTree *tree, string name, string toPlot, string title, string xtitle, string xunit, int nBins, double min, double max,
	string curCut, string cutSignal, string cutSideband, double nSig, double nSide, double drawLine1 = -9999, double drawLine2 = -9999)
{
    doSidebandPlot(tree, name, toPlot, title, xtitle, xunit, nBins, min, max, curCut+"&&"+cutSignal, curCut+"&&"+cutSideband, nSig, nSide, drawLine1, drawLine2);
}

void doSidebandPlot(TTree *datatree, TTree *mctree, string name, string toPlot, string title, string xtitle, string xunit, int nBins, double min, double max,
	string cutData, string cutMC, string cutSignal, string cutSideband, double nSig, double nBgr, double drawLine1 = -9999, double drawLine2 = -9999)
{
    doSidebandPlot(datatree, mctree, name, toPlot, title, xtitle, xunit, nBins, min, max, cutData+"&&"+cutSignal, cutData+"&&"+cutSideband, cutMC+"&&"+cutSignal, nSig, nBgr, drawLine1, drawLine2);
}

void doSidebandPlotLogY(TTree *tree, string name, string toPlot, string title, string xtitle, string xunit, int nBins, double min, double max, string cutSignal, string cutSideband, double nSig, double nSide)
{
    doSidebandPlot(tree, name, toPlot, title, xtitle, xunit, nBins, min, max, cutSignal, cutSideband, nSig, nSide);
    gPad->SetLogy();
}

void doSidebandPlotLogY(TTree *tree, string name, string toPlot, string title, string xtitle, string xunit, int nBins, double min, double max, string curCut, string cutSignal, string cutSideband, double nSig, double nSide)
{
    doSidebandPlotLogY(tree, name, toPlot, title, xtitle, xunit, nBins, min, max, curCut+"&&"+cutSignal, curCut+"&&"+cutSideband, nSig, nSide);
}

