#ifndef EFFICIENCYPLUTS_H_GUARD
#define EFFICIENCYPLUTS_H_GUARD

#include <utility>
#include <vector>
#include "TLegend.h"

// A collection of plot routines mainly for efficiency plots

void doPlot2d(TPad *pad, TTree *t, std::string hname, std::string todraw, 
	int nBinsX, double minX, double maxX,
	int nBinsY, double minY, double maxY,
        std::string title, std::string titleX, std::string titleY,
	std::string unitX, std::string unitY)
{
    const std::string plotstring = todraw + ">>" + hname + "("
	+ toString(nBinsX) + "," + toString(minX) + "," + toString(maxX) + ","
	+ toString(nBinsY) + "," + toString(minY) + "," + toString(maxY) + ")";
    cout << plotstring << endl;
    t->Draw(plotstring.c_str(),"","COLZ");
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(hname.c_str());
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h->GetYaxis()->SetTitle(unitY.size()>0 ? (titleY+" / "+unitY).c_str() : titleY.c_str());
    //const double binsize = (max-min)/nBins;
    //h->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    return; 
}

void doPlot2d(TPad *pad, TTree *t, std::string hname, std::string todraw, std::string cutstring,
	int nBinsX, double minX, double maxX,
	int nBinsY, double minY, double maxY,
        std::string title, std::string titleX, std::string titleY,
	std::string unitX, std::string unitY)
{
    const std::string plotstring = todraw + ">>" + hname + "("
	+ toString(nBinsX) + "," + toString(minX) + "," + toString(maxX) + ","
	+ toString(nBinsY) + "," + toString(minY) + "," + toString(maxY) + ")";
    cout << plotstring << endl;
    // Set up the first pad
    pad->cd();
    pad->SetTopMargin(0.15);
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(0.20);
    pad->SetRightMargin(0.30);
    pad->cd();
    // now draw the histo
    t->Draw(plotstring.c_str(),cutstring.c_str(),"COLZ");
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(hname.c_str());
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h->GetYaxis()->SetTitle(unitY.size()>0 ? (titleY+" / "+unitY).c_str() : titleY.c_str());
    h->DrawCopy("COLZ");
    //const double binsize = (max-min)/nBins;
    //h->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    return; 
}

void doPlot2d(TPad *pad, TTree *t, std::string hname, std::string todraw, std::string cutstring,
	std::vector<double> binsX,
	std::vector<double> binsY,
        std::string title, std::string titleX, std::string titleY,
	std::string unitX, std::string unitY)
{
    const std::string plotstring = todraw + ">>" + hname;
    // Set up the first pad
    pad->cd();
    pad->SetTopMargin(0.10);
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0.30);
    pad->cd();
    // now draw the histo
    TH2F *h = new TH2F(hname.c_str(),title.c_str(),binsX.size()-1,&binsX[0],binsY.size()-1,&binsY[0]);
    t->Draw(plotstring.c_str(),cutstring.c_str(),"COLZ");
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h->GetYaxis()->SetTitle(unitY.size()>0 ? (titleY+" / "+unitY).c_str() : titleY.c_str());
    h->DrawCopy("COLZ");
    //const double binsize = (max-min)/nBins;
    //h->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    return; 
}

// Divides the 2nd histo by the 1st
void doPlotRatio2d(TPad *pad, std::string htitle, std::string hname1, std::string hname2)
{
    TH2F *h1 = (TH2F*)gDirectory->GetList()->FindObject(hname1.c_str());
    TH2F *h2 = (TH2F*)gDirectory->GetList()->FindObject(hname2.c_str());
    pad->cd();
    pad->SetTopMargin(0.10);
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0.15);
    pad->SetLogz(0);
    pad->cd();
    h2->Divide(h1);
    h2->SetTitle(htitle.c_str());
    h2->SetStats(false);
    h2->GetZaxis()->SetTitle("Efficiency");
    h2->Draw("COLZ");
    return;
}

void plot1Dfrom2DforeachXbin(TPad *pad, std::string hlabel)
{
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(hlabel.c_str());
    pad->cd();
    std::vector<double> binsX, binsY;

    // extract bin boundaries
    for (int i=1; i<= h->GetNbinsX()+1; i++)
    {
	binsX.push_back(h->GetXaxis()->GetBinLowEdge(i));
    }

    for (int i=1; i<= h->GetNbinsY()+1; i++)
    {
	binsY.push_back(h->GetYaxis()->GetBinLowEdge(i));
    }

    // Prepare for the legend
    const double lineheight = .025;
    const double legendheight = h->GetNbinsX() * lineheight;
    const double legendwidth = .45;
    const double legendX = 0.4;
    const double legendY = 0.999;
    TLegend *legend = new TLegend(legendX,legendY,legendX+legendwidth,legendY-legendheight);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);

    // Set pad
    pad->SetTopMargin(legendheight+.02);
    pad->SetBottomMargin(0.20);
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0.15);
    pad->SetLogx();

    // now do plots
    for (int i=1; i<= h->GetNbinsX(); i++)
    {
	TH1F *h1 = new TH1F((hlabel+"_"+toString(i)).c_str(),"h",binsY.size()-1,&binsY[0]);
	h1->SetFillColor(0);
	h1->SetStats(false);
	for (int j=1; j<=h->GetNbinsY(); j++)
	{
	    h1->SetBinContent(j,h->GetCellContent(i,j));
	    h1->SetBinError(j,h->GetCellError(i,j));
	    cout << "(" << i << "," << j << ") " << h->GetCellContent(i,j) << " +/- " << h->GetCellError(i,j) << endl;
	}
	h1->SetLineColor(i);
	h1->SetMarkerColor(i);
	h1->GetXaxis()->SetTitle(h->GetYaxis()->GetTitle());
	h1->GetYaxis()->SetTitle(h->GetZaxis()->GetTitle());
	h1->Draw(i>1 ? "samep" : "p");
	const std::string legtit = toString(binsX[i-1]) + " < " + h->GetXaxis()->GetTitle() + " < " + toString(binsX[i]);
	legend->AddEntry(h1,legtit.c_str(),"l");
    }
    legend->Draw();
}

void setTH2params(TPad *pad, TH2F *h, bool logZ = false)
{
    h->GetXaxis()->SetNdivisions(509);
    h->Draw("COLZ");
    pad->SetLogz(logZ);
}

void doRatioPlot(TPad *pad, TTree *t1, TTree *t2, std::string hname, std::string todraw1, std::string todraw2, std::string cut1, std::string cut2,
	int nBins, double min, double max, std::string title, std::string titleX, std::string unitX)
{
    const int labelfont = 43;
    const int labelsize = 20;
    const double yaxistitleoffset = 1.6;
    const std::string plotstring = "(" + toString(nBins) + "," + toString(min) + "," + toString(max) + ")";
    const std::string hname1 = hname + "1";
    const std::string hname2 = hname + "2";
    const std::string plotstring1 = todraw1 + ">>" + hname1 + plotstring;
    const std::string plotstring2 = todraw2 + ">>" + hname2 + plotstring;
    // Set up the first pad
    pad->cd();
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetTopMargin(0.15);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.20);
    pad1->Draw();
    pad1->cd();

    // Draw the first histogram
    t1->Draw(plotstring1.c_str(), cut1.c_str());
    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject(hname1.c_str());
    h1->SetStats(0);
    h1->SetTitle(title.c_str());
    h1->GetYaxis()->SetTitleOffset(yaxistitleoffset);
    h1->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    const double binsize = (max-min)/nBins;
    h1->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    h1->SetLineColor(1);
    h1->SetMinimum(-0.05*h1->GetMaximum());
    // Draw the second histogram
    t2->Draw(plotstring2.c_str(), cut2.c_str());
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject(hname2.c_str());
    h2->SetStats(0);
    h2->GetXaxis()->SetLabelFont(labelfont); //font in pixels
    h2->GetXaxis()->SetLabelSize(labelsize); //in pixels
    h2->GetYaxis()->SetLabelFont(labelfont); //font in pixels
    h2->GetYaxis()->SetLabelSize(labelsize); //in pixels
    h2->GetYaxis()->SetTitleOffset(yaxistitleoffset);
    h2->SetTitle(title.c_str());
    h2->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h2->GetXaxis()->SetTitleFont(labelfont);
    h2->GetXaxis()->SetTitleSize(labelsize);
    h2->GetXaxis()->SetTitleOffset(4);
    h2->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    h2->SetLineColor(4);
    h1->Draw();
    h2->DrawCopy("same");
    // now draw the ratio
    pad->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.20);
    pad2->Draw();
    pad2->cd();
    h2->Sumw2();
    h2->SetStats(0);
    h2->Divide(h1);
    h2->SetMarkerStyle(21);
    h2->SetTitle("");
    h2->GetYaxis()->SetTitle("ratio");
    h2->GetYaxis()->SetTitleFont(labelfont);
    h2->GetYaxis()->SetTitleSize(labelsize);
    h2->GetYaxis()->SetTitleOffset(1.7*yaxistitleoffset);
    h2->GetYaxis()->SetNdivisions(505);
    h2->SetMinimum(0);
    h2->Draw("ep");

    return; 
}

void doRatioPlot(TPad *pad, TTree *t1, TTree *t2, std::string hname, std::string todraw1, std::string todraw2, 
	int nBins, double min, double max, std::string title, std::string titleX, std::string unitX)
{
    doRatioPlot(pad,t1,t2,hname,todraw1,todraw2,"","",nBins,min,max,title,titleX,unitX);
}

void doRatioPlot(TPad *pad, TTree *t1, TTree *t2, std::string hname, std::string todraw1, std::string todraw2, std::string cut1, std::string cut2,
	std::vector<double> bins, std::string title, std::string titleX, std::string unitX)
{
    const int labelfont = 43;
    const int labelsize = 20;
    const double yaxistitleoffset = 1.6;
    const std::string hname1 = hname + "1";
    const std::string hname2 = hname + "2";
    const std::string plotstring1 = todraw1 + ">>" + hname1;
    const std::string plotstring2 = todraw2 + ">>" + hname2;
    // Set up the first pad
    pad->cd();
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetTopMargin(0.15);
    pad1->SetBottomMargin(0);
    pad1->SetLeftMargin(0.20);
    pad1->Draw();
    pad1->cd();

    // Draw the first histogram
    TH1F *h1 = new TH1F(hname1.c_str(),title.c_str(),bins.size()-1,&bins[0]);
    t1->Draw(plotstring1.c_str(), cut1.c_str());
    h1->SetStats(0);
    h1->SetTitle(title.c_str());
    h1->GetYaxis()->SetTitleOffset(yaxistitleoffset);
    h1->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h1->GetYaxis()->SetTitle("entries per bin");
    h1->SetLineColor(1);
    h1->SetMinimum(-0.05*h1->GetMaximum());
    // Draw the second histogram
    TH1F *h2 = new TH1F(hname2.c_str(),title.c_str(),bins.size()-1,&bins[0]);
    t2->Draw(plotstring2.c_str(), cut2.c_str());
    h2->SetStats(0);
    h2->GetXaxis()->SetLabelFont(labelfont); //font in pixels
    h2->GetXaxis()->SetLabelSize(labelsize); //in pixels
    h2->GetYaxis()->SetLabelFont(labelfont); //font in pixels
    h2->GetYaxis()->SetLabelSize(labelsize); //in pixels
    h2->GetYaxis()->SetTitleOffset(yaxistitleoffset);
    h2->SetTitle(title.c_str());
    h2->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h2->GetXaxis()->SetTitleFont(labelfont);
    h2->GetXaxis()->SetTitleSize(labelsize);
    h2->GetXaxis()->SetTitleOffset(4);
    h2->GetYaxis()->SetTitle("entries per bin");
    h2->SetLineColor(4);
    h1->Draw();
    h2->DrawCopy("same");
    // now draw the ratio
    pad->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.20);
    pad2->Draw();
    pad2->cd();
    h2->Sumw2();
    h2->SetStats(0);
    h2->Divide(h1);
    h2->SetMarkerStyle(21);
    h2->SetTitle("");
    h2->GetYaxis()->SetTitle("ratio");
    h2->GetYaxis()->SetTitleFont(labelfont);
    h2->GetYaxis()->SetTitleSize(labelsize);
    h2->GetYaxis()->SetTitleOffset(1.7*yaxistitleoffset);
    h2->GetYaxis()->SetNdivisions(505);
    h2->SetMinimum(0);
    h2->Draw("ep");

    return; 
}

void repositionStatbox(std::string name)
{
    // Reposition and resize statbox
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(name.c_str());
    TPaveStats *st = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.53);
    st->SetY1NDC(0.72);
    st->SetX2NDC(0.99);
    st->SetY2NDC(0.99);
    st->Draw();
}

#endif

