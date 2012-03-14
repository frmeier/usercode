#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaletteAxis.h"
#include "TPaveStats.h"
#include "setTDRStyle_modified.C"
#include "cuts.C"

TCanvas *c, *c2;

/*
template <typename T>
std::string toString(T i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}
*/

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
    h2->Draw("ep");

    return; 
}

void doRatioPlot(TPad *pad, TTree *t1, TTree *t2, std::string hname, std::string todraw1, std::string todraw2, 
	int nBins, double min, double max, std::string title, std::string titleX, std::string unitX)
{
    doRatioPlot(pad,t1,t2,hname,todraw1,todraw2,"","",nBins,min,max,title,titleX,unitX);
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

void canvaspager(TCanvas *c, std::string filename, const int &ncd, int &cdcounter, int &pagecounter)
{
    cdcounter++;
    if (cdcounter > ncd)
    {
	c->SaveAs((filename + toString(pagecounter) + ".pdf").c_str());
	pagecounter++;
	cdcounter=1;
	c->Clear("D");
    }
}

void acceptance01(std::string fullPath)
{
    const std::string outpdf("acceptance01plot");
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    //c = new TCanvas("c1","c1",1000,600);
    c = new TCanvas("c1","c1",1400,1200);
    //c = new TCanvas("c1","c1",500,600);
    const unsigned int nPadX = 2;
    const unsigned int nPadY = 2;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTrees
    TTree* tgen = (TTree*) f->Get("genevents");
    if (tgen==0)
    {
	cout << "Tree genevents not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0) cout << "Got TTree genevents with " << tgen->GetEntries() << " entries" << endl;
    TTree* treco = (TTree*) f->Get("events");
    if (treco==0)
    {
	cout << "Tree events not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0) cout << "Got TTree evemts with " << treco->GetEntries() << " entries" << endl;

    // General cuts
    Cuts cuts;
    cuts.selectCut("an01","acc01","L1_0_01");
    //cuts.selectCut("an01","acc01","HLT_01");
    //cuts.selectCut("an01","acc01","L1_0_01","HLT_01");
    std::string cutsAnalysis = cuts.getCut();
    //std::string cutsgen = "TMath::Abs(etamu1)<2.5&&TMath::Abs(etamu2)<2.5&&TMath::Abs(etapr)<2.5&&TMath::Abs(etapi)<2.5&&ptmu1>3&&ptmu2>3";
    //std::string cutsreco = "isSig==1&&isMCmatch==1&&TMath::Abs(Seta1m)<2.5&&TMath::Abs(Seta2m)<2.5&&TMath::Abs(Setapr)<2.5&&TMath::Abs(Setapi)<2.5&&Spt1m>3&&Spt2m>3";
    //std::string cutsgen = "ptpi>.6&&ptpr>2";
    std::string cutsgen = "((TMath::Abs(etamu1)<1.3&&ptmu1>3.3)||(TMath::Abs(etamu1)>=1.3&&TMath::Abs(etamu1)<2.2&&pmu1>2.9)||(TMath::Abs(etamu1)>=2.2&&TMath::Abs(etamu1)<2.4&&ptmu1>0.8))";
    cutsgen += "&&((TMath::Abs(etamu2)<1.3&&ptmu2>3.3)||(TMath::Abs(etamu2)>=1.3&&TMath::Abs(etamu2)<2.2&&pmu2>2.9)||(TMath::Abs(etamu2)>=2.2&&TMath::Abs(etamu2)<2.4&&ptmu2>0.8))";
    //cutsgen += "&&vrl0>1&&vrl0<35&&TMath::Abs(vzl0)<100";
    cutsgen += "&&ptpi>.6&&ptpr>2";
    //std::string cutsreco = "isSig==1&&isMCmatch==1";
    std::string cutsreco = "isSig==1&&isMCmatch==1&&" + cutsAnalysis;

    // Do plots
    int canvasCdCounter(0), canvasPageCounter(0);
    const int nbinseta(13);
    const double etamax(2.6);
    const int nbinspt(20);
    //const int nbinseta(2);
    //const double etamax(2.4);
    //const int nbinspt(1);
    const bool etabetrag(true);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    if(etabetrag)
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetamu1", "TMath::Abs(etamu1)", "TMath::Abs(Seta1m)", cutsgen,cutsreco,nbinseta,0,etamax,"#eta(#mu_{1})","|#eta(#mu_{1})|","");
    else
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetamu1", "etamu1", "Seta1m", cutsgen,cutsreco,nbinseta,-etamax,etamax,"#eta(#mu_{1})","#eta(#mu_{1})","");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hptmu1", "ptmu1", "Spt1m", cutsgen,cutsreco,nbinspt,0,20,"p_{T}(#mu_{1})","p_{T}(#mu_{1})","GeV/c");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    if(etabetrag)
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetamu2", "TMath::Abs(etamu2)", "TMath::Abs(Seta2m)", cutsgen,cutsreco,nbinseta,0,etamax,"#eta(#mu_{2})","|#eta(#mu_{2})|","");
    else
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetamu2", "etamu2", "Seta2m", cutsgen,cutsreco,nbinseta,-etamax,etamax,"#eta(#mu_{2})","|#eta(#mu_{2}|","");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hptmu2", "ptmu2", "Spt2m", cutsgen,cutsreco,nbinspt,0,20,"p_{T}(#mu_{2})","p_{T}(#mu_{2})","GeV/c");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    if(etabetrag)
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetapr", "TMath::Abs(etapr)", "TMath::Abs(Setapr)", cutsgen,cutsreco,nbinseta,0,etamax,"#eta(p)","|#eta(p)|","");
    else
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetapr", "etapr", "Setapr", cutsgen,cutsreco,nbinseta,-etamax,etamax,"#eta(p)","|#eta(p|","");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hptpr", "ptpr", "Sptpr", cutsgen,cutsreco,nbinspt,0,20,"p_{T}(p)","p_{T}(p)","GeV/c");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    if(etabetrag)
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetapi", "TMath::Abs(etapi)", "TMath::Abs(Setapi)", cutsgen,cutsreco,nbinseta,0,etamax,"#eta(#pi)","|#eta(#pi)|","");
    else
	doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetapi", "etapi", "Setapi", cutsgen,cutsreco,nbinseta,-etamax,etamax,"#eta(#pi)","|#eta(#pi|","");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hptpi", "ptpi", "Sptpi", cutsgen,cutsreco,nbinspt,0,6,"p_{T}(#pi)","p_{T}(#pi)","GeV/c");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hvrl0", "vrl0", "vrl0", cutsgen,cutsreco,20,0,50,"r(#Lambda)","r(#Lambda)","cm");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hvzl0", "TMath::Abs(vzl0)", "TMath::Abs(vzl0)", cutsgen,cutsreco,14,0,140,"z(#Lambda)","|z|(#Lambda)","cm");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hetalb", "etalb", "etalb", cutsgen,cutsreco,nbinseta,-etamax,etamax,"#eta(#Lambda_{b})","#eta(#Lambda_{b})","");
    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    doRatioPlot((TPad*)c->cd(canvasCdCounter),tgen,treco,"hptlb", "ptlb", "ptlb", cutsgen,cutsreco,nbinspt,0,20,"p_{T}(#Lambda_{b})","p_{T}(#Lambda_{b})","GeV/c");

    // finalize current page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
}

