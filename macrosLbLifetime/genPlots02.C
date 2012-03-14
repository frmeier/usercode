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
#include "TMath.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TEventList.h"
#include "setTDRStyle_modified.C"

TCanvas *c, *c2;

template <typename T>
std::string toString(T i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}

void doPlot1d(TTree *t, std::string hname, std::string todraw, std::string cut,
	int nBins, double min, double max, std::string title, std::string titleX, std::string unitX)
{
    const std::string plotstring = todraw + ">>" + hname + "(" +
	toString(nBins) + "," + toString(min) + "," + toString(max) + ")";
    t->Draw(plotstring.c_str(), cut.c_str());
    TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(hname.c_str());
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    const double binsize = (max-min)/nBins;
    h->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    return; 
}

void doPlot1d(TTree *t, std::string hname, std::string todraw, 
	int nBins, double min, double max, std::string title, std::string titleX, std::string unitX)
{
    doPlot1d(t,hname,todraw,"",nBins,min,max,title,titleX,unitX);
}

double doPlot1dIntegral(TH1F *h, std::string hname, std::string title, std::string titleX, std::string unitX, double scale=0)
{
    const double min = h->GetXaxis()->GetXmin();
    const double max = h->GetXaxis()->GetXmax();
    const int nBins = h->GetNbinsX();
    TH1F *hi = new TH1F(hname.c_str(), title.c_str(), nBins, min, max);
    double sum = 0;
    for (int i=1; i!=h->GetNbinsX(); i++)
    {
	sum+=h->GetBinContent(i);
	hi->SetBinContent(i,sum);
    }
    hi->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    const double binsize = (max-min)/nBins;
    hi->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    if(scale==0)
	hi->Scale(1./sum);
    else
	hi->Scale(scale);
    hi->Draw();
    return sum;
}

void doPlot2dSym(TTree *t, std::string hname, std::string todraw, 
	int nBins, double min, double max, std::string title, std::string titleX, std::string titleY, std::string unit)
{
    const std::string plotstring = todraw + ">>" + hname + "("
	+ toString(nBins) + "," + toString(min) + "," + toString(max) + ","
	+ toString(nBins) + "," + toString(min) + "," + toString(max) + ")";
    cout << plotstring << endl;
    t->Draw(plotstring.c_str(),"","COLZ");
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(hname.c_str());
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unit.size()>0 ? (titleX+" / "+unit).c_str() : titleX.c_str());
    h->GetYaxis()->SetTitle(unit.size()>0 ? (titleY+" / "+unit).c_str() : titleY.c_str());
    //const double binsize = (max-min)/nBins;
    //h->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    return; 
}

void doPlot2dSym(TTree *t, std::string hname, std::string todraw, 
	std::vector<double> bins, std::string title, std::string titleX, std::string titleY, std::string unit)
{
    TH2F *h = new TH2F(hname.c_str(),title.c_str(),bins.size()-1,&bins[0],bins.size()-1,&bins[0]);
    const std::string plotstring = todraw + ">>" + hname;
    cout << plotstring << endl;
    t->Draw(plotstring.c_str(),"","COLZTEXT");
    h->GetXaxis()->SetTitle(unit.size()>0 ? (titleX+" / "+unit).c_str() : titleX.c_str());
    h->GetYaxis()->SetTitle(unit.size()>0 ? (titleY+" / "+unit).c_str() : titleY.c_str());
    return; 
}

void doPlot2d(TTree *t, std::string hname, std::string todraw, 
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

void doPlot2d(TTree *t, std::string hname, std::string todraw, 
	std::vector<double> binsX, std::vector<double> binsY,
        std::string title, std::string titleX, std::string titleY,
	std::string unitX, std::string unitY)
{
    TH2F *h = new TH2F(hname.c_str(),title.c_str(),binsX.size()-1,&binsX[0],binsY.size()-1,&binsY[0]);
    const std::string plotstring = todraw + ">>" + hname;
    cout << plotstring << endl;
    t->Draw(plotstring.c_str(),"","COLZTEXT");
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h->GetYaxis()->SetTitle(unitY.size()>0 ? (titleY+" / "+unitY).c_str() : titleY.c_str());
    //const double binsize = (max-min)/nBins;
    //h->GetYaxis()->SetTitle(("entries per "+toString(binsize)+" "+unitX).c_str());
    return; 
}

void repositionPalette(std::string name)
{
    // Reposition and resize palette
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(name.c_str());
    TPaletteAxis *pal;
    pal = (TPaletteAxis*)h->FindObject("palette");
    if (0==pal)
    {
	cout << "Palette not found - cannot reposition it" << endl;
	return;
    }
    pal->SetX1NDC(0.83);
    pal->SetY1NDC(0.14);
    pal->SetX2NDC(0.88);
    pal->SetY2NDC(0.60);
    pal->SetLabelSize(.04);
    cout << "Repo" << endl;
}

void drawTracker(TPad* pad)
{
	double r = 0, z =0;
	TLine *l;
	// PXB
	r = 4.4; l = new TLine(0,r,26.5,r); l->Draw();
	r = 7.3; l = new TLine(0,r,26.5,r); l->Draw();
	r = 10.2; l = new TLine(0,r,26.5,r); l->Draw();
	// PXE
	z = 34.5; l = new TLine(z,6,z,15); l->Draw();
	z = 46.5; l = new TLine(z,6,z,15); l->Draw();
	// TIB
	r = 25.5; l = new TLine(0,r,70,r); l->Draw();
	r = 33.9; l = new TLine(0,r,70,r); l->Draw();
	r = 41.9; l = new TLine(0,r,70,r); l->Draw();
	r = 49.8; l = new TLine(0,r,70,r); l->Draw();
	// TID (z-Positionen geraten aus Bild)
	z = 80; l = new TLine(z,20,z,50); l->Draw();
	z = 90; l = new TLine(z,20,z,50); l->Draw();
	z = 100; l = new TLine(z,20,z,50); l->Draw();
	// TOB
	r = 60.8; l = new TLine(0,r,109,r); l->Draw();
	r = 69.2; l = new TLine(0,r,109,r); l->Draw();
	r = 78.0; l = new TLine(0,r,109,r); l->Draw();
	r = 86.8; l = new TLine(0,r,109,r); l->Draw();
	r = 96.5; l = new TLine(0,r,109,r); l->Draw();
	r = 108; l = new TLine(0,r,109,r); l->Draw();
	// TOB (wiederum alles aus Zeichnung geschaetzt)
	z = 124; l = new TLine(z,22,z,113.5); l->Draw();
	z = 140; l = new TLine(z,22,z,113.5); l->Draw();
	z = 155; l = new TLine(z,22,z,113.5); l->Draw();
	z = 168; l = new TLine(z,33,z,113.5); l->Draw();
	z = 186; l = new TLine(z,33,z,113.5); l->Draw();
	z = 203; l = new TLine(z,33,z,113.5); l->Draw();
	z = 223; l = new TLine(z,40,z,113.5); l->Draw();
	z = 247; l = new TLine(z,40,z,113.5); l->Draw();
	z = 272; l = new TLine(z,52,z,113.5); l->Draw();
}

// draws the decay vertices of the Lambda_s
void genPlots01(std::string fullPath, int nOverlay = 500, bool custBinning = false)
{
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c1","c1",1000,600);
    const unsigned int nPadX = 1;
    const unsigned int nPadY = 1;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    for(unsigned int i=1; i<=nPads; i++)
    {
	TPad* pad= (TPad*)c->cd(i);
	pad->SetTopMargin(0.10);
	pad->SetRightMargin(0.20);
	pad->SetLeftMargin(0.15);
    }
    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTree
    TTree* t = (TTree*) f->Get("genevents");
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;
    // Do a cut, if needed
    //t->Draw(">>lst","chi2lb>.1&&mlb>5.61&&mlb<5.63");
    //t->Draw(">>lst","nRef2G==511");
    t->Draw(">>lst","ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<2.5&&TMath::Abs(etamu2)<2.5");
    TEventList *lst;
    lst = (TEventList*)gDirectory->Get("lst");
    t->SetEventList(lst);

    // Do plots
    c->cd(1);
    //doPlot2d(t,"hrzL0vtx", "vrl0:TMath::Abs(vzl0)",30,0,300,30,0,120,"Tit","|z|","r","cm","cm");
    if (custBinning)
    {
	double newbinsX[]={0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50};
	const int newbinsX_size = sizeof(newbinsX)/sizeof(double);
	std::vector<double> binvecX(newbinsX,newbinsX+newbinsX_size);
	//double newbinsY[]={0,1,2,3,4,5,6,7,8,9,10};
	double newbinsY[]={0,0.5,1,2,4,8,16,32};
	const int newbinsY_size = sizeof(newbinsY)/sizeof(double);
	std::vector<double> binvecY(newbinsY,newbinsY+newbinsY_size);
	doPlot2d(t,"hrzL0vtx", "vrl0:TMath::Abs(vzl0)",binvecX, binvecY,"#Lambda vertices","|z|","r","cm","cm");
    }
    else
    {
	doPlot2d(t,"hrzL0vtx", "vrl0:TMath::Abs(vzl0)",30,0,50,30,0,30,"#Lambda vertices","|z|","r","cm","cm");
    }

    // add tracker
    TPad* pad;
    pad = (TPad*)c->cd(1);
    pad->Modified();
    pad->Update();
    repositionPalette("hrzL0vtx");
    pad->Update();
    pad->SetLogz();
    drawTracker(pad);

    if (nOverlay<=0) return;
    int maxN = nOverlay;
    if (maxN > t->GetEntries()) maxN = t->GetEntries();
    double vrl0,vzl0,ppr,ppi,etapr,etapi;
    t->SetBranchAddress("vrl0",&vrl0);
    t->SetBranchAddress("vzl0",&vzl0);
    t->SetBranchAddress("ppr",&ppr);
    t->SetBranchAddress("etapr",&etapr);
    t->SetBranchAddress("ppi",&ppi);
    t->SetBranchAddress("etapi",&etapi);
    double scalepr = 4;
    double scalepi = 8;

    { // reference indicator
	const double x1pr = 0; const double y1pr = -3;
	const double x2pr = scalepr; const double y2pr = y1pr;
	const double versatz = 14;
	const double x1pi = x1pr+versatz; const double y1pi = -3;
	const double x2pi = x2pr+versatz+scalepi; const double y2pi = y1pi;
	TArrow *a;
	a = new TArrow(x1pr,y1pr,x2pr,y2pr,.01,">");
	a->SetLineColor(24);
	a->Draw();
	TLatex tl;
	tl.SetTextSize(20);
	tl.SetTextFont(4);
	tl.DrawLatex(x1pr,y2pr-1.2,"p(p) / 1 GeV");
	a = new TArrow(x1pi,y1pi,x2pi,y2pi,.01,">");
	a->SetLineColor(20);
	a->Draw();
	tl.SetTextSize(20);
	tl.SetTextFont(4);
	tl.DrawLatex(x1pi,y2pi-1.2,"p(#pi) / 1 GeV");
    }
    for (int i = 0; i!=maxN; i++)
    {
	t->GetEntry(i);
	const double thetapr = 2*TMath::ATan(TMath::Exp(-TMath::Abs(etapr)));
	const double thetapi = 2*TMath::ATan(TMath::Exp(-TMath::Abs(etapi)));
	const double x1=TMath::Abs(vzl0);
	const double y1=vrl0;
	const double x2pr=x1+scalepr*ppr*TMath::Cos(thetapr);
	const double y2pr=y1+scalepr*ppr*TMath::Sin(thetapr);
	const double x2pi=x1+scalepi*ppi*TMath::Cos(thetapi);
	const double y2pi=y1+scalepi*ppi*TMath::Sin(thetapi);
	TArrow *a;
        a = new TArrow(x1,y1,x2pr,y2pr,.01,">");
	a->SetLineColor(24);
	a->Draw();
        a = new TArrow(x1,y1,x2pi,y2pi,.01,">");
	a->SetLineColor(20);
	a->Draw();
	TMarker *m = new TMarker(x1,y1,7);
	m->SetMarkerColor(28);
	m->Draw();
    }
}

// draw the same thing but after reco
void genPlots02(std::string fullPath, int nOverlay = 500, bool custBinning = false)
{
    const int fVerbose(1);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c2","c2",1000,600);
    const unsigned int nPadX = 1;
    const unsigned int nPadY = 1;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    for(unsigned int i=1; i<=nPads; i++)
    {
	TPad* pad= (TPad*)c->cd(i);
	pad->SetTopMargin(0.10);
	pad->SetRightMargin(0.20);
	pad->SetLeftMargin(0.15);
    }
    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTree
    TTree* t = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;
    // Do a cut, if needed
    //t->Draw(">>lst","chi2lb>.1&&mlb>5.61&&mlb<5.63");
    //t->Draw(">>lst","chi2lb>.1&&isSig==1");
    t->Draw(">>lst","(rid1m&4)==4&&(rid2m&4)==4&&mjp>2.895&&mjp<3.295&&prob1m>0.1&&prob2m>0.1&&ptjp>2&&probjp>0.005&&ml0>1.101&&ml0<1.129&&probpr>0.02&&probpi>0.02&&rptpr>rptpi&&ptl0>3&&rptpr>1&&rptpi>0.5&&probl0>0.02&&alphal0<0.3&&d3l0>1&&d3l0/d3El0>10&&problb>0.001&&alphalb<0.3");
    TEventList *lst;
    lst = (TEventList*)gDirectory->Get("lst");
    t->SetEventList(lst);
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;

    // Do plots
    c->cd(1);
    //doPlot2d(t,"hrzL0vtx", "vrl0:TMath::Abs(vzl0)",30,0,300,30,0,120,"Tit","|z|","r","cm","cm");
    if (custBinning)
    {
	double newbinsX[]={0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50};
	const int newbinsX_size = sizeof(newbinsX)/sizeof(double);
	std::vector<double> binvecX(newbinsX,newbinsX+newbinsX_size);
	//double newbinsY[]={0,1,2,3,4,5,6,7,8,9,10};
	double newbinsY[]={0,0.5,1,2,4,8,16,32};
	const int newbinsY_size = sizeof(newbinsY)/sizeof(double);
	std::vector<double> binvecY(newbinsY,newbinsY+newbinsY_size);
	doPlot2d(t,"hrzL0vtxreco", "vrl0:TMath::Abs(vzl0)",binvecX, binvecY,"#Lambda vertices","|z|","r","cm","cm");
    }
    else
    {
	doPlot2d(t,"hrzL0vtxreco", "vrl0:TMath::Abs(vzl0)",30,0,50,30,0,30,"#Lambda vertices","|z|","r","cm","cm");
    }

    // add tracker
    TPad* pad;
    pad = (TPad*)c->cd(1);
    pad->Modified();
    pad->Update();
    repositionPalette("hrzL0vtxreco");
    pad->Update();
    pad->SetLogz();
    drawTracker(pad);

    if (nOverlay<=0) return;
    int maxN = nOverlay;
    if (maxN > t->GetEntries()) maxN = t->GetEntries();
    double vrl0,vzl0,ppr,ppi,etapr,etapi;
    t->SetBranchAddress("vrl0",&vrl0);
    t->SetBranchAddress("vzl0",&vzl0);
    t->SetBranchAddress("ppr",&ppr);
    t->SetBranchAddress("etapr",&etapr);
    t->SetBranchAddress("ppi",&ppi);
    t->SetBranchAddress("etapi",&etapi);
    double scalepr = 4;
    double scalepi = 8;

    { // reference indicator
	const double x1pr = 0; const double y1pr = -3;
	const double x2pr = scalepr; const double y2pr = y1pr;
	const double versatz = 14;
	const double x1pi = x1pr+versatz; const double y1pi = -3;
	const double x2pi = x2pr+versatz+scalepi; const double y2pi = y1pi;
	TArrow *a;
	a = new TArrow(x1pr,y1pr,x2pr,y2pr,.01,">");
	a->SetLineColor(24);
	a->Draw();
	TLatex tl;
	tl.SetTextSize(20);
	tl.SetTextFont(4);
	tl.DrawLatex(x1pr,y2pr-1.2,"p(p) / 1 GeV");
	a = new TArrow(x1pi,y1pi,x2pi,y2pi,.01,">");
	a->SetLineColor(20);
	a->Draw();
	tl.SetTextSize(20);
	tl.SetTextFont(4);
	tl.DrawLatex(x1pi,y2pi-1.2,"p(#pi) / 1 GeV");
    }
    for (int i = 0; i!=maxN; i++)
    {
	t->GetEntry(i);
	const double thetapr = 2*TMath::ATan(TMath::Exp(-TMath::Abs(etapr)));
	const double thetapi = 2*TMath::ATan(TMath::Exp(-TMath::Abs(etapi)));
	const double x1=TMath::Abs(vzl0);
	const double y1=vrl0;
	const double x2pr=x1+scalepr*ppr*TMath::Cos(thetapr);
	const double y2pr=y1+scalepr*ppr*TMath::Sin(thetapr);
	const double x2pi=x1+scalepi*ppi*TMath::Cos(thetapi);
	const double y2pi=y1+scalepi*ppi*TMath::Sin(thetapi);
	TArrow *a;
        a = new TArrow(x1,y1,x2pr,y2pr,.01,">");
	a->SetLineColor(24);
	a->Draw();
        a = new TArrow(x1,y1,x2pi,y2pi,.01,">");
	a->SetLineColor(20);
	a->Draw();
	TMarker *m = new TMarker(x1,y1,7);
	m->SetMarkerColor(28);
	m->Draw();
    }
}

// draw r distributiob
void genPlots03(std::string fullPath, int nOverlay = 500)
{
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c3","c3",1000,600);
    const unsigned int nPadX = 1;
    const unsigned int nPadY = 1;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    for(unsigned int i=1; i<=nPads; i++)
    {
	TPad* pad= (TPad*)c->cd(i);
	pad->SetTopMargin(0.10);
	pad->SetRightMargin(0.20);
	pad->SetLeftMargin(0.15);
	pad->SetLogy();
    }
    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTree
    TTree* t = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;
    // Do plots
    c->cd(1);
    //doPlot2d(t,"hrzL0vtx", "vrl0:TMath::Abs(vzl0)",30,0,300,30,0,120,"Tit","|z|","r","cm","cm");
    std::string histoName = "hrL0vtxreco";
    doPlot1d(t,histoName.c_str(), "vrl0",30,0,30,"#Lambda vertices","r","cm");
    // draw radii
    {
	TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(histoName.c_str());
	double r;
	const double min = h->GetMinimum();
	const double max = h->GetMaximum();

	TLine *l;
	r = 4.4; l = new TLine(r,min,r,max); l->Draw();
	r = 7.3; l = new TLine(r,min,r,max); l->Draw();
	r = 10.2; l = new TLine(r,min,r,max); l->Draw();
	r = 25.5; l = new TLine(r,min,r,max); l->Draw();
    }
}

void genPlotsRatio0102(std::string fullPath, bool custBinning = false)
{
    genPlots01(fullPath, 0, custBinning);
    TH2F *h1 = (TH2F*)gDirectory->GetList()->FindObject("hrzL0vtx");
    genPlots02(fullPath, 0, custBinning);
    TH2F *h2 = (TH2F*)gDirectory->GetList()->FindObject("hrzL0vtxreco");
    if (0==h1)
    {
	cout << "Histo hrzL0vtx not present" << endl;
	return;
    }
    if (0==h2)
    {
	cout << "Histo hrzL0vtxreco not present" << endl;
	return;
    }
    TH2F *hratio;
    if (custBinning)
    {
	double newbinsX[]={0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50};
	const int newbinsX_size = sizeof(newbinsX)/sizeof(double);
	double newbinsY[]={0,0.5,1,2,4,8,16,32};
	const int newbinsY_size = sizeof(newbinsY)/sizeof(double);
	hratio = new TH2F("hratio","hratio",newbinsX_size-1,newbinsX,newbinsY_size-1,newbinsY); 
    }
    else
    {
	hratio = new TH2F("hratio","hratio",30,0,50,30,0,30); 
    }
    c = new TCanvas("c3","c3",1000,600);
    const unsigned int nPadX = 1;
    const unsigned int nPadY = 1;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    for(unsigned int i=1; i<=nPads; i++)
    {
	TPad* pad= (TPad*)c->cd(i);
	pad->SetTopMargin(0.10);
	pad->SetRightMargin(0.20);
	pad->SetLeftMargin(0.15);
	//pad->SetLogy();
    }
    for(int i=1; i<=h1->GetNbinsX(); i++)
	for(int j=1; j<=h1->GetNbinsY(); j++)
	{
	    if (h1->GetBinContent(i,j) != 0.)
		hratio->SetBinContent(i,j,h2->GetBinContent(i,j) / h1->GetBinContent(i,j));
	    else
		hratio->SetBinContent(i,j,0);
	}
    hratio->Draw("COLZ");

    // add tracker
    TPad* pad;
    pad = (TPad*)c->cd(1);
    pad->Modified();
    pad->Update();
    repositionPalette("hratio");
    pad->Update();

    drawTracker(pad);
}
