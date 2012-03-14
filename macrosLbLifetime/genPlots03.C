#include <iostream>
#include <string>
#include <sstream>
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
    pal->SetX1NDC(0.83);
    pal->SetY1NDC(0.14);
    pal->SetX2NDC(0.88);
    pal->SetY2NDC(0.60);
    pal->SetLabelSize(.04);
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

void genPlots01(std::string fullPath)
{
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c1","c1",1000,600);
    const unsigned int nPadX = 4;
    const unsigned int nPadY = 2;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    for(unsigned int i=1; i<=nPads; i++)
    {
	TPad* pad= (TPad*)c->cd(i);
	pad->SetTopMargin(0.10);
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
    // General cuts
    //std::string cutsgen = "TMath::Abs(etamu1)<2.5";
    //std::string cutsgen = "TMath::Abs(etamu1)<2.5&&TMath::Abs(etamu2)<2.5";
    std::string cutsgen = "TMath::Abs(etamu1)<2.5&&TMath::Abs(etamu2)<2.5&&TMath::Abs(etapr)<2.5&&TMath::Abs(etapi)<2.5&&ptmu1>3&&ptmu2>3";
    // Do plots
    c->cd(1);
    doPlot1d(t,"hptmu1", "ptmu1", cutsgen,40,0,20,"p_{T} #mu_{1}","p_{T} #mu_{1}","GeV/c");
    c->Update();
    repositionStatbox("hptmu1");

    c->cd(1+nPadX);
    doPlot1d(t,"hetamu1", "etamu1", cutsgen, 40,0,4,"#eta #mu_{1}","#eta #mu_{1}","");
    c->Update();
    repositionStatbox("hetamu1");

    c->cd(2);
    doPlot1d(t,"hptmu2", "ptmu2", cutsgen, 40,0,20,"p_{T} #mu_{2}","p_{T} #mu_{2}","GeV/c");
    c->Update();
    repositionStatbox("hptmu2");
    
    c->cd(2+nPadX);
    doPlot1d(t,"hetamu2", "etamu2", cutsgen, 40,0,4,"#eta #mu_{2}","#eta #mu_{2}","");
    c->Update();
    repositionStatbox("hetamu2");
    
    c->cd(3);
    doPlot1d(t,"hptpr", "ptpr", cutsgen, 40,0,10,"p_{T} p","p_{T} p","GeV/c");
    c->Update();
    repositionStatbox("hptpr");
    
    c->cd(3+nPadX);
    doPlot1d(t,"hetapr", "etapr", cutsgen, 40,0,4,"#eta p","#eta p","");
    c->Update();
    repositionStatbox("hetapr");
    
    c->cd(4);
    doPlot1d(t,"hptpi", "ptpi", cutsgen, 40,0,5,"p_{T} #pi","p_{T} #pi","GeV/c");
    c->Update();
    repositionStatbox("hptpi");
    
    c->cd(4+nPadX);
    doPlot1d(t,"hetapi", "etapi", cutsgen, 40,0,4,"#eta #pi","#eta #pi","");
    c->Update();
    repositionStatbox("hetapi");
}

