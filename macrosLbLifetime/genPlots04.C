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
#include "cuts.C"

TCanvas *c, *c2;

//template <typename T>
//std::string toString(T i)
//{
//    ostringstream oss;
//    oss << i;
//    return oss.str();
//}

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
    st->SetX1NDC(0.63);
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

void genPlots04(std::string fullPath)
{
    const std::string outpdf("genPlots04plots");
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c1","c1",1000,600);
    int canvasCdCounter(0), canvasPageCounter(0);
    const unsigned int nPadX = 2;
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
    TTree* t = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;

    // General cuts
    Cuts cuts;
    cuts.selectCut("an03","acc01");
    std::string cutsAnalysis = cuts.getCut();

    //std::string cutsgen = "TMath::Abs(etamu1)<2.5";
    //std::string cutsgen = "TMath::Abs(etamu1)<2.5&&TMath::Abs(etamu2)<2.5";
    std::string cutsgen = cutsAnalysis + "&&isSig==1&&isMCmatch==1";
    // Do plots
    const double ptmax(1.0);
    const double etamax(.03);
    const double phimax(.03);
    const double ymax(.03);
    const int nBins(40);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hpt", "ptlb-ptgenlb", cutsgen,nBins,-ptmax,ptmax,"#Deltap_{T}(#Lambda_{b})","p_{T}(reco)-p_{T}(gen)","GeV/c");
    c->Update();
    repositionStatbox("hpt");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"heta", "etalb-etagenlb", cutsgen,nBins,-etamax,etamax,"#Delta#eta(#Lambda_{b})","#eta(reco)-#eta(gen)","");
    c->Update();
    repositionStatbox("heta");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hphi", "philb-phigenlb", cutsgen,nBins,-phimax,phimax,"#Delta#phi(#Lambda_{b})","#phi(reco)-#phi(gen)","");
    c->Update();
    repositionStatbox("hphi");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hy", "ylb-ygenlb", cutsgen,nBins,-ymax,ymax,"#Deltay(#Lambda_{b})","y(reco)-y(gen)","");
    c->Update();
    repositionStatbox("hy");

    // finalize current page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());
}

