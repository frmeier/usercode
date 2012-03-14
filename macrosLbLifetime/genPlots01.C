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

void genPlots01(std::string fullPath)
{
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c1","c1",1000,600);
    const unsigned int nPadX = 3;
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
    // Do plots
    c->cd(1);
    doPlot1d(t,"hanprpi", "anprpi", 20,0,1,"Opening angle p#pi","Opening angle","rad");
    c->cd(2);
    doPlot1d(t,"hanmumu", "anmumu", 20,0,1,"Opening angle #mu#mu","Opening angle","rad");
    c->cd(3);
    doPlot1d(t,"hanl0jp", "anl0lb", 20,0,1,"Opening angle #Lambda^{0}J/#psi","Opening angle","rad");
    c->cd(4);
    doPlot1d(t,"hanl0mumin", "anl0mumin", 20,0,1,"Opening angle #Lambda^{0}#mu_{closest}","Opening angle","rad");
    c->cd(5);
    doPlot1d(t,"hanl0muPt", "anl0muPt", 20,0,1,"Opening angle #Lambda^{0}J/#mu_{highpT}","Opening angle","rad");
    c->cd(6);
    doPlot2dSym(t,"hptmumu", "ptmu1:ptmu2", 20,0,10,"p_{#perp} of #mu#mu","p_{#perp}","p_{#perp}","GeV/c");
}

// make plots of signal yield as fct of eta and pt
// the rebinned versions cover detector acceptances with special trigger features
void genPlots02(std::string fullPath)
{
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c1","c1",1000,600);
    const unsigned int nPadX = 3;
    const unsigned int nPadY = 2;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    for(unsigned int i=1; i<=nPads; i++)
    {
	TPad* pad= (TPad*)c->cd(i);
	pad->SetTopMargin(0.10);
	pad->SetRightMargin(0.20);
	pad->SetLeftMargin(0.20);
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
    // Do plots
    TPad* pad;
    // -- Plots eta/eta
    pad = (TPad*)c->cd(1);
    doPlot2dSym(t,"hetamumu", "TMath::Abs(etamu1):TMath::Abs(etamu2)", 25,0,2.5,"|#eta| of #mu#mu","|#eta| #mu_{1}","|#eta| #mu_{2}","");
    pad->Modified();
    pad->Update();
    repositionPalette("hetamumu");
    pad = (TPad*)c->cd(1+nPadX);
    {
	double newbins[]={0,1.6,2.1,2.5};
	const int newbins_size = sizeof(newbins)/sizeof(double);
	std::vector<double> binvec(newbins,newbins+newbins_size);
	doPlot2dSym(t,"hetamumuRebin", "TMath::Abs(etamu1):TMath::Abs(etamu2)", binvec,"#eta of #mu#mu","|#eta| #mu_{1}","|#eta| #mu_{2}","");
	pad->Modified();
	pad->Update();
	repositionPalette("hetamumuRebin");
    }
    // -- Plots pt/eta 1
    pad = (TPad*)c->cd(2);
    doPlot2d(t,"hptetamu1", "ptmu1:TMath::Abs(etamu1)", 25,0,2.5,10,0,10,"#eta vs p_{#perp} of #mu_{1}","|#eta|","p_{#perp}","","GeV/c");
    pad->Modified();
    pad->Update();
    repositionPalette("hptetamu1");
    pad = (TPad*)c->cd(2+nPadX);
    {
	double newbinsX[]={0,1.6,2.1,2.5};
	const int newbinsX_size = sizeof(newbinsX)/sizeof(double);
	std::vector<double> binvecX(newbinsX,newbinsX+newbinsX_size);
	double newbinsY[]={0,1,2,3,4,5,6,7,8,9,10};
	const int newbinsY_size = sizeof(newbinsY)/sizeof(double);
	std::vector<double> binvecY(newbinsY,newbinsY+newbinsY_size);
	doPlot2d(t,"hetamu1Rebin", "ptmu1:TMath::Abs(etamu1)", binvecX, binvecY,"#eta vs p_{#perp} of #mu_{1}","|#eta|","p_{#perp}","","GeV/c");
	pad->Modified();
	pad->Update();
	repositionPalette("hetamu1Rebin");
    }
    // -- Plots pt/eta 2
    pad = (TPad*)c->cd(3);
    doPlot2d(t,"hptetamu2", "ptmu2:TMath::Abs(etamu2)", 25,0,2.5,10,0,10,"#eta vs p_{#perp} of #mu_{2}","|#eta|","p_{#perp}","","GeV/c");
    pad->Modified();
    pad->Update();
    repositionPalette("hptetamu2");
    pad = (TPad*)c->cd(3+nPadX);
    {
	double newbinsX[]={0,1.6,2.1,2.5};
	const int newbinsX_size = sizeof(newbinsX)/sizeof(double);
	std::vector<double> binvecX(newbinsX,newbinsX+newbinsX_size);
	double newbinsY[]={0,1,2,3,4,5,6,7,8,9,10};
	const int newbinsY_size = sizeof(newbinsY)/sizeof(double);
	std::vector<double> binvecY(newbinsY,newbinsY+newbinsY_size);
	doPlot2d(t,"hetamu2Rebin", "ptmu2:TMath::Abs(etamu2)", binvecX, binvecY,"#eta vs p_{#perp} of #mu_{2}","|#eta|","p_{#perp}","","GeV/c");
	pad->Modified();
	pad->Update();
	repositionPalette("hetamu2Rebin");
    }
}

// plots if signal yield vs. eta and pt
void genPlots03(std::string fullPath)
{
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c1","c1",1000,600);
    const unsigned int nPadX = 3;
    const unsigned int nPadY = 2;
    c->Divide(nPadX,nPadY);
    const unsigned int nPads=nPadX*nPadY;
    for(unsigned int i=1; i<=nPads; i++)
    {
	TPad* pad= (TPad*)c->cd(i);
	pad->SetTopMargin(0.10);
	pad->SetRightMargin(0.20);
	pad->SetLeftMargin(0.20);
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
    // Do plots
    TPad* pad;
    // -- Plots pt/eta 1
    pad = (TPad*)c->cd(1);
    doPlot1d(t,"hetamu1", "TMath::Max(TMath::Abs(etamu1),TMath::Abs(etamu2))",25,0,2.5,"#eta vs p_{t} of #mu_{1}","|#eta|","");
    pad = (TPad*)c->cd(1+nPadX);
    std::string xTitle="max(|#eta_{1}|,|#eta_{2}|)";
    const double sum = doPlot1dIntegral((TH1F*)gDirectory->GetList()->FindObject("hetamu1"),"hetamu1int","Integrated signal yield, no p_{t} cut",xTitle.c_str(),"");
    pad = (TPad*)c->cd(2);
    doPlot1d(t,"hetamu1pt3", "TMath::Max(TMath::Abs(etamu1),TMath::Abs(etamu2))","ptmu1>3&&ptmu2>3",25,0,2.5,"#eta vs p_{t} of #mu_{1}","|#eta|","");
    pad = (TPad*)c->cd(2+nPadX);
    doPlot1dIntegral((TH1F*)gDirectory->GetList()->FindObject("hetamu1pt3"),"hetamu1pt3int","Integrated signal yield, p_{t} > 3 GeV/c",xTitle.c_str(),"",1./sum);
    pad = (TPad*)c->cd(3);
    doPlot1d(t,"hetamu1pt4", "TMath::Max(TMath::Abs(etamu1),TMath::Abs(etamu2))","ptmu1>4&&ptmu2>4",25,0,2.5,"#eta vs p_{t} of #mu_{1}","|#eta|","");
    pad = (TPad*)c->cd(3+nPadX);
    doPlot1dIntegral((TH1F*)gDirectory->GetList()->FindObject("hetamu1pt4"),"hetamu1pt4int","Integrated signal yield, p_{t} > 4 GeV/c",xTitle.c_str(),"",1./sum);
    // now superimpose them
    c2 = new TCanvas("c2","c2",1000,600);
    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("hetamu1int");
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject("hetamu1pt3int");
    TH1F *h3 = (TH1F*)gDirectory->GetList()->FindObject("hetamu1pt4int");
    h1->SetStats(false);
    h2->SetStats(false);
    h3->SetStats(false);
    h1->SetTitle("Integrated signal yield for #Lambda_{b}");
    h1->SetFillColor(5);
    h2->SetFillColor(3);
    h3->SetFillColor(4);
    h1->Draw();
    h2->Draw("same");
    h3->Draw("same");
    TLegend *leg = new TLegend(0.4,0.6,0.89,0.89);
    leg->AddEntry(h1,"no cut on p_{t}","lf");
    leg->AddEntry(h2,"p_{t} > 3 GeV/c","lf");
    leg->AddEntry(h3,"p_{t} > 5 GeV/c","lf");
    leg->Draw();
}



