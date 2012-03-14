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
#include "TEventList.h"
#include "TMath.h"
#include "setTDRStyle_modified.C"
#include "cuts.C"
#include "utils.h"

TCanvas *c, *c2;

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

void doPlot2d(TTree *t, std::string hname, std::string todraw, std::string cut,
	int nBinsX, double minX, double maxX, int nBinsY, double minY, double maxY,
	std::string title, std::string titleX, std::string unitX, std::string titleY, std::string unitY)
{
    const std::string plotstring = todraw + ">>" + hname + "(" +
	toString(nBinsX) + "," + toString(minX) + "," + toString(maxX) + "," +
	toString(nBinsY) + "," + toString(minY) + "," + toString(maxY) + ")";
    t->Draw(plotstring.c_str(), cut.c_str(),"COLZ");
    TH2F *h = (TH2F*)gDirectory->GetList()->FindObject(hname.c_str());
    h->SetTitle(title.c_str());
    h->GetXaxis()->SetTitle(unitX.size()>0 ? (titleX+" / "+unitX).c_str() : titleX.c_str());
    h->GetYaxis()->SetTitle(unitY.size()>0 ? (titleY+" / "+unitY).c_str() : titleY.c_str());
    return; 
}

void setPadParams(TPad *pad)
{
    pad->SetTopMargin(0.1);
    pad->SetBottomMargin(0.10);
    pad->SetLeftMargin(0.10);
    pad->SetRightMargin(0.30);
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

double calcEff(TH1F* h)
{
    if (h->GetBinContent(0) != 0) cout << "Caution! Underflow bin has enrties." << endl;
    if (h->GetBinContent(h->GetNbinsX()+1) != 0) cout << "Caution! Overflow bin has enrties." << endl;
    double invsum(0);
    for (int i=1; i<= h->GetNbinsX(); i++)
    {
	//cout << "bin " << i << ": " << h->GetBinCenter(i) << " " << h->GetBinContent(i) << endl;
	invsum += h->GetBinContent(i) / h->GetBinCenter(i);
    }
    return h->GetEntries()/invsum;
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

void efficiencyLb_01(std::string fullPath, bool isMC = false)
{
    const bool doAddltPlots(true);
    const bool doAddltPlots2d(true);

    const std::string outpdf("efficiencyLb_01plots");
    const int fVerbose(0);
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Canvas
    c = new TCanvas("c1","c1",1200,900);
    int canvasCdCounter(0), canvasPageCounter(0);
    const unsigned int nPadX = 4;
    const unsigned int nPadY = doAddltPlots ? 4 : 2;
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
    if(fVerbose>-10)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTree
    TTree* t = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << t->GetEntries() << " entries" << endl;

    // General cuts
    Cuts cuts;
    //cuts.selectCut("an03","acc03","HLT_matched_01");
    cuts.selectCut("an05","acc03","HLT_matched_01");
    std::string cutsAnalysis = cuts.getCut();
    // if it is MC make sure we only take events in Lb->JpL0

    if (isMC)
    {
	cutsAnalysis+= "&&isSig==1";
	cout << "MC treatment requested, require isSig==1" << endl;
    }

    //std::string cutssel = cutsAnalysis + "&&isSig==1&&isMCmatch==1&&HLTmatch==1";
    //std::string cutssel = cutsAnalysis + "&&HLTmatch==1";
    std::string cutssel = cutsAnalysis;
    cout << "Full cutstring for the record:" << endl;
    cout << cutssel << endl;
    // Do plots
    const double ptmax(1.0);
    const double etamax(.03);
    const double phimax(.03);
    const double ymax(.03);
    const int nBins(40);

    // Apply cuts

    t->Draw(">>lst",cutssel.c_str());
    TEventList* lst;
    lst = (TEventList*)gDirectory->Get("lst");
    t->SetEventList(lst);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hmlb2", "mlb", "",nBins,5.22,6.04,"m(J/#psi#Lambda)","m","GeV/c^{2}");
    c->Update();
    repositionStatbox("hmlb2");

    // now we calculate the overall efficiency within a given mass range
    const double mass = 5.623;
    const double mass_sigma = 0.022;
    const double nsigmas = 2.5;
    const double mass_lo = mass-nsigmas*mass_sigma;
    const double mass_hi = mass+nsigmas*mass_sigma;
    std::string masscut = "mlb>" + toString(mass_lo) + "&&mlb<" + toString(mass_hi);

    t->Draw(">>lst2",masscut.c_str());
    TEventList* lst2;
    lst2 = (TEventList*)gDirectory->Get("lst2");
    t->SetEventList(lst2);

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"hmlb3", "mlb", "",nBins,5.22,6.04,"m(J/#psi#Lambda), window used for eff","m","GeV/c^{2}");
    c->Update();
    repositionStatbox("hmlb3");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"heff", "Eff", "Eff>0&&Eff<=1",nBins,0,1,"Overall efficiency","Efficiency","");
    c->Update();
    repositionStatbox("heff");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"heffErr", "EffErr", "Eff>0&&Eff<=1",nBins,0,1,"Error of overall efficiency","Eff error","");
    c->Update();
    repositionStatbox("heffErr");

    if (doAddltPlots)
    {
	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	c->cd(canvasCdCounter);
	doPlot1d(t,"heffZoom", "Eff", "Eff>0&&Eff<=1",nBins,0,.2,"Overall efficiency (zoomed in)","Efficiency","");
	c->Update();
	repositionStatbox("heff");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	c->cd(canvasCdCounter);
	doPlot1d(t,"heffErrZoom", "Eff", "Eff>0&&Eff<=1",nBins,0,.2,"Error of overall efficiency (zoomed)","Eff error","");
	c->Update();
	repositionStatbox("heff");
    }

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"heffRelErr", "EffErr/Eff", "Eff>0&&Eff<=1",nBins,0,3,"rel. Error of ov. eff.","rel. error","");
    c->Update();
    repositionStatbox("heffRelErr");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"heffMu", "EffMu", "Eff>0&&Eff<=1",nBins,0,1,"Overall muon efficiency","Efficiency","");
    c->Update();
    repositionStatbox("heffMu");


    if (doAddltPlots)
    {
	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	c->cd(canvasCdCounter);
	doPlot1d(t,"heffMuId1", "EffMuId1", "Eff>0&&Eff<=1",nBins,0,1,"MuId1","Efficiency","");
	c->Update();
	repositionStatbox("heffMuId1");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	c->cd(canvasCdCounter);
	doPlot1d(t,"heffMuTrk1", "EffMuTrk1", "Eff>0&&Eff<=1",nBins,0,1,"MuTrk1","Efficiency","");
	c->Update();
	repositionStatbox("heffMuTrk1");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	c->cd(canvasCdCounter);
	doPlot1d(t,"heffMuTrg1", "EffMuTrg1", "Eff>0&&Eff<=1",nBins,0,1,"MuTrg1","Efficiency","");
	c->Update();
	repositionStatbox("heffMuTrg1");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	c->cd(canvasCdCounter);
	doPlot1d(t,"heffMuId2", "EffMuId2", "Eff>0&&Eff<=1",nBins,0,1,"MuId2","Efficiency","");
	c->Update();
	repositionStatbox("heffMuId2");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	c->cd(canvasCdCounter);
	doPlot1d(t,"heffMuTrk2", "EffMuTrk2", "Eff>0&&Eff<=1",nBins,0,1,"MuTrk2","Efficiency","");
	c->Update();
	repositionStatbox("heffMuTrk2");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	c->cd(canvasCdCounter);
	doPlot1d(t,"heffMuTrg2", "EffMuTrg2", "Eff>0&&Eff<=1",nBins,0,1,"MuTrg2","Efficiency","");
	c->Update();
	repositionStatbox("heffMuTrg2");
    }

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"heffL0", "EffL0", "Eff>0&&Eff<=1",nBins,0,1,"Lambda reco efficiency","Efficiency","");
    c->Update();
    repositionStatbox("heffL0");

    canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
    c->cd(canvasCdCounter);
    doPlot1d(t,"heffCuts", "EffCuts", "Eff>0&&Eff<=1",nBins,0,1,"Cut efficiency","Efficiency","");
    c->Update();
    repositionStatbox("heffCuts");

    if (doAddltPlots2d)
    {
	const int nBins(10); // local change of bin count

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	setPadParams((TPad*)c->cd(canvasCdCounter));
	doPlot2d(t,"heffMuId1_MuId2", "EffMuId2:EffMuId1", "Eff>0&&Eff<=1",nBins,0,1,nBins,0,1,"MuId2 vs. MuId1","MuId1","","MuId2","");
	c->Update();
	repositionPalette("heffMuId1_MuId2");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	setPadParams((TPad*)c->cd(canvasCdCounter));
	doPlot2d(t,"heffMuTrk1_MuTrk2", "EffMuTrk2:EffMuTrk1", "Eff>0&&Eff<=1",nBins,0,1,nBins,0,1,"MuTrk2 vs. MuTrk1","MuTrk1","","MuTrk2","");
	c->Update();
	repositionPalette("heffMuTrk1_MuTrk2");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	setPadParams((TPad*)c->cd(canvasCdCounter));
	doPlot2d(t,"heffMuTrg1_MuTrg2", "EffMuTrg2:EffMuTrg1", "Eff>0&&Eff<=1",nBins,0,1,nBins,0,1,"MuTrg2 vs. MuTrg1","MuTrg1","","MuTrg2","");
	c->Update();
	repositionPalette("heffMuTrg1_MuTrg2");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	setPadParams((TPad*)c->cd(canvasCdCounter));
	doPlot2d(t,"heffMu_L0", "EffL0:EffMu", "Eff>0&&Eff<=1",nBins,0,1,nBins,0,1,"L0 vs. Mu","Mu","","L0","");
	c->Update();
	repositionPalette("heffMu_L0");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	setPadParams((TPad*)c->cd(canvasCdCounter));
	doPlot2d(t,"heffMu_Cuts", "EffCuts:EffMu", "Eff>0&&Eff<=1",nBins,0,1,nBins,0,1,"Cuts vs. Mu","Mu","","Cuts","");
	c->Update();
	repositionPalette("heffMu_Cuts");

	canvaspager(c,outpdf,nPads,canvasCdCounter,canvasPageCounter);
	setPadParams((TPad*)c->cd(canvasCdCounter));
	doPlot2d(t,"heffL0_Cuts", "EffCuts:EffL0", "Eff>0&&Eff<=1",nBins,0,1,nBins,0,1,"Cuts vs. L0","L0","","Cuts","");
	c->Update();
	repositionPalette("heffL0_Cuts");
    }

    //t->SetBranchAddress("mlb", &mlb);
    double lbeff, lbeffmin(1), lbeffmax(0);
    double lbeffErr;
    double lbeffL0, lbeffMu, lbeffCuts;
    double lbeffSum(0);
    double lbeffErrSumSq(0);
    int N(0);
    t->SetBranchAddress("Eff", &lbeff);
    t->SetBranchAddress("EffErr", &lbeffErr);
    t->SetBranchAddress("EffL0", &lbeffL0);
    t->SetBranchAddress("EffMu", &lbeffMu);
    t->SetBranchAddress("EffCuts", &lbeffCuts);
    double pt1m, eta1m, phi1m;
    t->SetBranchAddress("Spt1m", &pt1m);
    t->SetBranchAddress("Seta1m", &eta1m);
    t->SetBranchAddress("Sphi1m", &phi1m);
    double pt2m, eta2m, phi2m;
    t->SetBranchAddress("Spt2m", &pt2m);
    t->SetBranchAddress("Seta2m", &eta2m);
    t->SetBranchAddress("Sphi2m", &phi2m);
    double ptpr, etapr, phipr;
    t->SetBranchAddress("Sptpr", &ptpr);
    t->SetBranchAddress("Setapr", &etapr);
    t->SetBranchAddress("Sphipr", &phipr);
    double ptpi, etapi, phipi;
    t->SetBranchAddress("Sptpi", &ptpi);
    t->SetBranchAddress("Setapi", &etapi);
    t->SetBranchAddress("Sphipi", &phipi);
    double ptl0, etal0, phil0;
    t->SetBranchAddress("ptl0", &ptl0);
    t->SetBranchAddress("etal0", &etal0);
    t->SetBranchAddress("phil0", &phil0);
    // now loop on all events select by the analysis cuts and the mass cut
    for(int i=0; i!=lst2->GetN(); i++)
    {
	t->GetEntry(lst2->GetEntry(i));
	if (lbeff>0. && lbeff <=1.)
	{
	    N++;
	    if (lbeff<lbeffmin) lbeffmin = lbeff;
	    if (lbeff>lbeffmax) lbeffmax = lbeff;
	    //cout << "value " << i << ": " << lbeff << endl;
	    lbeffSum += 1./lbeff;
	    lbeffErrSumSq += square(lbeffErr/lbeff);
	}
	else
	{
	    cout << "Efficiency for entry " << i 
		 << " out of bounds: " << lbeff 
		 << " EffL0: " << lbeffL0
		 << " EffMu: " << lbeffMu
		 << " EffCut: " << lbeffCuts << endl;
	    cout << "  m1: (" << pt1m << "," << eta1m << "," << phi1m << ")"
	         << " m2: (" << pt2m << "," << eta2m << "," << phi2m << ")"
	         << " pr: (" << ptpr << "," << etapr << "," << phipr << ")"
	         << " pi: (" << ptpi << "," << etapi << "," << phipi << ")"
	         << " l0: (" << ptl0 << "," << etal0 << "," << phil0 << ")" << endl;
	}
    }
    const double totaleff = N/lbeffSum;
    const double totaleffErr = N/square(lbeffSum)*TMath::Sqrt(lbeffErrSumSq);
    cout << "No. of events after cuts: " << lst2->GetN() << endl;
    cout << "No. of events whith a valid efficiency value: " << N << endl;
    cout << "Minimum efficiency: " << lbeffmin << endl;
    cout << "Maximum efficiency: " << lbeffmax << endl;
    cout << "Total efficiency: " << totaleff << endl;
    cout << "Error on total efficiency: " << totaleffErr << endl;

    // finalize current page
    c->SaveAs((outpdf + toString(canvasPageCounter) + ".pdf").c_str());

    // Calculate efficiency from bins in histo
    cout << "Total efficiency from histogrammed data" << endl;
    cout << "heff: " << calcEff((TH1F*)gDirectory->GetList()->FindObject("heff")) << endl;
    cout << "heffZoom: " << calcEff((TH1F*)gDirectory->GetList()->FindObject("heffZoom")) << endl;
    /*
    {
	cout << "------------------" << endl;
	double invsum(0);
	TH1F *h = (TH1F*)gDirectory->GetList()->FindObject("heffZoom");
	for (int i=1; i!= h->GetNbinsX(); i++)
	{
	    cout << "bin " << i << ": " << h->GetBinCenter(i) << " " << h->GetBinContent(i) << endl;
	    invsum += h->GetBinContent(i) / h->GetBinCenter(i);
	}
	cout << "Invsum: " << invsum << " " << h->GetEntries()/invsum << endl;
    }*/
}

