#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TPaletteAxis.h"
#include <iostream>
#include "setTDRStyle_modified.C"

TFile *fdat;
TTree *tgen, *tdat;
TCanvas *c1;
TH2D *h1, *h2, *h3, *h4;
const double MPROTON = 0.9383;
const double MPION = 0.1396;

void anal02(int ndraw=-1) {
    setTDRStyle();
    gStyle->SetOptStat(112211);
    gStyle->SetPalette(1);
    // Objekte
    c1 = new TCanvas();
    c1->Divide(2,2);
    //fdat = TFile::Open("/media/DR-1/run023.root");
    fdat = TFile::Open("/scratch/frmeier/run026.root");
    tgen = (TTree*)fdat->Get("genevents");  
    tdat = (TTree*)fdat->Get("events");   

    // Histo
    int nbins = 40;
    //double min=0, max=4;
    double min=1.0, max=1.2;
    h1 = new TH2D("h1","cands with cuts&&ptpr>ptpi",nbins,min,max,nbins,min,max);
    h2 = new TH2D("h2","chi2lb/ndoflb<chicut",nbins,min,max,nbins,min,max);
    h3 = new TH2D("h3","chi2l0/ndofl0<chicut",nbins,min,max,nbins,min,max);
    h4 = new TH2D("h4","chi2jp/ndofjp<chicut",nbins,min,max,nbins,min,max);
    Int_t nentries = tdat->GetEntries();
    if(ndraw>=0) nentries = ndraw;
    // Treestuff
    double fml0; tdat->SetBranchAddress("ml0",&fml0);
    double fptpr; tdat->SetBranchAddress("ptpr",&fptpr);
    double fetapr; tdat->SetBranchAddress("etapr",&fetapr);
    double fphipr; tdat->SetBranchAddress("phipr",&fphipr);
    double fptpi; tdat->SetBranchAddress("ptpi",&fptpi);
    double fetapi; tdat->SetBranchAddress("etapi",&fetapi);
    double fphipi; tdat->SetBranchAddress("phipi",&fphipi);
    // for cuts
    int fisSig; tdat->SetBranchAddress("isSig",&fisSig);
    int fid1m; tdat->SetBranchAddress("id1m",&fid1m);
    int fid2m; tdat->SetBranchAddress("id2m",&fid2m);
    double fmjp; tdat->SetBranchAddress("mjp",&fmjp);
    double falphal0; tdat->SetBranchAddress("alphal0",&falphal0);
    double fptl0; tdat->SetBranchAddress("ptl0",&fptl0);
    double fd3l0; tdat->SetBranchAddress("d3l0",&fd3l0);
    double fd3lb; tdat->SetBranchAddress("d3lb",&fd3lb);
    double fd3Elb; tdat->SetBranchAddress("d3Elb",&fd3Elb);
    double fchi2lb; tdat->SetBranchAddress("chi2lb",&fchi2lb);
    double fchi2l0; tdat->SetBranchAddress("chi2l0",&fchi2l0);
    double fchi2jp; tdat->SetBranchAddress("chi2jp",&fchi2jp);
    double fndoflb; tdat->SetBranchAddress("ndoflb",&fndoflb);
    double fndofl0; tdat->SetBranchAddress("ndofl0",&fndofl0);
    double fndofjp; tdat->SetBranchAddress("ndofjp",&fndofjp);

    // tlv
    TLorentzVector tlvPr;
    TLorentzVector tlvPi;
    //Loop
    for (Int_t i = 0; i<nentries; i++)
    {
	// cut (id1m&4)==4&&(id2m&4)==4&&mjp>3.03&&mjp<3.16&&alphal0<0.01&&ml0<1.15&&ptl0>3.0&&d3l0>5&&d3lb/d3Elb>1.&&ptpr>ptpi
	tdat->GetEntry(i);
	tlvPr.SetPtEtaPhiM(fptpr,fetapr,fphipr,MPROTON);
	tlvPi.SetPtEtaPhiM(fptpi,fetapi,fphipi,MPION);
	const TLorentzVector tlvL0 = tlvPr + tlvPi;
	//std::cout << i << " " << fml0 << " " << tlvL0.M() << std::endl;
	if((fid1m&4)==4&(fid2m&4)==4&&fmjp>3.03&&fmjp<3.16&&falphal0<0.01&&fml0<1.15&&fptl0>3.0&&fd3l0>5&&fd3lb/fd3Elb>1.&&fptpr>fptpi)
	{
	    h1->Fill(fml0,tlvL0.M());
	    const double chicut = 10;
	    if(fchi2lb/fndoflb<chicut) h2->Fill(fml0,tlvL0.M());
	    if(fchi2l0/fndofl0<chicut) h3->Fill(fml0,tlvL0.M());
	    if(fchi2jp/fndofjp<chicut) h4->Fill(fml0,tlvL0.M());
	}
    }
    c1->cd(1); h1->Draw("colz");
    c1->cd(2); h2->Draw("colz");
    c1->cd(3); h3->Draw("colz");
    c1->cd(4); h4->Draw("colz");
    c1->Update();

    // Reposition and resize palette
    TPaletteAxis *pal1,*pal2,*pal3,*pal4;
    pal1 = (TPaletteAxis*)h1->GetListOfFunctions()->FindObject("palette");
    pal1->SetX1NDC(0.83); pal1->SetY1NDC(0.14); pal1->SetX2NDC(0.88); pal1->SetY2NDC(0.60);
    pal1->SetLabelSize(.04);
    pal2 = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
    pal2->SetX1NDC(0.83); pal2->SetY1NDC(0.14); pal2->SetX2NDC(0.88); pal2->SetY2NDC(0.60);
    pal2->SetLabelSize(.04);
    pal3 = (TPaletteAxis*)h3->GetListOfFunctions()->FindObject("palette");
    pal3->SetX1NDC(0.83); pal3->SetY1NDC(0.14); pal3->SetX2NDC(0.88); pal3->SetY2NDC(0.60);
    pal3->SetLabelSize(.04);
    pal4 = (TPaletteAxis*)h4->GetListOfFunctions()->FindObject("palette");
    pal4->SetX1NDC(0.83); pal4->SetY1NDC(0.14); pal4->SetX2NDC(0.88); pal4->SetY2NDC(0.60);
    pal4->SetLabelSize(.04);
}

