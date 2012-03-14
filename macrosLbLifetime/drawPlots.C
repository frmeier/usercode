#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include <string>
#include "setTDRStyle_modified.C"

using std::endl;
using std::cout;

void drawOnePlot(std::string name, TCanvas *c)
{
    TH1F *h = (TH1F*)gDirectory->Get(name.c_str());
    h->SetStats(true);
    h->Draw();
    c->Update();
    TPaveStats *st = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
    st->SetX1NDC(0.63);
    st->SetY1NDC(0.72);
    st->SetX2NDC(0.99);
    st->SetY2NDC(0.99);
    st->Draw();
}

void drawplots (std::string filename)
{
    const std::string outpdf = "drawplots";

    // Open file
    TFile* file = TFile::Open(filename.c_str());
    if (file ==0)
    {
        cout << "File " << filename << " not found. Exiting." << endl;
        return;
    }
    cout << "File " << filename << " succesfully opened." << endl;

    setTDRStyle();
    gStyle->SetOptStat(111111);

    TCanvas *c = new TCanvas("c","c",1050,1435);
    //TCanvas *c = new TCanvas("c","c",2100,2970);
    c->Divide(2,2);
    c->cd(1); drawOnePlot("h20",c);
    c->cd(2); drawOnePlot("h21",c);
    c->cd(3); drawOnePlot("h22",c);
    c->cd(4); drawOnePlot("h23",c);
    /*
    c->Divide(3,5);
    c->cd(1); drawOnePlot("htm_vtxdeltaR",c);
    c->cd(2); drawOnePlot("htm_vtxdR",c);
    c->cd(3); drawOnePlot("htm_vtxratioR",c);
    c->cd(4); drawOnePlot("htm_vtx2d",c);
    c->cd(5); drawOnePlot("htm_vtx3d",c);
    c->cd(6); drawOnePlot("htm_prdpt",c);
    c->cd(7); drawOnePlot("htm_prdphi",c);
    c->cd(8); drawOnePlot("htm_prdeta",c);
    c->cd(9); drawOnePlot("htm_prdeltaR",c);
    c->cd(10); drawOnePlot("htm_pidpt",c);
    c->cd(11); drawOnePlot("htm_pidphi",c);
    c->cd(12); drawOnePlot("htm_pideta",c);
    c->cd(13); drawOnePlot("htm_pideltaR",c);
    */
    c->SaveAs((outpdf + ".pdf").c_str());
}
