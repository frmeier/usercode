#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TLegend.h"

#include "utils.h"
#include "setTDRStyle_modified.C"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int nPlotsDrawn(0);

struct FileEntry
{
    FileEntry(string fn, string p, string n, string ti, Color_t mc) : filename(fn), particle(p), name(n), title(ti), markercol(mc) {};
    string filename;
    string particle;
    string name;
    string title;
    Color_t markercol;
    TFile * file;
    TTree * tree;
};

int recoEff01(string filename, string name, string toPlot, int nBins, double lo, double hi, string unit, string cut = "")
{
    setTDRStyle();
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
	return -1;
    }
    TTree *tree = (TTree*)f->Get("genevents");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
	return -2;
    }
    TCanvas *c1 = new TCanvas("c1","c1",400,400);
    //tree->Draw("ctB0>>h1(20,0,16e-12)","");
    string plotstring = "(" + toString(nBins) + "," + toString(lo) + "," + toString(hi) + ")";
    tree->Draw((toPlot+">>h1"+plotstring).c_str(),cut.c_str());
    TH1F *h1 = (TH1F*)gDirectory->GetList()->FindObject("h1");
    h1->Sumw2();
    h1->SetMarkerStyle(8);
    h1->SetMinimum(0.0);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.08);
    h1->SetTitle("");
    h1->GetXaxis()->SetTitle((toPlot+" (MC)"+(unit.size()>0 ? " ["+unit+"]" : "")).c_str());
    h1->GetXaxis()->SetLabelOffset(+0.006);
    h1->GetXaxis()->SetTitleOffset(1.00);
    h1->GetXaxis()->SetNdivisions(10505);
    Double_t binwidth = (hi-lo)/nBins;
    h1->GetYaxis()->SetTitle((unit.size()>0 ? "Entries per " + toString(binwidth) + " " + unit : "Entries per bin").c_str());
    TGaxis::SetMaxDigits(4);
    c1->SaveAs(("recoEff01plot_"+name+"_MC.pdf").c_str());

    TCanvas *c2 = new TCanvas("c2","c2",300,800);
    tree->Draw((toPlot+">>h2"+plotstring).c_str(), (cut.size()==0 ? "hasCand==1" : ("hasCand==1&&"+cut).c_str()));
    TH1F *h2 = (TH1F*)gDirectory->GetList()->FindObject("h2");
    h2->Sumw2();
    h2->Divide(h1);
    h2->SetMinimum(0.0);
    h2->SetMaximum(0.3);
    h2->SetMarkerStyle(7);
    h2->Draw("pe");
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.05);
    gPad->SetTopMargin(0.04);
    h2->SetTitle("");
    h2->GetXaxis()->SetTitle((toPlot+" (MC)"+(unit.size()>0 ? " ["+unit+"]" : "")).c_str());
    h2->GetXaxis()->SetLabelOffset(-0.016);
    h2->GetXaxis()->SetTitleOffset(0.40);
    h2->GetXaxis()->SetNdivisions(10505);
    h2->GetYaxis()->SetTitle("Efficiency MC/reco");
    c2->SaveAs(("recoEff01plot_"+name+"_eff.pdf").c_str());

    nPlotsDrawn++;
    //delete c1;
    //delete c2;
    return 0;
}

int recoEff01multi(vector<FileEntry> fileVec, string name, string toPlot, int nBins, double lo, double hi, string unit, string cutAll = "", string cutPass = "")
{
    setTDRStyle();
    TCanvas *c1 = new TCanvas("c1","c1",400,400); // distribution of all (=100%)
    TCanvas *c2 = new TCanvas("c2","c2",300,800); // efficiency
    string plotstring = "(" + toString(nBins) + "," + toString(lo) + "," + toString(hi) + ")";
    int i(0);
    Double_t binwidth = (hi-lo)/nBins;
    double h0Entries(0);
    double maxValue(0);
    // prepare the legend
    TLegend *legend = new TLegend(0.40,0.75,0.83,0.90);
    legend->SetFillStyle(1000);
    legend->SetBorderSize(1.);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);
    TLegend *legendEff = new TLegend(0.20,0.76,0.82,0.90);
    legendEff->SetFillStyle(1000);
    legendEff->SetBorderSize(1.);
    legendEff->SetTextSize(0.04);
    legendEff->SetFillColor(0);
    // do the plots
    for(vector<FileEntry>::iterator it = fileVec.begin(); it!=fileVec.end(); it++)
    {
	// first draw the 100% plot
	c1->cd();
	it->tree->Draw((toPlot+it->particle+">>h"+toString(i)+plotstring).c_str(), cutAll.c_str(), (i!=0 ? "same" : ""));
	TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(("h"+toString(i)).c_str());
	h->SetMarkerStyle(8);
	h->SetMarkerColor(it->markercol);
	h->Sumw2();
	h->SetMinimum(0.0);
	h->SetTitle("");
	h->GetXaxis()->SetTitle((toPlot+" (MC)"+(unit.size()>0 ? " ["+unit+"]" : "")).c_str());
	h->GetXaxis()->SetLabelOffset(+0.006);
	h->GetXaxis()->SetTitleOffset(1.00);
	h->GetXaxis()->SetNdivisions(10505);
	h->GetYaxis()->SetTitle((unit.size()>0 ? "Entries per " + toString(binwidth) + " " + unit : "Entries per bin").c_str());
	legend->AddEntry(h, (it->title+" ("+toString(h->GetEntries())+" entries)").c_str(),"p");
	// now for the pass plot
	c2->cd();
	it->tree->Draw((toPlot+it->particle+">>heff"+toString(i)+plotstring).c_str(), cutPass.c_str(), (i!=0 ? "same" : ""));
	TH1F *heff = (TH1F*)gDirectory->GetList()->FindObject(("heff"+toString(i)).c_str());
	heff->Sumw2();
	heff->Divide(h);
	heff->SetMinimum(0.0);
	heff->SetMaximum(1.0);
	heff->SetMarkerStyle(7);
	heff->SetMarkerColor(it->markercol);
	if (i==0)
	    heff->Draw("pe");
	else
	    heff->Draw("samepe");
	heff->SetTitle("");
	heff->GetXaxis()->SetTitle((toPlot+" (MC)"+(unit.size()>0 ? " ["+unit+"]" : "")).c_str());
	heff->GetXaxis()->SetLabelOffset(-0.016);
	heff->GetXaxis()->SetTitleOffset(0.40);
	heff->GetXaxis()->SetNdivisions(10505);
	heff->GetYaxis()->SetTitle("Efficiency MC/reco");
	legendEff->AddEntry(h, (it->title+" ("+toString(h->GetEntries())+" entries)").c_str(),"p");
	// revisiting h and detect scaling factor
	if (i==0)
	    h0Entries = h->GetEntries();
	else
	{
	    const double scalefactor = h0Entries / h->GetEntries();
	    h->Scale(scalefactor);
	}
	maxValue = std::max(h->GetBinContent(h->GetMaximumBin()),maxValue);
	// counter increase
	i++;
    }
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.05);
    gPad->SetTopMargin(0.04);
    const double scaleMax(1.35);
    for(i=0; i!=fileVec.size(); i++)
	((TH1F*)gDirectory->GetList()->FindObject(("h"+toString(i)).c_str()))->SetMaximum(maxValue*scaleMax);
    TGaxis::SetMaxDigits(4);
    c1->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.08);
    legend->Draw();
    c2->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.05);
    gPad->SetTopMargin(0.04);
    legendEff->Draw();
    c1->SaveAs(("recoEff01plotMulti_"+name+"_MC.pdf").c_str());
    c2->SaveAs(("recoEff01plotMulti_"+name+"_eff.pdf").c_str());

    return 0;
}

void doSomePlotsB0(string filename)
{
    nPlotsDrawn = 0;
    recoEff01(filename, "ctB0"      , "ctB0",  20,    0, 14e-12,"s");
    recoEff01(filename, "ctB0zoom"  , "ctB0",  20,    0, 5e-12,"s");
    recoEff01(filename, "pB0"       , "pB0",   20,    0, 80,    "GeV/c");
    recoEff01(filename, "pB0zoom"   , "pB0",   20,    0, 30,    "GeV/c");
    recoEff01(filename, "pB0inv"    , "1/pB0", 20,    0, .3,    "c/GeV");
    recoEff01(filename, "pB0invzoom", "1/pB0", 20,    0, .1,    "c/GeV");
    recoEff01(filename, "d3dB0"     , "d3dB0", 20,    0, 2.0,   "cm");
    recoEff01(filename, "d3dB0zoom" , "d3dB0", 20,    0, 0.4,   "cm");
    recoEff01(filename, "etaB0"     , "etaB0", 20, -2.5, 2.5,   "");
    recoEff01(filename, "ptB0"      , "ptB0",  20,    0, 50,    "GeV/c");
    recoEff01(filename, "ptB0zoom"  , "ptB0",  20,    0, 20,    "GeV/c");
}

void doSomePlotsLb(string filename)
{
    nPlotsDrawn = 0;
    recoEff01(filename, "ctlb"      , "ctlb",  20,    0, 14e-12,"s");
    recoEff01(filename, "ctlbzoom"  , "ctlb",  20,    0, 5e-12,"s");
    recoEff01(filename, "plb"       , "plb",   20,    0, 80,    "GeV/c");
    recoEff01(filename, "plbzoom"   , "plb",   20,    0, 30,    "GeV/c");
    recoEff01(filename, "plbinv"    , "1/plb", 20,    0, .3,    "c/GeV");
    recoEff01(filename, "plbinvzoom", "1/plb", 20,    0, .1,    "c/GeV");
    recoEff01(filename, "d3dlb"     , "d3dlb", 20,    0, 2.0,   "cm");
    recoEff01(filename, "d3dlbzoom" , "d3dlb", 20,    0, 0.4,   "cm");
    recoEff01(filename, "etalb"     , "etalb", 20, -2.5, 2.5,   "");
    recoEff01(filename, "ptlb"      , "ptlb",  20,    0, 50,    "GeV/c");
    recoEff01(filename, "ptlbzoom"  , "ptlb",  20,    0, 20,    "GeV/c");
}

void doListOfPlots()
{
    vector<FileEntry> fileList;
    fileList.push_back(FileEntry("../data/run276.root", "B0", "B0MC", "B^{0}", 1));
    fileList.push_back(FileEntry("../data/run279.root", "B0", "B0longtauMC", "B^{0} long #tau",4));
    fileList.push_back(FileEntry("../data/run282.root", "lb", "lbMC", "#Lambda_{b}", 2));

    // open files and check their existence
    for(vector<FileEntry>::iterator it = fileList.begin(); it!=fileList.end(); it++)
    {
	it->file = TFile::Open(it->filename.c_str());
	if (it->file == 0)
	{
	    cout << "Problem: File \"" << it->filename << "\" not found. Exiting..." << endl;
	    return;
	}
	it->tree = (TTree*)it->file->Get("genevents");
	if (it->tree==0)
	{
	    cout << "Unable to get tree genevents from file " << it->filename << ". Exiting..." << endl;
	    return;
	}
    }

    // Alter Akzeptanzcut
    //const string cut1m("((TMath::Abs(etamu1)<1.3&&ptmu1>3.3)||(TMath::Abs(etamu1)>=1.3&&TMath::Abs(etamu1)<2.2&&pmu1>2.9)||(TMath::Abs(etamu1)>2.2&&TMath::Abs(etamu1)<2.4&&ptmu1>0.8))");
    //const string cut2m("((TMath::Abs(etamu2)<1.3&&ptmu2>3.3)||(TMath::Abs(etamu2)>=1.3&&TMath::Abs(etamu2)<2.2&&pmu2>2.9)||(TMath::Abs(etamu2)>2.2&&TMath::Abs(etamu2)<2.4&&ptmu2>0.8))");
    //const string cutAll(cut1m+"&&"+cut2m+"&&(vrKs>1&&vrKs<35&&TMath::Abs(vzKs)<100)"); // 8tung: funzt nur bei B0
    // Neuer Akzeptanzcut "tracker50"
    //const string cut1m("((TMath::Abs(etamu1)<1.2&&ptmu1>3.5)||(TMath::Abs(etamu1)>=1.2&&TMath::Abs(etamu1)<1.6&&ptmu1>8-3.75*TMath::Abs(etamu1))||(TMath::Abs(etamu1)>1.6&&TMath::Abs(etamu1)<2.4&&ptmu1>2.0))");
    //const string cut2m("((TMath::Abs(etamu2)<1.3&&ptmu2>3.5)||(TMath::Abs(etamu2)>=1.2&&TMath::Abs(etamu2)<1.6&&ptmu2>8-3.75*TMath::Abs(etamu2))||(TMath::Abs(etamu2)>1.6&&TMath::Abs(etamu2)<2.4&&ptmu2>2.0))");
    // Neuer Akzeptanzcut "global50"
    //const string cut1m("((TMath::Abs(etamu1)<1.2&&ptmu1>4.6-1.45*TMath::Abs(etamu1))||(TMath::Abs(etamu1)>=1.2&&TMath::Abs(etamu1)<1.6&&ptmu1>8.59-4.17*TMath::Abs(etamu1))||(TMath::Abs(etamu1)>1.6&&TMath::Abs(etamu1)<2.4&&ptmu1>3.8-.75*TMath::Abs(etamu1)))");
    //const string cut2m("((TMath::Abs(etamu2)<1.2&&ptmu2>4.6-1.45*TMath::Abs(etamu2))||(TMath::Abs(etamu2)>=1.2&&TMath::Abs(etamu2)<1.6&&ptmu2>8.59-4.17*TMath::Abs(etamu2))||(TMath::Abs(etamu2)>1.6&&TMath::Abs(etamu2)<2.4&&ptmu2>3.8-.75*TMath::Abs(etamu2)))");
    //const string cutAll(cut1m+"&&"+cut2m);
    // simuliert Barrel trigger Akzeptanz
    const string cutAll("(TMath::Abs(etamu1)<1.2&&ptmu1>4.6-1.45*TMath::Abs(etamu1))&&(TMath::Abs(etamu2)<1.2&&ptmu1>4.6-1.45*TMath::Abs(etamu2))");
    const string cutPass(cutAll + "&&hasCand==1");
    /*
    const string cutAll("");
    const string cutPass("hasCand==1");
    */
    recoEff01multi(fileList, "ct"      , "ct",  20,    0, 14e-12,"s", cutAll, cutPass);
    recoEff01multi(fileList, "ctzoom"  , "ct",  20,    0, 5e-12, "s", cutAll, cutPass);
    recoEff01multi(fileList, "p"       , "p",   20,    0, 80,    "GeV/c", cutAll, cutPass);
    recoEff01multi(fileList, "pzoom"   , "p",   20,    0, 30,    "GeV/c", cutAll, cutPass);
    recoEff01multi(fileList, "pinv"    , "1/p", 20,    0, .3,    "c/GeV", cutAll, cutPass);
    recoEff01multi(fileList, "pinvzoom", "1/p", 20,    0, .1,    "c/GeV", cutAll, cutPass);
    recoEff01multi(fileList, "d3d"     , "d3d", 20,    0, 2.0,   "cm", cutAll, cutPass);
    recoEff01multi(fileList, "d3dzoom" , "d3d", 20,    0, 0.4,   "cm", cutAll, cutPass);
    recoEff01multi(fileList, "eta"     , "eta", 20, -2.5, 2.5,   "", cutAll, cutPass);
    recoEff01multi(fileList, "pt"      , "pt",  20,    0, 50,    "GeV/c", cutAll, cutPass);
    recoEff01multi(fileList, "ptzoom"  , "pt",  20,    0, 20,    "GeV/c", cutAll, cutPass);
}

