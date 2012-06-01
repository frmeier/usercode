#include<string>
#include<vector>
#include<memory>

#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TEventList.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"

#include "utils.h"
#include "setTDRStyle_modified.C"
#include "Cuts.C"

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::auto_ptr;

bool suppressTitle(true);

class OneHisto
{
    public:
	OneHisto() {};
	OneHisto(const OneHisto& oh)
	{
	    filename = oh.filename;
	    file = oh.file;
	    tree = oh.tree;
	    histo = oh.histo;
	    title = oh.title;
	    cut = oh.cut;
	    treeName = oh.treeName;
	    color = oh.color;
	    scale = oh.scale;
	};
	OneHisto(string f, string tn, string t, string c, Color_t cl, double sc) : filename(f), treeName(tn), title(t), cut(c), color(cl), scale(sc) {};
	
	string filename;
	TFile* file;
	TTree* tree;
	TH1F *histo;
	string treeName;
	string title;
	string cut;
	Color_t color;
	double scale; // ==0 then scale to getEntries, else to this factor w.r.t. first entry
};

TH1F* mkSoverSBhisto(TH1F* hSig, TH1F* hBgr, bool LtoR)
{
    string newname = (string)hSig->GetName() + "_" + (string)hBgr->GetName() + "_SoverSB";
    TH1F* ret = (TH1F*)hSig->Clone(newname.c_str());
    const int nBins = hSig->GetNbinsX();
    for (int i=1; i!=nBins; i++)
    {
	// calculate integrals including overflow bins
	const double S = (LtoR ? hSig->Integral(0,i) : hSig->Integral(i,nBins+1) );
	const double B = (LtoR ? hBgr->Integral(0,i) : hBgr->Integral(i,nBins+1) );
	const double SB = S+B;
	const double SoverSB = SB > 0 ? S/sqrt(SB) : 0;
	ret->SetBinContent(i, SoverSB);
    }
    return ret;
}

// General routine to plot a set of plots superimposed, normalized to the entires of the first one
void plotMatchNoMatch(vector<OneHisto> hVec, const string title, const string saveAs,  const string toPlot, const string name,
	const int nBins, const double lo, const double hi, const string axisTitle, const string cut = "")
{
    setTDRStyle();
    gStyle->SetOptStat(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleBorderSize(0); 
    double allMax = 0;
    int i = 0;
    bool doNormEntries(false);
    double normTo = 0;
    for (vector<OneHisto>::iterator it=hVec.begin(); it!=hVec.end(); it++, i++)
    {
	it->file = TFile::Open(it->filename.c_str());
	if (it->file==0) { cout << "File " << it->filename << " doesn't exist - exiting :-(" << endl; return; }
	it->tree = (TTree*)it->file->Get("events");
	if (it->tree==0) { cout << "Tree " << it->treeName << "not contained in file " << it->filename << " - exiting :-(" << endl; return; }
	const string curName = name+toString(i);
	it->histo = new TH1F(curName.c_str(), suppressTitle ? "" : title.c_str(), nBins, lo, hi);
	const string curCut = cut + ((cut.size()!=0 && it->cut.size()!=0) ? "&&" : "") + it->cut;
	it->tree->Draw((toPlot+">>"+curName).c_str(), curCut.c_str(), i==0 ? "" : "same");
	it->histo->SetLineColor(it->color);
	it->histo->SetFillColor(0);
	it->histo->SetFillStyle(0);
	it->histo->GetXaxis()->SetTitle(axisTitle.c_str());
	it->histo->GetYaxis()->SetTitle(entriesPerBin(nBins, lo, hi, "").c_str());
	it->histo->GetYaxis()->SetTitleOffset(1.5);
	it->histo->GetXaxis()->SetNdivisions(509);
	if (i==0)
	{
	    if (it->scale < 0.)
	    {
		doNormEntries = true;
		normTo = it->histo->GetEntries();
	    }
	    else
	    {
		doNormEntries = false;
		normTo = it->scale;
	    }
	    allMax = it->histo->GetMaximum();
	    it->histo->SetLineWidth(4);
	}
	else
	{
	    if(doNormEntries)
		it->histo->Scale( normTo / (double)it->histo->GetEntries());
	    else
		it->histo->Scale( normTo * it->scale);
	    if(it->histo->GetMaximum() > allMax) allMax = it->histo->GetMaximum();
	    it->histo->SetLineWidth(2);
	}
	it->histo->SetMinimum(0);
    }

    hVec[0].histo->Draw("same");

    gPad->SetTopMargin(suppressTitle ? 0.02 : 0.08);
    gPad->SetLeftMargin(0.18);
    gPad->Update();
    const double legWidth(.25);
    const double legHeightPerEntry(.04);
    //TLegend *leg = new TLegend(0.65,0.80,0.90,0.90);
    TLegend *leg = new TLegend(0.90-legWidth,0.90-i*legHeightPerEntry, 0.90, 0.90);
    for (vector<OneHisto>::iterator it=hVec.begin(); it!=hVec.end(); it++, i++)
    {
	it->histo->SetMaximum(allMax*1.1);
	leg->AddEntry(it->histo, it->title.c_str(), "l");
    }
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetLineStyle(0);
    leg->SetLineColor(0);
    leg->SetShadowColor(0);
    leg->Draw();

    // save pad if requested
    if (saveAs.size() != 0) gPad->SaveAs(saveAs.c_str());

    // garbage collection
    delete leg;
    for (vector<OneHisto>::iterator it=hVec.begin(); it!=hVec.end(); it++, i++)
    {
	it->histo->Delete();
	it->tree->Delete();
	it->file->Delete();
    }

}

void plotMatchNoMatch(string filename, const string toPlot, const string name, const int nBins, const double lo, const double hi, const string cut = "")
{
    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
    }
    TTree *tree = (TTree*)f->Get("events");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
    }
    tree->Draw(">>lst", cut.c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");
    cout << "Selection has " << lst->GetN() << " entries" << endl;
    tree->SetEventList(lst);

    const string title = "truth matched comparison: " + toPlot;
    TH1F *hgood = new TH1F("hgood", title.c_str(), nBins, lo, hi);
    TH1F *hbad = new TH1F("hbad", "hbad", nBins, lo, hi);

    tree->Draw((toPlot+">>hgood").c_str(), "isMCmatch==1");
    tree->Draw((toPlot+">>hbad").c_str(), "isMCmatch==0");
    //TH1F *hSoverSB = mkSoverSBhisto(hgood, hbad, LtoR);

    hgood->Draw();
    hbad->Draw("same");
    hgood->SetLineColor(3);
    hbad->SetLineColor(2);
    const double scale = hgood->GetEntries()/hbad->GetEntries();
    hbad->Scale(scale);
    if (hbad->GetMaximum() > hgood->GetMaximum()) hbad->SetMaximum(hbad->GetMaximum());
    double max = (hbad->GetMaximum() > hgood->GetMaximum()) ? hbad->GetMaximum() : hgood->GetMaximum();
    hgood->SetMaximum(max);
    hbad->SetMaximum(max);
    hgood->SetMinimum(0);
    hbad->SetMinimum(0);
    cout << "hbad->GetMaximum(): " << hbad->GetMaximum() << endl;
    cout << "hgood->GetMaximum(): " << hgood->GetMaximum() << endl;

    TLegend *leg = new TLegend(0.4,0.5,0.65,0.6);
    leg->AddEntry(hgood, "matched", "l");
    leg->AddEntry(hbad,("unmatched (scaled: "+roundToString(scale,2)+")").c_str(), "l");
    leg->Draw();

    gPad->Update();
    gPad->SaveAs(("plotMatchNoMatch_"+name+".pdf").c_str());

    // draw SoverSB
    /*
    hSoverSB->Draw("same");
    hSoverSB->SetLineColor(1);
    hSoverSB->Scale(hgood->GetMaximum()/hSoverSB->GetMaximum());
    */

    gPad->Update();
    gPad->SaveAs(("plotMatchNoMatch_"+name+".pdf").c_str());
}

/*
void plotMatchNoMatch2d(string filename, const string toPlot, const string name, const int nBins, const double lo, const double hi, const string cut = "")
{
    // data
    TFile *f = TFile::Open(filename.c_str());
    if (f==0)
    {
	cout << "Problem: File \"" << filename << "\" not found. Exiting..." << endl;
    }
    TTree *tree = (TTree*)f->Get("events");
    if (tree==0)
    {
	cout << "Unable to get tree from file. Exitng..." << endl;
    }
    tree->Draw(">>lst", cut.c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");
    cout << "Selection has " << lst->GetN() << " entries" << endl;
    tree->SetEventList(lst);

    const string title = "truth matched comparison: " + toPlot;
    TH1F *hgood = new TH1F("hgood", title.c_str(), nBins, lo, hi);
    TH1F *hbad = new TH1F("hbad", "hbad", nBins, lo, hi);

    tree->Draw((toPlot+">>hgood").c_str(), "isMCmatch==1");
    tree->Draw((toPlot+">>hbad").c_str(), "isMCmatch==0");

    hgood->Draw();
    hbad->Draw("same");
    hgood->SetLineColor(3);
    hbad->SetLineColor(2);
    const double scale = hgood->GetEntries()/hbad->GetEntries();
    hbad->Scale(scale);
    if (hbad->GetMaximum() > hgood->GetMaximum()) hbad->SetMaximum(hbad->GetMaximum());
    double max = (hbad->GetMaximum() > hgood->GetMaximum()) ? hbad->GetMaximum() : hgood->GetMaximum();
    hgood->SetMaximum(max);
    hbad->SetMaximum(max);
    hgood->SetMinimum(0);
    hbad->SetMinimum(0);
    cout << "hbad->GetMaximum(): " << hbad->GetMaximum() << endl;
    cout << "hgood->GetMaximum(): " << hgood->GetMaximum() << endl;

    TLegend *leg = new TLegend(0.4,0.5,0.65,0.6);
    leg->AddEntry(hgood, "matched", "l");
    leg->AddEntry(hbad,("unmatched (scaled: "+roundToString(scale,2)+")").c_str(), "l");
    leg->Draw();
}
*/

void doPlots(const vector<OneHisto> &myHVec, const bool isB0, const string cut = "", const string pdfAddPrefix = "", int nBins = 100)
{
    const double mL0_lo(1.06568), mL0_hi(1.16568);
    const double mKs_lo(0.4476),  mKs_hi(0.5476);

    const string pTitle(isB0  ? "B^{0}" : "#Lambda_{b}");
    const string pName(isB0  ? "B0" : "lb");
    const string jpTitle("J/#psi");
    const string v0Title(isB0 ? "K_{s}" : "#Lambda^{0}");
    const string v0name(isB0  ? "Ks"    : "l0");
    const string otherv0Title(!isB0 ? "K_{s}" : "#Lambda^{0}");
    const string otherv0name(!isB0  ? "Ks"    : "L0");
    const string pdfPrefix("plotMatchNoMatch_"+pName+"_"+(pdfAddPrefix.size() != 0 ? pdfAddPrefix+"_" : ""));

    plotMatchNoMatch(myHVec, "mass "+pTitle, pdfPrefix+"m"+pName+".pdf", "m"+pName, "m"+pName, nBins, 5, 6, valueWithUnit("m("+pTitle+")","GeV/c^{2}"), cut);
    // cut observables
    plotMatchNoMatch(myHVec, "alpha "+v0Title, pdfPrefix+"alphaV0.pdf", "alpha"+v0name, "alpha"+v0name, nBins, 0.0, 0.02, valueWithUnit("#alpha("+v0Title+")","rad"), cut);
    plotMatchNoMatch(myHVec, "d3 "+v0Title, pdfPrefix+"d3V0.pdf", "d3"+v0name, "d3"+v0name, nBins, 0.0, isB0 ? 3.0 : 10.0, valueWithUnit("d_{3}("+v0Title+")","cm"), cut);
    plotMatchNoMatch(myHVec, "d3 significance "+v0Title, pdfPrefix+"d3"+v0name+"_sig.pdf", "d3"+v0name+"/d3E"+v0name, "d3"+v0name+"ig", nBins, 0, 100, valueWithUnit("significance d_{3}("+v0Title+")","#sigmas"), cut);
    plotMatchNoMatch(myHVec, "prob "+v0Title, pdfPrefix+"probV0.pdf", "prob"+v0name, "prob"+v0name, nBins, 0.0, 1.0, valueWithUnit("prob("+v0Title+")",""), cut);
    plotMatchNoMatch(myHVec, "alpha "+pTitle, pdfPrefix+"alpha"+pName+".pdf", "alpha"+pName, "alpha"+pName, nBins, 0.0, 0.4, valueWithUnit("#alpha("+pTitle+")","rad"), cut);
    // B0/Lb observables
    plotMatchNoMatch(myHVec, "pT "+pTitle, pdfPrefix+"pt"+pName+".pdf", "pt"+pName, "pt"+pName, nBins, 0.0, 40, valueWithUnit("p_{T}("+pTitle+")","GeV/c"), cut);
    plotMatchNoMatch(myHVec, "maxdoca "+pTitle, pdfPrefix+"maxdoca"+pName+".pdf", "maxdoca"+pName, "maxdoca"+pName, nBins, 0.0, 0.1, valueWithUnit("maxdoca("+pTitle+")","cm"), cut);
    plotMatchNoMatch(myHVec, "eta "+pTitle, pdfPrefix+"eta"+pName+".pdf", "eta"+pName, "eta"+pName, nBins, -2.5, 2.5, valueWithUnit("#eta("+pTitle+")",""), cut);
    plotMatchNoMatch(myHVec, "prob "+pTitle, pdfPrefix+"prob"+pName+".pdf", "prob"+pName, "prob"+pName, nBins, 0, 1, valueWithUnit("prob("+pTitle+")",""), cut);
    // Jp observables
    plotMatchNoMatch(myHVec, "mass "+jpTitle, pdfPrefix+"mjp.pdf", "mjp", "mjp", nBins, 2.7, 3.5, valueWithUnit("m("+jpTitle+")","GeV/c^{2}"), cut);
    plotMatchNoMatch(myHVec, "pT "+jpTitle, pdfPrefix+"ptjp.pdf", "ptjp", "ptjp", nBins, 0.0, 40, valueWithUnit("p_{T}("+jpTitle+")","GeV/c"), cut);
    plotMatchNoMatch(myHVec, "pT #mu_{1}", pdfPrefix+"ptmu1.pdf", "rpt1m", "rpt1m", nBins, 0.0, 30, valueWithUnit("p_{T}(#mu_{1})","GeV/c"), cut);
    plotMatchNoMatch(myHVec, "pT #mu_{2}", pdfPrefix+"ptmu2.pdf", "rpt2m", "rpt2m", nBins, 0.0, 30, valueWithUnit("p_{T}(#mu_{2})","GeV/c"), cut);
    plotMatchNoMatch(myHVec, "maxdoca "+jpTitle, pdfPrefix+"maxdocajp.pdf", "maxdocajp", "maxdocajp", nBins, 0.0, 0.03, valueWithUnit("maxdoca("+jpTitle+")","cm"), cut);
    plotMatchNoMatch(myHVec, "dR mumu", pdfPrefix+"dRmumu.pdf", "dRmumu", "dRmumu", nBins, 0.0, 1.2, valueWithUnit("dR(#mu#mu)","#DeltaR"), cut);
    plotMatchNoMatch(myHVec, "eta Jp", pdfPrefix+"etaJp.pdf", "etajp", "etajp", nBins, -2.5, 2.5, valueWithUnit("#eta("+jpTitle+")",""), cut);
    // V0 observables
    plotMatchNoMatch(myHVec, "mass "+v0Title, pdfPrefix+"mV0.pdf", "m"+v0name, "m"+v0name, nBins, isB0 ? mKs_lo : mL0_lo, isB0 ? mKs_hi : mL0_hi, valueWithUnit("m("+v0Title+")","GeV/c^{2}"), cut);
    plotMatchNoMatch(myHVec, "pT "+v0Title, pdfPrefix+"ptv0.pdf", "pt"+v0name, "pt"+v0name, nBins, 0.0, 20, valueWithUnit("p_{T}("+v0Title+")","GeV/c"), cut);
    plotMatchNoMatch(myHVec, "maxdoca "+v0Title, pdfPrefix+"maxdocaV0.pdf", "maxdoca"+v0name, "maxdoca"+v0name, nBins, 0.0, 0.1, valueWithUnit("maxdoca("+v0Title+")","cm"), cut);
    plotMatchNoMatch(myHVec, "dR "+v0Title+jpTitle, pdfPrefix+"dRV0jp.pdf", "dR"+v0name+"jp", "dR"+v0name+"jp", nBins, 0.0, 1.2, valueWithUnit("dR("+v0Title+jpTitle+")","#DeltaR"), cut);
    plotMatchNoMatch(myHVec, toString("dR ")+(isB0 ? "pipi" : "prpi"), pdfPrefix+"dRV0pp.pdf", toString("dR")+(isB0 ? "pipi" : "prpi"), toString("dR")+(isB0 ? "pipi" : "prpi"), nBins, 0.0, 1.0, valueWithUnit(toString("dR(")+(isB0 ? "#pi#pi" : "p#pi")+")","#DeltaR"), cut);
    if (isB0)
    {
	plotMatchNoMatch(myHVec, "pi_1", pdfPrefix+"ptpi1.pdf", "(rptpi1>rptpi2)?rptpi1:rptpi2", "rptpi1", nBins, 0.0, 15.0, valueWithUnit("p_{T}(#pi_{1})","GeV/c"), cut);
	plotMatchNoMatch(myHVec, "pi_2", pdfPrefix+"ptpi2.pdf", "(rptpi1>rptpi2)?rptpi2:rptpi1", "rptpi2", nBins, 0.0, 15.0, valueWithUnit("p_{T}(#pi_{2})","GeV/c"), cut);
    }
    else
    {
	plotMatchNoMatch(myHVec, "pi", pdfPrefix+"ptpi.pdf", "rptpi", "rptpi", nBins, 0.0, 4.0, valueWithUnit("p_{T}(#pi)","GeV/c"), cut);
	plotMatchNoMatch(myHVec, "pr", pdfPrefix+"ptpr.pdf", "rptpr", "rptpr", nBins, 0.0, 15.0, valueWithUnit("p_{T}(p)","GeV/c"), cut);
    }
    // other observables
    //plotMatchNoMatch(myHVec, "isoDocatrk ", pdfPrefix+"isoDocatrk.pdf", "isoDocatrk", "isoDocatrk", nBins, 0.0, 0.1, valueWithUnit("isoDocatrk","cm"), cut);
    plotMatchNoMatch(myHVec, "PvLip ", pdfPrefix+"PvLip.pdf", "PvLip", "PvLip", nBins, 0.0, 0.03, valueWithUnit("PvLip","cm"), cut);
    plotMatchNoMatch(myHVec, "mass hypothesis "+otherv0Title, pdfPrefix+"mV0hypo.pdf", otherv0name+"hypo", otherv0name+"hypo", nBins, !isB0 ? mKs_lo : mL0_lo , !isB0 ? mKs_hi : mL0_hi , valueWithUnit("m_{hypo}("+otherv0Title+")","GeV/c^{2}"), cut);
}

void plotSome()
{
    const string path("../data/");
    vector<OneHisto> myLbVec, myB0Vec;
    // Lb files
    myLbVec.push_back(OneHisto(path+"run385__run387.root", "events", "#Lambda_{b} matched", "isMCmatch==1", 1, 1.0));
    myLbVec.push_back(OneHisto(path+"run385__run387.root", "events", "#Lambda_{b} unmatched", "isMCmatch==0", 2, 13.1));
    myLbVec.push_back(OneHisto(path+"run391.root", "events", "Bgr B^{0}", "", 4, 15.0));
    myLbVec.push_back(OneHisto(path+"run389.root", "events", "Bgr B_{s}", "", 7, 127));
    myLbVec.push_back(OneHisto(path+"run390.root", "events", "Bgr B^{+}", "", 8, 36.4));
    myLbVec.push_back(OneHisto(path+"run388.root", "events", "Bgr J/#psi prompt", "", 6, 18.2));
    // B0 files
    myB0Vec.push_back(OneHisto(path+"run398__run399.root", "events", "B^{0} matched", "isMCmatch==1", 1, 1.0));
    myB0Vec.push_back(OneHisto(path+"run398__run399.root", "events", "B^{0} unmatched", "isMCmatch==0", 2, 18.9));
    myB0Vec.push_back(OneHisto(path+"run403.root", "events", "Bgr #Lambda_{b}", "", 4, 17.2));
    myB0Vec.push_back(OneHisto(path+"run401.root", "events", "Bgr B_{s}", "", 7, 9.44));
    myB0Vec.push_back(OneHisto(path+"run402.root", "events", "Bgr B^{+}", "", 8, 3.78));
    myB0Vec.push_back(OneHisto(path+"run400.root", "events", "Bgr J/#psi prompt", "", 6, 3.56));

    //const string mCutLb = "mlb>5.55&&mlb<5.70";
    //const string mCutB0 = "mB0>5.20&&mB0<5.40";
    const string mCutLb = "mlb>5.00&&mlb<6.00";
    const string mCutB0 = "mB0>5.00&&mB0<6.00";
    Cuts anaCutTrgBarrel, anaCutTrgDispl, anaCutLb, anaCutB0;
    anaCutTrgBarrel.selectCut("HLT_jpsiBarrel", "HLT_matched");
    anaCutTrgDispl.selectCut("HLT_jpsiDispl", "HLT_matched");
    anaCutLb.selectCut("acc05Lb", "muSoft", "lb11");
    anaCutB0.selectCut("acc05B0", "muSoft", "B005");

    const string simBarrelTrg = "ptjp>10&&TMath::Abs(reta1m)<1.25&&TMath::Abs(reta2m)<1.25";
    const string simDisplTrg  = "ptjp>6.9";

    vector<OneHisto> myLbVecBarrel = myLbVec;
    myLbVecBarrel[0].cut += "&&"+anaCutTrgBarrel.getCut();
    myLbVecBarrel[1].cut += "&&"+anaCutTrgBarrel.getCut();
    myLbVecBarrel[2].cut += simBarrelTrg;
    myLbVecBarrel[3].cut += simBarrelTrg;
    myLbVecBarrel[4].cut += simBarrelTrg;
    myLbVecBarrel[5].cut += simBarrelTrg;

    vector<OneHisto> myLbVecDispl = myLbVec;
    myLbVecDispl[0].cut += "&&"+anaCutTrgDispl.getCut();
    myLbVecDispl[1].cut += "&&"+anaCutTrgDispl.getCut();
    myLbVecDispl[2].cut += simDisplTrg;
    myLbVecDispl[3].cut += simDisplTrg;
    myLbVecDispl[4].cut += simDisplTrg;
    myLbVecDispl[5].cut += simDisplTrg;

    vector<OneHisto> myB0VecBarrel = myB0Vec;
    myB0VecBarrel[0].cut += "&&"+anaCutTrgBarrel.getCut();
    myB0VecBarrel[1].cut += "&&"+anaCutTrgBarrel.getCut();
    myB0VecBarrel[2].cut += simBarrelTrg;
    myB0VecBarrel[3].cut += simBarrelTrg;
    myB0VecBarrel[4].cut += simBarrelTrg;
    myB0VecBarrel[5].cut += simBarrelTrg;

    vector<OneHisto> myB0VecDispl = myB0Vec;
    myB0VecDispl[0].cut += "&&"+anaCutTrgDispl.getCut();
    myB0VecDispl[1].cut += "&&"+anaCutTrgDispl.getCut();
    myB0VecDispl[2].cut += simDisplTrg;
    myB0VecDispl[3].cut += simDisplTrg;
    myB0VecDispl[4].cut += simDisplTrg;
    myB0VecDispl[5].cut += simDisplTrg;

    doPlots(myLbVec, false, "", "", 100);
    //doPlots(myLbVec, false, mCutLb, "massWindow", 100);
    doPlots(myLbVecBarrel, false, mCutLb+"&&"+anaCutLb.getCut(), "anaCutBarrel", 50);
    doPlots(myLbVecDispl, false, mCutLb+"&&"+anaCutLb.getCut(), "anaCutDispl", 50);

    doPlots(myB0Vec, true, "", "");
    //doPlots(myB0Vec, true, mCutB0, "massWindow", 100);
    doPlots(myB0VecBarrel, true, mCutB0+"&&"+anaCutB0.getCut(), "anaCutBarrel", 50);
    doPlots(myB0VecDispl, true, mCutB0+"&&"+anaCutB0.getCut(), "anaCutDispl", 50);
}


