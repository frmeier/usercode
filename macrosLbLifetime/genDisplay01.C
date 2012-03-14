#include <iostream>
#include <iomanip>
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
#include "TLine.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TEventList.h"
#include "setTDRStyle_modified.C"
#include "TLorentzVector.h"

TCanvas *crz, *crphi, *c2;
TH1F *tmph1, *tmph2, *tmph3;
TH2F *tmph4;

const double scalemu = 2;
const double scalepr = 4;
const double scalepi = 4;

template <typename T>
std::string toString(T i)
{
    ostringstream oss;
    oss << i;
    return oss.str();
}

template <typename T>
T sign(T x)
{
    return (x > 0) - (x < 0);
}

int findRecoEvents(TTree *t, int run, int ls, int evt, int matchno = 0)
{
    int ret = -1;
    int curmatch = 0;
    int trun, tls, tevt;
    t->SetBranchAddress("run", &trun);
    t->SetBranchAddress("LS", &tls);
    t->SetBranchAddress("event", &tevt);
    for (int i=0; i!=t->GetEntries(); i++)
    {
	t->GetEntry(i);
	if(run==trun && ls==tls && evt==tevt)
	{
	    cout << "found" << endl;
	    if (matchno == curmatch)
	    {
		ret = i;
		break;
	    }
	    curmatch++;
	}
    }
    return ret;
}

void drawTrackerRZ(TPad* pad)
{
    pad->cd();
	double r = 0, z =0;
	TLine *l;
	// PXB
	z = 26.5;
	r = 4.4; l = new TLine(-z,r,+z,r); l->Draw();
	r = 7.3; l = new TLine(-z,r,+z,r); l->Draw();
	r = 10.2; l = new TLine(-z,r,+z,r); l->Draw();
	r = -4.4; l = new TLine(-z,r,+z,r); l->Draw();
	r = -7.3; l = new TLine(-z,r,+z,r); l->Draw();
	r = -10.2; l = new TLine(-z,r,+z,r); l->Draw();
	// PXE
	z = 34.5; l = new TLine(z,+6,z,+15); l->Draw();
	z = 46.5; l = new TLine(z,+6,z,+15); l->Draw();
	z = 34.5; l = new TLine(z,-6,z,-15); l->Draw();
	z = 46.5; l = new TLine(z,-6,z,-15); l->Draw();
	z = -34.5; l = new TLine(z,+6,z,+15); l->Draw();
	z = -46.5; l = new TLine(z,+6,z,+15); l->Draw();
	z = -34.5; l = new TLine(z,-6,z,-15); l->Draw();
	z = -46.5; l = new TLine(z,-6,z,-15); l->Draw();
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

void drawEllipse(double x, double y, double r)
{
    TEllipse *e;
    e = new TEllipse(x,y,r);
    e->SetFillStyle(4000);
    e->Draw();
}

void drawTrackerRPhi(TPad* pad)
{
    pad->cd();
    double r = 0, x = 0, y = 0;
    TEllipse *e;
    // PXB
    drawEllipse(x,y,4.4);
    drawEllipse(x,y,7.3);
    drawEllipse(x,y,10.2);
    // TIB
    drawEllipse(x,y,25.5);
    drawEllipse(x,y,33.9);
    drawEllipse(x,y,41.9);
    drawEllipse(x,y,49.8);
}

double etaToTheta(double eta)
{
    return 2*TMath::ATan(TMath::Exp(-eta));
}

void drawEventZR(double vrlb, double vzlb, double vrl0, double vzl0, double pmu1, double etamu1, double pmu2, double etamu2, double ppr, double etapr, double ppi, double etapi, int colors = 0)
{
    Color_t colMu1(1), colMu2(1), colPr(1), colPi(1), colL0(1), colLb(1);
    switch(colors)
    {
	case 0:
	    colLb=1; colL0=2; colMu1=2; colMu2=2; colPr=3; colPi=4; break;
	case 1:
	    colLb=11; colL0=50; colMu1=50; colMu2=50; colPr=8; colPi=9; break;
    }
    cout << etamu1 << " " << etamu2 << " " << etapr << " " << etapi << endl;
    const double thetamu1 = 2*TMath::ATan(TMath::Exp(-etamu1));
    const double thetamu2 = 2*TMath::ATan(TMath::Exp(-etamu2));
    const double thetapr = 2*TMath::ATan(TMath::Exp(-etapr));
    const double thetapi = 2*TMath::ATan(TMath::Exp(-etapi));
    //const double thetamu1 = 2*TMath::ATan(TMath::Exp(-TMath::Abs(etamu1)))*sign(etamu1);
    //const double thetamu2 = 2*TMath::ATan(TMath::Exp(-TMath::Abs(etamu2)))*sign(etamu2);
    //const double thetapr = 2*TMath::ATan(TMath::Exp(-TMath::Abs(etapr)))*sign(etapr);
    //const double thetapi = 2*TMath::ATan(TMath::Exp(-TMath::Abs(etapi)))*sign(etapi);
    // draw the vertices
    TMarker *m;
    const double xlb = vzlb;
    const double ylb = vrlb;
    m = new TMarker(xlb,ylb,7);
    m->SetMarkerColor(colLb);
    m->Draw();
    const double xl0 = vzl0;
    const double yl0 = vrl0;
    m = new TMarker(xl0,yl0,7);
    m->SetMarkerColor(colL0);
    m->Draw();
    // draw the l0 flight line
    TLine *l;
    l = new TLine(vzlb,vrlb,vzl0,vrl0);
    l->SetLineColor(colL0);
    l->SetLineStyle(2);
    l->Draw();
    // draw the muons
    TArrow *a;
    const double xmu1 = xlb + scalemu * pmu1 * TMath::Cos(thetamu1);
    const double ymu1 = ylb + scalemu * pmu1 * TMath::Sin(thetamu1);
    a = new TArrow(xlb,ylb,xmu1,ymu1,.01,">");
    a->SetLineColor(colMu1);
    a->Draw();
    const double xmu2 = xlb + scalemu * pmu2 * TMath::Cos(thetamu2);
    const double ymu2 = ylb + scalemu * pmu2 * TMath::Sin(thetamu2);
    a = new TArrow(xlb,ylb,xmu2,ymu2,.01,">");
    a->SetLineColor(colMu2);
    a->Draw();
    // draw the p and pi
    const double xpr = xl0 + scalepr * ppr * TMath::Cos(thetapr);
    const double ypr = yl0 + scalepr * ppr * TMath::Sin(thetapr);
    a = new TArrow(xl0,yl0,xpr,ypr,.01,">");
    a->SetLineColor(colPr);
    a->Draw();
    const double xpi = xl0 + scalepi * ppi * TMath::Cos(thetapi);
    const double ypi = yl0 + scalepi * ppi * TMath::Sin(thetapi);
    a = new TArrow(xl0,yl0,xpi,ypi,.01,">");
    a->SetLineColor(colPi);
    a->Draw();
}

void drawEventRPhi(double vxlb, double vylb, double vxl0, double vyl0, double pmu1, double phimu1, double pmu2, double phimu2, double ppr, double phipr, double ppi, double phipi, int colors = 0)
{
    Color_t colMu1(1), colMu2(1), colPr(1), colPi(1), colL0(1), colLb(1);
    switch(colors)
    {
	case 0:
	    colLb=1; colL0=2; colMu1=2; colMu2=2; colPr=3; colPi=4; break;
	case 1:
	    colLb=11; colL0=50; colMu1=50; colMu2=50; colPr=8; colPi=9; break;
    }
    // draw the vertices
    TMarker *m;
    const double xlb = vxlb;
    const double ylb = vylb;
    m = new TMarker(xlb,ylb,7);
    m->SetMarkerColor(colLb);
    m->Draw();
    const double xl0 = vxl0;
    const double yl0 = vyl0;
    m = new TMarker(xl0,yl0,7);
    m->SetMarkerColor(colL0);
    m->Draw();
    // draw the l0 flight line
    TLine *l;
    l = new TLine(vxlb,vylb,vxl0,vyl0);
    l->SetLineColor(colL0);
    l->SetLineStyle(2);
    l->Draw();
    // draw the muons
    TArrow *a;
    const double xmu1 = xlb + scalemu * pmu1 * TMath::Cos(phimu1);
    const double ymu1 = ylb + scalemu * pmu1 * TMath::Sin(phimu1);
    a = new TArrow(xlb,ylb,xmu1,ymu1,.01,">");
    a->SetLineColor(colMu1);
    a->Draw();
    const double xmu2 = xlb + scalemu * pmu2 * TMath::Cos(phimu2);
    const double ymu2 = ylb + scalemu * pmu2 * TMath::Sin(phimu2);
    a = new TArrow(xlb,ylb,xmu2,ymu2,.01,">");
    a->SetLineColor(colMu2);
    a->Draw();
    // draw the p and pi
    const double xpr = xl0 + scalepr * ppr * TMath::Cos(phipr);
    const double ypr = yl0 + scalepr * ppr * TMath::Sin(phipr);
    a = new TArrow(xl0,yl0,xpr,ypr,.01,">");
    a->SetLineColor(colPr);
    a->Draw();
    const double xpi = xl0 + scalepi * ppi * TMath::Cos(phipi);
    const double ypi = yl0 + scalepi * ppi * TMath::Sin(phipi);
    a = new TArrow(xl0,yl0,xpi,ypi,.01,">");
    a->SetLineColor(colPi);
    a->Draw();
}

double NDCtoUserX(TPad* pad, double x)
{
    //pad->Update();
    return (x*(pad->GetX2()-pad->GetX1())+pad->GetX1());
}

double NDCtoUserY(TPad* pad, double y)
{
    //pad->Update();
    return (y*(pad->GetY2()-pad->GetY1())+pad->GetY1());
}

// write some TLatex text. Return value is new position after the text
double writeLatex(TPad* pad, double x, double y, std::string text, double scale = 1.)
{
    //pad->Update();
    pad->cd();
    TLatex *tl = new TLatex(NDCtoUserX(pad,x), NDCtoUserY(pad,y), text.c_str());
    tl->SetTextSize(.040*scale);
    tl->Draw();
    cout << "latex: " << tl->GetXsize() << " " << tl->GetX() << " " << tl->GetYsize()
         << " pad: " << pad->GetX1() << " " << pad->GetX2() << endl;
    return (x+tl->GetXsize()/(pad->GetX2()-pad->GetX1()));
}

// draw a horizontal arrow of a certain length in pad coordinates at xy in NDC
// returns the new x in NDC
double drawArrowLength(TPad* pad, double x, double y, double length, Color_t color = 1)
{
    const double x1 = NDCtoUserX(pad,x);
    const double y1 = NDCtoUserY(pad,y);
    const double x2 = x1 + length;
    const double y2 = y1;
    TArrow *a = new TArrow(x1,y1,x2,y2,.01,">");
    a->SetLineColor(color);
    a->Draw();
    //pad->Update();
    cout << "Arrow: " << length << " " << pad->GetX2()-pad->GetX1() << endl;
    return x + (length/(pad->GetX2()-pad->GetX1()));
}

// draw a legend for the momenta of the arrows
void drawArrowLegend(TPad* pad, double scale = 1.)
{
    pad->Update();

    double ndcx=.1;
    double ndcy=.1;
    double space = .01*scale;
    double textheight = .03*scale;

    ndcx = space + writeLatex(pad, ndcx, ndcy, "#mu:", scale);
    ndcx = space + drawArrowLength(pad, ndcx, ndcy+textheight*.5, scale*scalemu*10, 2);
    ndcx = space + writeLatex(pad, ndcx, ndcy, toString(10*scale) + " GeV/c", scale);

    ndcx += space;
    ndcx = space + writeLatex(pad, ndcx, ndcy, "p:", scale);
    ndcx = space + drawArrowLength(pad, ndcx, ndcy+textheight*.5, scale*scalepr, 3);
    ndcx = space + writeLatex(pad, ndcx, ndcy, toString(1*scale) + " GeV/c", scale);

    ndcx += space;
    ndcx = space + writeLatex(pad, ndcx, ndcy, "#pi:", scale);
    ndcx = space + drawArrowLength(pad, ndcx, ndcy+textheight*.5, scale*scalepi, 4);
    ndcx = space + writeLatex(pad, ndcx, ndcy, toString(1*scale) + " GeV/c", scale);
}

void setPadMargins(TPad* pad, double t, double b, double l, double r)
{
    pad->SetTopMargin(t);
    pad->SetBottomMargin(b);
    pad->SetRightMargin(l);
    pad->SetLeftMargin(r);
}

// draws the decay vertices of the Lambda_s
void genDisplay(std::string fullPath, int iGen, int iReco = -1, bool savePdf = false)
{
    const int fVerbose(0);
    setTDRStyle();
    //gStyle->SetOptStat(112211);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    // Canvases
    crz = new TCanvas("crz","crz",1000,600);
    TPad* padrz = (TPad*)crz->cd(1);
    setPadMargins(padrz, 0.02, 0.05, 0.02, 0.05);
    crphi = new TCanvas("crphi","crphi",600,600);
    TPad* padrphi = (TPad*)crphi->cd(1);
    setPadMargins(padrphi, 0.02, 0.05, 0.02, 0.05);
    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTree with GenEvents
    TTree* tgen = (TTree*) f->Get("genevents");
    if(fVerbose>0) cout << "Got TTree with " << tgen->GetEntries() << " entries" << endl;
    // Do a cut, if needed
    tgen->Draw(">>lst","ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<2.5&&TMath::Abs(etamu2)<2.5");
    TEventList *lst;
    lst = (TEventList*)gDirectory->Get("lst");
    tgen->SetEventList(lst);

    // Get TTree with iRecoEvents
    TTree* treco = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << treco->GetEntries() << " entries" << endl;

    // Do plots
    const int nbins(30);
    const double rangeR(30), rangeZ(100); // standard
    //const double rangeR(60), rangeZ(200); // zoomed out
    //const double rangeR(100), rangeZ(200); // zoomed in
    //const double rangeR(12), rangeZ(30); // zoomed further in
    padrz->cd();
    TH2F *hrz = new TH2F("hdisp","", nbins,-rangeZ,rangeZ,nbins,-rangeR,rangeR);
    hrz->Draw();
    padrphi->cd();
    TH2F *hrphi = new TH2F("hdisp","", nbins,-rangeR,rangeR,nbins,-rangeR,rangeR);
    hrphi->Draw();

    // add tracker
    padrz->Modified();
    padrz->Update();
    drawTrackerRZ(padrz);
    padrphi->Modified();
    padrphi->Update();
    drawTrackerRPhi(padrphi);

    // set branch addresses for gentree
    double vrl0,vxl0,vyl0,vzl0,vrlb,vxlb,vylb,vzlb;
    double pmu1,phimu1,etamu1,pmu2,phimu2,etamu2;
    double ppr,ppi,phipr,phipi,etapr,etapi;
    tgen->SetBranchAddress("vrl0",&vrl0);
    tgen->SetBranchAddress("vxl0",&vxl0);
    tgen->SetBranchAddress("vyl0",&vyl0);
    tgen->SetBranchAddress("vzl0",&vzl0);
    tgen->SetBranchAddress("vrlb",&vrlb);
    tgen->SetBranchAddress("vxlb",&vxlb);
    tgen->SetBranchAddress("vylb",&vylb);
    tgen->SetBranchAddress("vzlb",&vzlb);
    tgen->SetBranchAddress("pmu1",&pmu1);
    tgen->SetBranchAddress("etamu1",&etamu1);
    tgen->SetBranchAddress("phimu1",&phimu1);
    tgen->SetBranchAddress("pmu2",&pmu2);
    tgen->SetBranchAddress("etamu2",&etamu2);
    tgen->SetBranchAddress("phimu2",&phimu2);
    tgen->SetBranchAddress("ppr",&ppr);
    tgen->SetBranchAddress("etapr",&etapr);
    tgen->SetBranchAddress("phipr",&phipr);
    tgen->SetBranchAddress("ppi",&ppi);
    tgen->SetBranchAddress("etapi",&etapi);
    tgen->SetBranchAddress("phipi",&phipi);
    int genrun, genls, genevt;
    tgen->SetBranchAddress("run",&genrun);
    tgen->SetBranchAddress("LS",&genls);
    tgen->SetBranchAddress("event",&genevt);

    // set branch addresses for recotree
    double rvrl0,rvxl0,rvyl0,rvzl0,rvrlb,rvxlb,rvylb,rvzlb;
    double rpmu1,retamu1,rphimu1,rpmu2,retamu2,rphimu2;
    double rppr,rppi,retapr,retapi,rphipr,rphipi;
    treco->SetBranchAddress("vrl0",&rvrl0);
    treco->SetBranchAddress("vxl0",&rvxl0);
    treco->SetBranchAddress("vyl0",&rvyl0);
    treco->SetBranchAddress("vzl0",&rvzl0);
    treco->SetBranchAddress("vrlb",&rvrlb);
    treco->SetBranchAddress("vxlb",&rvxlb);
    treco->SetBranchAddress("vylb",&rvylb);
    treco->SetBranchAddress("vzlb",&rvzlb);
    treco->SetBranchAddress("rpt1m",&rpmu1);
    treco->SetBranchAddress("reta1m",&retamu1);
    treco->SetBranchAddress("rphi1m",&rphimu1);
    treco->SetBranchAddress("rpt2m",&rpmu2);
    treco->SetBranchAddress("reta2m",&retamu2);
    treco->SetBranchAddress("rphi2m",&rphimu2);
    treco->SetBranchAddress("rptpr",&rppr);
    treco->SetBranchAddress("retapr",&retapr);
    treco->SetBranchAddress("rphipr",&rphipr);
    treco->SetBranchAddress("rptpi",&rppi);
    treco->SetBranchAddress("retapi",&retapi);
    treco->SetBranchAddress("rphipi",&rphipi);

    drawArrowLegend(padrz);
    drawArrowLegend(padrphi,.8);
    // write the event indices
    writeLatex(padrz,.1,.9,("iGen: " + toString(iGen)).c_str());
    writeLatex(padrphi,.1,.9,("iGen: " + toString(iGen)).c_str());
    if (iReco >= 0)
    {
	writeLatex(padrz,.1,.85,("iReco: " + toString(iReco)).c_str());
	writeLatex(padrphi,.1,.85,("iReco: " + toString(iReco)).c_str());
    }

    if (iGen<0) return;
    if (iGen>tgen->GetEntries()) return;
    if (iReco>treco->GetEntries()) return;
    tgen->GetEntry(iGen);
    cout << "iGen: " << iGen << " - " << genrun << " " << genls << " " << genevt << endl;
    padrz->cd();
    drawEventZR(vrlb, vzlb, vrl0, vzl0, pmu1, etamu1, pmu2, etamu2, ppr, etapr, ppi, etapi, 0);
    padrphi->cd();
    drawEventRPhi(vxlb, vylb, vxl0, vyl0, pmu1, phimu1, pmu2, phimu2, ppr, phipr, ppi, phipi, 0);
    if(iReco>=0)
    {
        treco->GetEntry(iReco);
        cout << "iReco: " << iReco << endl;
	padrz->cd();
        drawEventZR(rvrlb, rvzlb, rvrl0, rvzl0, rpmu1, retamu1, rpmu2, retamu2, rppr, retapr, rppi, retapi, 1);
	padrphi->cd();
        drawEventRPhi(rvxlb, rvylb, rvxl0, rvyl0, rpmu1, rphimu1, rpmu2, rphimu2, rppr, rphipr, rppi, rphipi, 1);
	cout << phimu1 << " " << rphimu1 << endl;
    }

    if(savePdf)
    {
	padrz->Update();
	padrphi->Update();
	std::string filename = "genDisplay_" + toString(iGen);
	if(iReco >= 0) filename += "_" + toString(iReco);
	crz->SaveAs((filename+"_rz.pdf").c_str());
	crphi->SaveAs((filename+"_rphi.png").c_str()); // weiss der Geier warum auf PDF die Kreise und die Skala fehlen...
    }

    return;
}

void listMatchedEvents(std::string fullPath, bool matchesonly = false)
{
    // Essentially a sort of diff for reco and gen tree
    // May get slow for large trees as the search for the matching entry in the second tree
    // traverses the whole tree. This is because for merged trees the order is not guaranteed
    const int fVerbose(1);
    // Open file
    TFile *f = TFile::Open(fullPath.c_str());
    if (f==0)
    {
	cout << "File " << fullPath << " not found -- exiting" << endl;
	return;
    }
    if(fVerbose>0)
	cout << "Succesfully opened file " << fullPath << endl;
    // Get TTree with GenEvents
    TTree* tgen = (TTree*) f->Get("genevents");
    if(fVerbose>0) cout << "Got TTree with " << tgen->GetEntries() << " entries" << endl;
    // Do a cut, if needed
    //tgen->Draw(">>lst","ptmu1>3&&ptmu2>3&&TMath::Abs(etamu1)<2.5&&TMath::Abs(etamu2)<2.5");
    tgen->Draw(">>lst","");
    TEventList *lst;
    lst = (TEventList*)gDirectory->Get("lst");
    tgen->SetEventList(lst);
    if(fVerbose>0) cout << " After cuts: " << lst->GetN() << " entries" << endl;

    // Get TTree with iRecoEvents
    TTree* treco = (TTree*) f->Get("events");
    if(fVerbose>0) cout << "Got TTree with " << treco->GetEntries() << " entries" << endl;

    // set branch addresses
    int genrun, genls, genevt;
    tgen->SetBranchAddress("run",&genrun);
    tgen->SetBranchAddress("LS",&genls);
    tgen->SetBranchAddress("event",&genevt);

    int recorun, recols, recoevt;
    treco->SetBranchAddress("run",&recorun);
    treco->SetBranchAddress("LS",&recols);
    treco->SetBranchAddress("event",&recoevt);

    int mcmatch;
    treco->SetBranchAddress("isMCmatch",&mcmatch);

    double ptmu1, ptmu2, ptpr, ptpi;
    tgen->SetBranchAddress("ptmu1",&ptmu1);
    tgen->SetBranchAddress("ptmu2",&ptmu2);
    tgen->SetBranchAddress("ptpr",&ptpr);
    tgen->SetBranchAddress("ptpi",&ptpi);
    double etamu1, etamu2, etapr, etapi;
    tgen->SetBranchAddress("etamu1",&etamu1);
    tgen->SetBranchAddress("etamu2",&etamu2);
    tgen->SetBranchAddress("etapr",&etapr);
    tgen->SetBranchAddress("etapi",&etapi);
    double phimu1, phimu2, phipr, phipi;
    tgen->SetBranchAddress("phimu1",&phimu1);
    tgen->SetBranchAddress("phimu2",&phimu2);
    tgen->SetBranchAddress("phipr",&phipr);
    tgen->SetBranchAddress("phipi",&phipi);

    cout << "#      run         LS      event     genidx    recoidx" << endl;
    cout << "---------- ---------- ---------- ---------- ----------" << endl;

    TLorentzVector tlvmu1, tlvmu2, tlvpr, tlvpi;
    tmph1 = new TH1F ("tmph1","tmphisto",20,3.096,3.098);
    tmph2 = new TH1F ("tmph2","tmphisto",20,1.112,1.120);
    tmph3 = new TH1F ("tmph3","tmphisto",20,5.61,5.64);
    tmph4 = new TH2F ("tmph4","tmphisto",20,0,5,20,-2,8);

    // outer loop on gen events
    for (int i=0; i!=lst->GetN(); i++)
    {
	tgen->GetEntry(lst->GetEntry(i));
	tlvmu1.SetPtEtaPhiM(ptmu1,etamu1,phimu1,0.105658367);
	tlvmu2.SetPtEtaPhiM(ptmu2,etamu2,phimu2,0.105658367);
	tlvpr.SetPtEtaPhiM( ptpr, etapr, phipr ,0.938272013);
	tlvpi.SetPtEtaPhiM( ptpi, etapi, phipi ,0.13957018);
	tmph1->Fill((tlvmu1+tlvmu2).M());
	tmph2->Fill((tlvpr+tlvpi).M());
	tmph3->Fill((tlvmu1+tlvmu2+tlvpr+tlvpi).M());
	tmph4->Fill((tlvpr+tlvpi).P(),tlvpr.P()-tlvpi.P());
	//cout << (tlvmu1+tlvmu2).M() << " " << (tlvpr+tlvpi).M() << " " << (tlvmu1+tlvmu2+tlvpr+tlvpi).M() << endl;
	if (!matchesonly)
	    cout << setw(10) << genrun << " " << setw(10) << genls << " " << setw(10) << genevt << " " << setw(10) << lst->GetEntry(i) << " ";
	// inner loop on reco events
	int nfound(0);
	for (int j=0; j!=treco->GetEntries(); j++)
	{
	    treco->GetEntry(j);
	    if(genrun==recorun && genls==recols && genevt==recoevt)
	    {
		nfound++;
		if (matchesonly && nfound==1)
		    cout << setw(10) << genrun << " " << setw(10) << genls << " " << setw(10) << genevt << " " << setw(10) << lst->GetEntry(i) << " ";
		if(nfound>1) cout << "                                 " << setw(10) << lst->GetEntry(i) << " ";
		cout << setw(10) << j << " match: " << mcmatch << endl;
	    }
	}
	if (!matchesonly && nfound==0) cout << "         -" << endl;
    }
    c2 = new TCanvas();
    c2->Divide(2,2);
    c2->cd(1);
    tmph1->Draw();
    c2->cd(2);
    tmph2->Draw();
    c2->cd(3);
    tmph3->Draw();
    c2->cd(4);
    tmph4->Draw();
}

