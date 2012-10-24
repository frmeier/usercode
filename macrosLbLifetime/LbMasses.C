#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TText.h"
#include <string>
#include <vector>
#include "setTDRStyle_modified.C"

using std::string;
using std::vector;

struct Msrmt
{
    Msrmt(double t, double stMi, double stMa, double syMi, double syMa, string year, string exp, string ch) :
	t_(t), statMin_(stMi), statMax_(stMa), systMin_(syMi), systMax_(syMa), year_(year), experiment_(exp), channel_(ch) {};
    double t_;
    double statMin_, statMax_;
    double systMin_, systMax_;
    string year_, experiment_, channel_;
};

void LbMasses()
{
    // Data
    vector<Msrmt> msrmts;
    //                     val   --stat--  --syst--   yrs    experiment  channel
    msrmts.push_back(Msrmt(1.110,0.18,0.19,0.05,0.05,"91-94","DELPHI","#Lambda_{c}^{+}#font[12]{l}"));
    msrmts.push_back(Msrmt(1.180,0.12,0.13,0.03,0.03,"91-95","ALEPH","#Lambda_{c}^{+}#font[12]{l}"));
    msrmts.push_back(Msrmt(1.300,0.21,0.26,0.04,0.04,"91-95","ALEPH","#Lambda#font[12]{l}^{-}#font[12]{l}^{+}"));
    msrmts.push_back(Msrmt(1.320,0.15,0.15,0.07,0.07,"91-95","CDF1","#Lambda_{c}^{+}#font[12]{l}"));
    msrmts.push_back(Msrmt(1.290,0.22,0.24,0.06,0.06,"90-95","OPAL","#Lambda_{c}^{+}#font[12]{l}, #Lambda#font[12]{l}^{-}#font[12]{l}^{+}"));
    msrmts.push_back(Msrmt(1.401,.046,.046,.035,.035,"02-06","CDF2","#Lambda_{c}^{+}#pi"));
    msrmts.push_back(Msrmt(1.290,.110,.119,.091,.087,"02-06","D0","#Lambda_{c}^{+}#mu"));
    msrmts.push_back(Msrmt(1.218,.115,.130,.042,.042,"02-06","D0","J/#psi#Lambda"));
    msrmts.push_back(Msrmt(1.537,.045,.045,.014,.014,"02-09","CDF2","J/#psi#Lambda"));
    msrmts.push_back(Msrmt(1.303,.075,.075,.035,.035,"02-11","D0","J/#psi#Lambda"));
    msrmts.push_back(Msrmt(1.449,.036,.036,.017,.017,"2011","ATLAS","J/#psi#Lambda"));
    msrmts.push_back(Msrmt(1.506,.053,.053,.030,.030,"2011","CMS prel","J/#psi#Lambda"));
    //msrmts.push_back(Msrmt());

    const int msrmtSize = msrmts.size();
    double y[msrmtSize];
    double t[msrmtSize], tmin[msrmtSize], tmax[msrmtSize];
    double tmin2[msrmtSize], tmax2[msrmtSize]; // for adding stat+syst
    for(int i = 0; i!=msrmtSize; i++)
    {
	y[i] = (double)i+1;
	t[i] = msrmts[i].t_;
	tmin[i] = msrmts[i].statMin_;
	tmax[i] = msrmts[i].statMax_;
	tmin2[i] = sqrt(msrmts[i].statMin_*msrmts[i].statMin_+msrmts[i].systMin_*msrmts[i].systMin_);
	tmax2[i] = sqrt(msrmts[i].statMax_*msrmts[i].statMax_+msrmts[i].systMax_*msrmts[i].systMax_);
    }

    // prepare canvas
    setTDRStyle();
    TCanvas *c = new TCanvas("c","c",700,600);
    c->SetRightMargin(.42);
    c->SetLeftMargin(.02);
    c->SetTopMargin(.10);

    // Draw plot
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(msrmtSize,t,y,tmin,tmax,0,0);
    gr->SetTitle("#Lambda_{b} lifetime");
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleSize(0.06);

    gr->GetYaxis()->SetBinLabel(1,"");
    gr->GetXaxis()->SetTitle("#tau [ps]");
    gr->GetXaxis()->SetNdivisions(50205);

    // write the labels
    TLatex tl;
    tl.SetTextAlign(11);
    tl.SetTextSize(0.04);
    gr->SetMarkerStyle(21);
    gr->Draw("AP");
    const double xpos = 1.7;
    for(int i = 0; i!=msrmtSize; i++)
    {
	const double ypos = 1+i*1;
	tl.DrawLatex(xpos,ypos,msrmts[i].experiment_.c_str());
	tl.DrawLatex(xpos+0.17,ypos,("("+msrmts[i].year_+")").c_str());
	tl.DrawLatex(xpos+0.34,ypos,msrmts[i].channel_.c_str());
    }
    const double yspace = 0.30;
    const double ypos1 = 0.30, ypos2 = -0.90;
    tl.DrawLatex(xpos, ypos1-0*yspace,"#scale[.6]{errors in black: statistical only}");
    tl.DrawLatex(xpos, ypos1-1*yspace,"#scale[.6]{errors in #color[14]{grey}: syst. added in quadrature}");
    tl.DrawLatex(xpos, ypos1-2*yspace,"#scale[.6]{#color[592]{band}: current best value}");
    tl.DrawLatex(xpos, ypos2-0*yspace,"#scale[.6]{data from arXiv:1010.1589}");
    tl.DrawLatex(xpos, ypos2-1*yspace,"#scale[.6]{prepared for PDG2011}");

    // Draw plot with syst errors
    TGraphAsymmErrors *grs = new TGraphAsymmErrors(msrmtSize,t,y,tmin2,tmax2,0,0);
    grs->SetTitle("#Lambda_b lifetime (data from arXiv:1010.1589, for PDG2011)");
    grs->SetLineColor(15);
    grs->GetYaxis()->SetBinLabel(1,"");
    grs->GetXaxis()->SetTitle("#tau [ps]");

    // Draw box for PDG value
    double pdg_mean(1.425), pdg_sigma(0.032);
    double pdg_lo(pdg_mean-pdg_sigma), pdg_hi(pdg_mean+pdg_sigma);
    double pdg_y((double)msrmtSize/2+.5), pdgye((double)msrmtSize/2);
    TGraphAsymmErrors *pdg_gr = new TGraphAsymmErrors(1,&pdg_mean,&pdg_y,&pdg_sigma,&pdg_sigma,&pdgye,&pdgye);
    pdg_gr->SetFillColor(590);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(pdg_gr,"2");
    mg->Add(grs,"P");
    mg->Add(gr,"AP");
    mg->Draw();

}


