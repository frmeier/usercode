#include <string>
#include <iostream>

#include "TROOT.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TPad.h"
#include "TCanvas.h"

#include "utils.h"
#include "setTDRStyle_modified.C"

using std::cout;
using std::endl;
using std::string;

// case/I:amount/D:t_lbbar:tE_lbbar:m_lbbar:mE_lbbar:N_lbbar:t_lb:tE_lb:m_lb:mE_lb:N_lb:t_lbboth:tE_lbboth:m_lbboth:mE_lbboth:N_lbboth:t_b0:tE_b0:m_b0:mE_b0:N_b0

string getCaseName(int i)
{
    switch(i)
    {
	case 0: return "std"; break;
	case 1: return "#Deltar/r"; break;
	case 2: return "#Deltaz/r"; break;
	case 3: return "r#Delta#varphi/r"; break;
	case 4: return "#Deltar/z"; break;
	case 5: return "#Deltaz/z"; break;
	case 6: return "r#Delta#varphi/z"; break;
	case 7: return "#Deltar/r#Delta#varphi"; break;
	case 8: return "#Deltaz/r#Delta#varphi"; break;
	case 9: return "r#Delta#varphi/r#Delta#varphi"; break;
	case 10: return "surf"; break;
	case 11: return "layer1"; break;
    }
    return "N/A";
}

string getCaseName(int i, double val)
{
    if (val>0) return ("+ "+getCaseName(i));
    if (val<0) return ("- "+getCaseName(i));
    return getCaseName(i);
}

void draw1D(TTree *t, string hname, string valuename, string errorname, string valuetitle, string unit)
{
    cout << "hname" << endl;
    int type_val;
    double value, error, amount_val;
    t->SetBranchAddress("case", &type_val);
    t->SetBranchAddress("amount", &amount_val);
    t->SetBranchAddress(valuename.c_str(), &value);
    t->SetBranchAddress(errorname.c_str(), &error);
    const int N = t->GetEntries();
    TH1D *h = new TH1D(hname.c_str(), "",  N, 0.0, (double)N);
    for (int i=0; i!=t->GetEntries(); i++)
    {
	t->GetEntry(i);
	if(value!=0)
	{
	    h->SetBinContent(i+1, value);
	    h->SetBinError(i+1, error);
	}
	cout << value << " " << error << endl;
	h->GetXaxis()->SetBinLabel(i+1, getCaseName(type_val, amount_val).c_str());
	h->GetXaxis()->LabelsOption("v");
    }
    h->Draw();
    h->GetYaxis()->SetTitle(valueWithUnit(valuetitle, unit).c_str());
    gPad->SetBottomMargin(.15);
    gPad->SaveAs(("ali_"+hname+".pdf").c_str());
}

void draw1D(TTree *t, string hname, string valuename, string valuetitle, string unit)
{
    cout << "hname" << endl;
    int type_val;
    double value, error, amount_val;
    t->SetBranchAddress("case", &type_val);
    t->SetBranchAddress("amount", &amount_val);
    t->SetBranchAddress(valuename.c_str(), &value);
    const int N = t->GetEntries();
    TH1D *h = new TH1D(hname.c_str(), "",  N, 0.0, (double)N);
    for (int i=0; i!=t->GetEntries(); i++)
    {
	t->GetEntry(i);
	if(value!=0)
	{
	    h->SetBinContent(i+1, value);
	    h->SetBinError(i+1, sqrt(value));
	}
	cout << value << " " << error << endl;
	h->GetXaxis()->SetBinLabel(i+1, getCaseName(type_val, amount_val).c_str());
	h->GetXaxis()->LabelsOption("v");
    }
    h->SetMinimum(0);
    h->Draw();
    h->GetYaxis()->SetTitle(valueWithUnit(valuetitle, unit).c_str());
    gPad->SetBottomMargin(.15);
    gPad->SaveAs(("ali_"+hname+".pdf").c_str());
}

void doAliStudies(int i)
{
    setTDRStyle();
    TTree *t = new TTree();
    t->ReadFile("../data/aliSyst2.dat");
    TCanvas *c = new TCanvas("c","c",800,400);
    if (i==1) draw1D(t, "t_lbboth", "t_lbboth", "tE_lbboth", "t(#Lambda_{b}+#bar{#Lambda}_{b})", "s");
    if (i==2) draw1D(t, "m_lbboth", "m_lbboth", "mE_lbboth", "m(#Lambda_{b}+#bar{#Lambda}_{b})", "GeV/c^{2}");
    if (i==3) draw1D(t, "N_lbboth", "N_lbboth", "N(#Lambda_{b}+#bar{#Lambda}_{b})", "");

    if (i==4) draw1D(t, "t_lb", "t_lb", "tE_lb", "t(#Lambda_{b})", "s");
    if (i==5) draw1D(t, "m_lb", "m_lb", "mE_lb", "m(#Lambda_{b})", "GeV/c^{2}");
    if (i==6) draw1D(t, "N_lb", "N_lb", "N(#Lambda_{b})", "");

    if (i==7) draw1D(t, "t_lbbar", "t_lbbar", "tE_lbbar", "t(#bar{#Lambda}_{b})", "s");
    if (i==8) draw1D(t, "m_lbbar", "m_lbbar", "mE_lbbar", "m(#bar{#Lambda}_{b})", "GeV/c^{2}");
    if (i==9) draw1D(t, "N_lbbar", "N_lbbar", "N(#bar{#Lambda}_{b})", "");

    if (i==10) draw1D(t, "t_b0", "t_b0", "tE_b0", "t(B^{0})", "s");
    if (i==11) draw1D(t, "m_b0", "m_b0", "mE_b0", "m(B^{0})", "GeV/c^{2}");
    if (i==12) draw1D(t, "N_b0", "N_b0", "N(B^{0})", "");
}

