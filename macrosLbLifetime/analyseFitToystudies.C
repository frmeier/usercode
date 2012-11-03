#include "utils.h"

string fileprefix = "fitToyStudies";

void plotHisto(TTree *t, string hname, string valtitle, string valunit, string toPlot, int nBins, double lo, double hi)
{
    string plotstring = makePlotsString(toPlot, hname, nBins, lo, hi);
    cout << "Thisplot: " << hname << endl;
    cout << plotstring << endl;
    t->Draw(plotstring.c_str());

    TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(hname.c_str());
    h->SetTitle(("Distribution of " + valtitle).c_str());
    h->GetXaxis()->SetTitle(valueWithUnit(valtitle,valunit).c_str());
    h->Fit("gaus");
    gStyle->SetOptFit(1);
    gPad->SaveAs((fileprefix+"_"+hname+".pdf").c_str());
}

void plotPulls(TTree *t, string hname, string valtitle, string toPlot, double trueVal, int nBins, double lo, double hi)
{
    string plotstring = makePlotsString("("+toPlot+"-"+toString(trueVal)+")/"+toPlot+"_err", hname, nBins, lo, hi);
    cout << plotstring << endl;
    t->Draw(plotstring.c_str());
    cout << "Thisplot: " << hname << endl;

    TH1F *h = (TH1F*)gDirectory->GetList()->FindObject(hname.c_str());
    h->SetTitle(("Pull distribution of " + valtitle).c_str());
    h->GetXaxis()->SetTitle(("pull "+valtitle).c_str());
    h->Fit("gaus");
    gStyle->SetOptFit(1);
    gPad->SaveAs((fileprefix+"_"+hname+".pdf").c_str());
}

void plotPulls(TTree *t, string hname, string valtitle, string toPlot, double trueVal, int nBins = 100, double range = 5.0)
{
    plotPulls(t, hname, valtitle, toPlot, trueVal, nBins, -range, range);
}

void analyseFitToystudiesLb(string filename, bool justTau = false)
{
    fileprefix = "fitToyStudiesLb";
    TTree *t = new TTree;
    //t->ReadFile(filename.c_str(),"tau:tau_err:mass_peak:mass_peak_err:m_sigma1:m_sigma1_err:m_sigma2:m_sigma2_err:frac_m_gauss:frac_m_gauss_err:core_frac:core_frac_err:tau_bk:tau_bk_err:n_nonprompt:n_nonprompt_err:n_prompt:n_prompt_err:n_signal:n_signal_err:prompt_p1:prompt_p1_err:prompt_sigma_core:prompt_sigma_core_err:prompt_sigma_tail:prompt_sigma_tail_err"); // alter Fit
    //t->ReadFile(filename.c_str(),"tau:tau_err:mass_peak:mass_peak_err:m_sigma1:m_sigma1_err:m_sigma2:m_sigma2_err:frac_m_gauss:frac_m_gauss_err:core_frac:core_frac_err:tau_bk:tau_bk_err:n_nonprompt:n_nonprompt_err:n_prompt:n_prompt_err:n_signal:n_signal_err:prompt_p1:prompt_p1_err:prompt_sigma_core:prompt_sigma_core_err:prompt_sigma_tail:prompt_sigma_tail_err:npr_reso_sigma:npr_reso_sigma_err:sig_reso_sigma:sig_reso_sigma_err"); // Fit ohne PEE
    //t->ReadFile(filename.c_str(),"tau:tau_err:mass_peak:mass_peak_err:m_sigma1:m_sigma1_err:m_sigma2:m_sigma2_err:frac_m_gauss:frac_m_gauss_err:tau_bk:tau_bk_err:n_nonprompt:n_nonprompt_err:n_prompt:n_prompt_err:n_signal:n_signal_err:prompt_p1:prompt_p1_err:reso_bias:reso_bias_err:reso_sigma:reso_sigma_err"); // Fit mit PEE 1 Gauss
    t->ReadFile(filename.c_str(),"tau:tau_err:mass_peak:mass_peak_err:m_sigma1:m_sigma1_err:m_sigma2:m_sigma2_err:frac_m_gauss:frac_m_gauss_err:tau_bk:tau_bk_err:n_nonprompt:n_nonprompt_err:n_prompt:n_prompt_err:n_signal:n_signal_err:prompt_p1:prompt_p1_err:reso_bias:reso_bias_err:reso_sigma:reso_sigma_err:npr_reso_bias:npr_reso_bias_err:npr_reso_sigma:npr_reso_sigma_err:sig_reso_bias:sig_reso_bias_err:sig_reso_sigma:sig_reso_sigma_err"); // Fit mit PEE 3 Gauss
    //t.Draw("tau>>htau");
    //t->Draw("(tau-1.507)/tau_err>>htau_pulls(100,-5,5)");
    if (!justTau)
    {
	const int Nevents = 5000;
	const double bgrfrac = 0.8;
	const double nonprfrac = 0.25;
	const int Nsig = Nevents * (1-bgrfrac);
	const int Npr = Nevents * bgrfrac * ( 1-nonprfrac);
	const int Nnpr = Nevents * bgrfrac * nonprfrac;
	plotPulls(t, "tau_pulls", "#tau(#Lambda_{b})", "tau", 1.507, 100, 5.0);
	plotPulls(t, "mass_pulls", "m(#Lambda_{b})", "mass_peak", 5.6202, 100, 5.0);
	plotPulls(t, "tau_bk_pulls", "#tau_{bgr}(#Lambda_{b})", "tau_bk", 0.750, 100, 5.0);
	plotPulls(t, "nsig_pulls", "n_{sig}", "n_signal", Nsig, 100, 5.0);
	plotPulls(t, "npr_pulls", "n_{prompt}", "n_prompt", Npr, 100, 10.0);
	plotPulls(t, "nnpr_pulls", "n_{nonprompt}", "n_nonprompt", Nnpr, 100, 10.0);
    }

    plotHisto(t, "tau_histo", "#tau(#Lambda_{b})", "ps", "tau", 100, 1.30, 1.70);
    if (!justTau)
    {
	plotHisto(t, "mass_histo", "m(#Lambda_{b})", "GeV", "mass_peak", 100, 5.615, 5.625);
	plotHisto(t, "tau_bk_histo", "#tau_{bgr}(#Lambda_{b})", "ps", "tau_bk", 100, 0.5, 0.9);
	plotHisto(t, "nsig_histo", "n_{sig}", "", "n_signal", 100, .8*Nsig, 1.2*Nsig);
	plotHisto(t, "npr_histo", "n_{prompt}", "", "n_prompt", 100, .8*Npr, 1.2*Npr);
	plotHisto(t, "nnpr_histo", "n_{nonprompt}", "", "n_nonprompt", 100, 0.8*Nnpr, 1.2*Nnpr);
    }
}

void analyseFitToystudiesB0(string filename, bool justTau = false)
{
    fileprefix = "fitToyStudiesB0";
    TTree *t = new TTree;
    t->ReadFile(filename.c_str(),"tau:tau_err:mass_peak:mass_peak_err:m_sigma1:m_sigma1_err:m_sigma2:m_sigma2_err:frac_m_gauss:frac_m_gauss_err:core_frac:core_frac_err:tau_bk:tau_bk_err:n_nonprompt:n_nonprompt_err:n_prompt:n_prompt_err:n_signal:n_signal_err:prompt_p1:prompt_p1_err:prompt_sigma_core:prompt_sigma_core_err:prompt_sigma_tail:prompt_sigma_tail_err");

    if (!justTau)
    {
	plotPulls(t, "tau_pulls", "#tau(B^{0})", "tau", 1.529, 100, 5.0);
	plotPulls(t, "mass_pulls", "m(B^{0})", "mass_peak", 5.27950, 100, 5.0);
	plotPulls(t, "tau_bk_pulls", "#tau_{bgr}(B^{0})", "tau_bk", 0.750, 100, 5.0);
	plotPulls(t, "nsig_pulls", "n_{sig}", "n_signal", 6750, 100, 5.0);
	plotPulls(t, "npr_pulls", "n_{prompt}", "n_prompt", 3375, 100, 10);
	plotPulls(t, "nnpr_pulls", "n_{nonprompt}", "n_nonprompt", 1125, 100, 10);
    }

    plotHisto(t, "tau_histo", "#tau(B^{0})", "ps", "tau", 100, 1.30, 1.70);
    if (!justTau)
    {
	plotHisto(t, "mass_histo", "m(B^{0})", "GeV", "mass_peak", 100, 5.2785, 5.281);
	plotHisto(t, "tau_bk_histo", "#tau_{bgr}(B^{0})", "ps", "tau_bk", 100, 0.5, 0.9);
	plotHisto(t, "nsig_histo", "n_{sig}", "", "n_signal", 100, 6500, 7000);
	plotHisto(t, "npr_histo", "n_{prompt}", "", "n_prompt", 100, 2900, 3900);
	plotHisto(t, "nnpr_histo", "n_{nonprompt}", "", "n_nonprompt", 100, 950, 1450);
    }
}


