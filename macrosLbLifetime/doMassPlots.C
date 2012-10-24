#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <memory>

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TEventList.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include "utils.h"
#include "Cuts.C"
#include "doMassFitLb01.C"
#include "setTDRStyle_modified.C"

void doMassPlot(string filename, string title)
{
    setTDRStyle();
    TFile *f = TFile::Open(filename.c_str()); // Data
    TTree *t= (TTree*)f->Get("fittree");

    doMassFitLb01_fitresults res = doMassFitLb01_DG(t, 0, title, false, true, false, false, "mass");
    cout << "Results: sig: " << res.sig << " bgr: " << res.bgr << endl;
}

void doSomeResoPlots()
{
    const string path = "../data/";
    std::vector< std::pair<string,string> > filelist;

    /*
    const int n_bin_ptha1 = 5;
    double bin_ptha1[n_bin_ptha1]   = { 2.0, 1.9, 1.8, 1.7, 1.6 };
    for (int i=0; i!=n_bin_ptha1; i++)
    {
	filelist.push_back(std::make_pair("#Lambda_{b}+#bar{#Lambda}_{b} pt>"+toString(bin_ptha1[i]), path+"vrt_r480_data_lb_acc_barrel_ptpr_lball_"+toString(i)+".root"));
	filelist.push_back(std::make_pair("#Lambda_{b} pt>"+toString(bin_ptha1[i]), path+"vrt_r480_data_lb_acc_barrel_ptpr_lb_"+toString(i)+".root"));
	filelist.push_back(std::make_pair("#bar{#Lambda}_{b} pt>"+toString(bin_ptha1[i]), path+"vrt_r480_data_lb_acc_barrel_ptpr_lbbar_"+toString(i)+".root"));
    }
    */

    for (int i=0; i!=4; i++)
	for (int j=0; j!=4; j++)
	    for (int k=0; k!=4; k++)
	    {
		filelist.push_back(std::make_pair("#Lambda_{b} "+toString(i)+"_"+toString(j)+"_"+toString(k), path+"vrt_r480_data_lb_acc_barrel_cutbin_"+toString(i)+"_"+toString(j)+"_"+toString(k)+".root"));
	    }

    //filelist.push_back(std::make_pair("", path+"vrt_r480_data_lb_acc_barrel_ptpr_lbbar_4"));

    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    const int nX(2), nY(2);
    c->Divide(nX, nY);
    int cdcounter = 1;

    for (unsigned int i = 0; i!=filelist.size(); i++)
    {
	cout << "---- " << filelist[i].first << endl;
	c->cd(cdcounter);
	doMassPlot(filelist[i].second, filelist[i].first);

	cdcounter++;
	if (cdcounter > nX*nY)
	{
	    c->SaveAs(("c_"+toString(i)+".pdf").c_str());
	    cdcounter = 1;
	}
    }
    c->SaveAs(("c_"+toString(filelist.size())+".pdf").c_str());
}

