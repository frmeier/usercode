#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TRandom3.h"

#include "Cut.h"
#include "Cuts.h"

using std::string;
using std::cout;
using std::endl;

void mkDecTree(string infile, string outfile, string obj = "lb", string addCut = "")
{
    cout << "Infile: " << infile << endl;
    cout << "Outfile: " << outfile << endl;

    // source tree
    TFile* tfIn = TFile::Open(infile.c_str(), "READ");
    if (tfIn == 0)
    {
	cout << "Problem opening infile \"" << infile << "\" - exting";
	return;
    }
    TTree* intree = (TTree*) tfIn->Get("events");
    // cut selection
    Cuts cut;
    //cut.selectCut("isMCmatch");
    if (obj == "lb") cut.selectCut("acc04Lb","lb07");
    if (obj == "B0") cut.selectCut("acc04B0","B001");
    //if (obj == "B0") cut.selectCut("acc04B0");
    //if (obj == "B0") cut.selectCut("acc03B0","B001exp");

    if (addCut.size() > 0)
	if (cut.getCut().size() > 0)
	    intree->Draw(">>lst", (cut.getCut()+"&&"+addCut).c_str());
	else
	    intree->Draw(">>lst", addCut.c_str());
    else
	intree->Draw(">>lst", cut.getCut().c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");

    double inmass, intau, intauE, intauTruth;
    intree->SetBranchAddress(("m"+obj).c_str(),&inmass);
    intree->SetBranchAddress(("ct3d"+obj).c_str(),&intau);
    intree->SetBranchAddress(("ct3d"+obj+"E").c_str(),&intauE);
    intree->SetBranchAddress(("ct"+obj+"truth").c_str(),&intauTruth);

    // destination tree
    TFile* tfOut = new TFile(outfile.c_str(),"RECREATE");
    if (tfOut == 0)
    {
	cout << "Problem opening outfile \"" << outfile << "\" - exting";
	return;
    }

    double outmass, outtau, outtauE, outtauTruth, outtauDiff;
    TTree * outtree = new TTree("fittree","Tree with very reduced selection");
    outtree->Branch("mass", &outmass, "mass/D");
    outtree->Branch("t",    &outtau,  "t/D");
    outtree->Branch("tE",   &outtauE, "tE/D");
    outtree->Branch("tTruth",&outtauTruth,  "tTruth/D");
    outtree->Branch("tDiff",&outtauDiff,  "tDiff/D");

    cout << "Selection has " << lst->GetN() << " entries" << endl;
    for(Long64_t i = 0; i!=lst->GetN(); i++)
    {
	intree->GetEntry(lst->GetEntry(i));
	outmass = inmass;
	outtau = intau;
	outtauE = intauE;
	outtauTruth = intauTruth;
	outtauDiff = intau - intauTruth;
	outtree->Fill();
    }
    tfOut->Write();
    tfOut->Close();
    tfIn->Close();
    delete tfIn;
    delete tfOut;
}

void mkRandomTree(string outfile, Long64_t N, double bgrFrac = .5, UInt_t seed = 4357)
{
    // destination tree
    TFile* tfOut = new TFile(outfile.c_str(),"RECREATE");
    if (tfOut == 0)
    {
	cout << "Problem opening outfile \"" << outfile << "\" - exting";
	return;
    }

    double outmass, outtau, outtauE, outtauTruth;
    TTree * outtree = new TTree("fittree","Tree with very reduced selection");
    outtree->Branch("mass", &outmass, "mass/D");
    outtree->Branch("t",    &outtau,  "t/D");
    outtree->Branch("tE",   &outtauE, "tE/D");
    outtree->Branch("tTruth",&outtauTruth,  "tTruth/D");

    TRandom3 *ran = new TRandom3(seed);

    //const double mass_truth(5.27950); // B0
    //const double width_truth(0.015);
    //const double tau_truth(1.525e-12);
    const double mass_truth(5.6202); // Lb
    const double mass_window(0.2);
    const double width_truth(0.016);
    const double tau_truth(1.391e-12);
    const double taureso_truth(0.8e-12);
    const double tauresoE_truth(.1e-12);
    //const double tau_bgr_truth(3e-12);
    const double tau_bgr_truth(tau_truth);
    for(Long64_t i = 0; i!=N; i++)
    {
	if (ran->Uniform(1) > bgrFrac)
	{   // we do a signal events
	    outmass = ran->Gaus(mass_truth,width_truth);
	    outtauTruth = ran->Exp(tau_truth);
	    outtau = ran->Gaus(outtauTruth,taureso_truth);
	    outtauE = ran->Gaus(taureso_truth,tauresoE_truth);
	    outtree->Fill();
	}
	else
	{   // we do a background event
	    outmass = mass_truth-mass_window+ran->Uniform(2*mass_window);
	    outtauTruth = ran->Exp(tau_truth);
	    outtau = ran->Gaus(outtauTruth,taureso_truth);
	    outtauE = ran->Gaus(taureso_truth,tauresoE_truth);
	    outtree->Fill();
	}
    }
    tfOut->Write();
    tfOut->Close();
    delete tfOut;
}

