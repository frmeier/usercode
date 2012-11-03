#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TRandom3.h"

#include "utils.h"

using std::string;
using std::cout;
using std::endl;

const double timeFactor(1e12); // choose 1e12 if you want to have all data in ps, choose 1 if s

//============================================================================================
void mkRandomTree(string outfile, Long64_t N, double bgrFrac = .5, double nonpromptFrac = .5, UInt_t seed = 4357)
{
    // destination tree
    TFile* tfOut = new TFile(outfile.c_str(),"RECREATE");
    if (tfOut == 0)
    {
	cout << "Problem opening outfile \"" << outfile << "\" - exting";
	return;
    }

    double outmass, outtau, outtauRed, outtauE, outtauTruth, outtauDiff, outp, outd3sig;
    TTree * outtree = new TTree("fittree","Tree with very reduced selection");
    outtree->Branch("mass", &outmass, "mass/D");
    outtree->Branch("t",    &outtau,  "t/D");
    outtree->Branch("tRed", &outtauRed,"tRed/D");
    outtree->Branch("tE",   &outtauE, "tE/D");
    outtree->Branch("tTruth",&outtauTruth,  "tTruth/D");
    outtree->Branch("tDiff",&outtauDiff,  "tDiff/D");
    outtree->Branch("p",    &outp, "p/D");
    outtree->Branch("d3sig",&outd3sig, "d3sig/D");

    TRandom3 *ran = new TRandom3(seed);

    bool isB0(false);

    // signal
    const double mass_truth(isB0 ? 5.27950 : 5.6202);
    const double width_truth(0.015);

    const double tau_truth(isB0 ? 1.536e-12*timeFactor : 1.507e-12*timeFactor);
    const double taureso_truth(0.2e-12*timeFactor);
    const double tauresoE_truth(.04e-12*timeFactor);

    const double reducedSigma(3);

    // background
    const double mass_lo(isB0 ? 5.16 : 5.4);
    const double mass_hi(isB0 ? 5.75 : 6.0);
    const double mass_window(mass_hi-mass_lo);

    const double bgr_prompt_width(0.2e-12);
    const double bgr_nonprompt_tau1(1.15e-12*timeFactor);
    const double bgr_nonprompt_tau2(0.75e-12*timeFactor);
    const double bgr_nonprompt_frac1(.99);
    const double bgr_taureso(0.2e-12*timeFactor);

    // now dice the tree
    for(Long64_t i = 0; i!=N; i++)
    {
	if (ran->Uniform(1) > bgrFrac)
	{   // we do a signal events
	    outmass = ran->Gaus(mass_truth,width_truth);
	    outtauTruth = ran->Exp(tau_truth);
	    outtau = ran->Gaus(outtauTruth,taureso_truth);
	    outtauE = ran->Gaus(taureso_truth,tauresoE_truth);
	    outtauRed = outtau-reducedSigma*outtauE;
	    outtauDiff = outtau - outtauTruth;
	    outp = -9999; // just there to make trees mergable
	    outd3sig = (outtauE!=0) ? outtauE/outtau : -9999.0; // not really compatible with what we have in data
	    outtree->Fill();
	}
	else
	{   // we do a background event
	    outmass = mass_lo+ran->Uniform(mass_window);
	    if (ran->Uniform(1) > nonpromptFrac)
	    {	// we make a prompt bgr event
		outtauTruth = 0.0;
		outtau = ran->Gaus(0.0, bgr_prompt_width);
		outtauE = ran->Gaus(bgr_taureso, tauresoE_truth);
	    }
	    else
	    {   // we make a non-prompt bgr event
		if (ran->Uniform(1) > bgr_nonprompt_frac1)
		    outtauTruth = ran->Exp(bgr_nonprompt_tau1);
		else
		    outtauTruth = ran->Exp(bgr_nonprompt_tau2);
		outtau = ran->Gaus(outtauTruth, bgr_taureso);
		outtauE = ran->Gaus(bgr_taureso, tauresoE_truth);
	    }
	    outtauRed = outtau-reducedSigma*outtauE;
	    outtauDiff = outtau - outtauTruth;
	    outp = -9999; // just there to make trees mergable
	    outd3sig = (outtauE!=0) ? outtauE/outtau : -9999.0; // not really compatible with what we have in data
	    outtree->Fill();
	}
    }
    tfOut->Write();
    tfOut->Close();
    delete tfOut;
}

void mkSomeRandomTrees(int nFiles, string outfileStem, Long64_t N, double bgrFrac = .5, double nonpromptFrac = .5, UInt_t seed = 4357)
{
    for (int i = 0; i != nFiles; i++)
    {
	string filename = outfileStem + toString(i) + ".root";
	cout << "Working on file " << filename << endl;
	UInt_t curSeed = seed+i;
	mkRandomTree(filename, N, bgrFrac, nonpromptFrac, curSeed);
    }
}

