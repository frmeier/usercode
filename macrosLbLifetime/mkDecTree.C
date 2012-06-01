#include <string>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TRandom3.h"

#include "Cuts.C"
#include "utils.h"

using std::string;
using std::cout;
using std::endl;

void mkDecTree(string infile, string outfile, string obj = "lb", string addCut = "", string addCut1 = "HLT_jpsiBarrel", string addCut2 = "HLT_matched", double reducedSigma = 3.0)
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
    //if (obj == "B0") cut.selectCut("acc05B0", "muSoft", "B004", "HLT_jpsiBarrel", "HLT_matched");
    //if (obj == "lb") cut.selectCut("acc05Lb", "muSoft", "lb11", addCut1, addCut2);
    //if (obj == "B0") cut.selectCut("acc05B0", "muSoft", "B005", addCut1, addCut2);
    if (obj == "lb") cut.selectCut("acc06Lb", "muSoft", "lb12", addCut1, addCut2);
    if (obj == "B0") cut.selectCut("acc06B0", "muSoft", "B006", addCut1, addCut2);

    if (addCut.size() > 0)
	if (cut.getCut().size() > 0)
	    intree->Draw(">>lst", (cut.getCut()+"&&"+addCut).c_str());
	else
	    intree->Draw(">>lst", addCut.c_str());
    else
	intree->Draw(">>lst", cut.getCut().c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");

    obj = "bc"; // Im Moment och so geloest damit allenfalls alte reduced trees noch einfach lesbar sind. TODO: obj-String unten fallenlassen und durch bc ersetzen

    double inmass, inp, intau, intauE, intauTruth, ind3, ind3E;
    intree->SetBranchAddress(("m"+obj).c_str(),&inmass);
    intree->SetBranchAddress(("p"+obj).c_str(),&inp);
    intree->SetBranchAddress(("ct3d"+obj).c_str(),&intau);
    intree->SetBranchAddress(("ct3d"+obj+"E").c_str(),&intauE);
    intree->SetBranchAddress(("ct"+obj+"truth").c_str(),&intauTruth);
    intree->SetBranchAddress(("d3"+obj).c_str(),&ind3);
    intree->SetBranchAddress(("d3E"+obj).c_str(),&ind3E);

    // destination tree
    TFile* tfOut = new TFile(outfile.c_str(),"RECREATE");
    if (tfOut == 0)
    {
	cout << "Problem opening outfile \"" << outfile << "\" - exting";
	return;
    }

    double outmass, outp, outtau, outtauRed, outtauE, outtauTruth, outtauDiff, outd3sig;
    TTree * outtree = new TTree("fittree","Tree with very reduced selection");
    outtree->Branch("mass", &outmass, "mass/D");
    outtree->Branch("t",    &outtau,  "t/D");
    outtree->Branch("tRed", &outtauRed,"tRed/D");
    outtree->Branch("tE",   &outtauE, "tE/D");
    outtree->Branch("tTruth",&outtauTruth,  "tTruth/D");
    outtree->Branch("tDiff",&outtauDiff,  "tDiff/D");
    outtree->Branch("p",    &outp, "p/D");
    outtree->Branch("d3sig",&outd3sig, "d3sig/D");

    cout << "Selection has " << lst->GetN() << " entries" << endl;
    for(Long64_t i = 0; i!=lst->GetN(); i++)
    {
	intree->GetEntry(lst->GetEntry(i));
	outmass = inmass;
	outp = inp;
	outtau = intau;
	outtauRed = intau-reducedSigma*intauE;
	outtauE = intauE;
	outtauTruth = intauTruth;
	outtauDiff = intau - intauTruth;
	outd3sig = (ind3E!=0) ? ind3/ind3E : -9999.0;
	outtree->Fill();
    }
    tfOut->Write();
    tfOut->Close();
    tfIn->Close();
    delete tfIn;
    delete tfOut;
}

void mkSomeDecTrees()
{
    /*
    mkDecTree("../data/run380__run384.root", "../data/vrt_r380__384_Lb_data_lb11_barrel.root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run380__run384.root", "../data/vrt_r380__384_Lb_data_lb11_displ.root", "lb", "", "HLT_jpsiDispl", "HLT_matched", 3.0);
    mkDecTree("../data/run385__run387.root", "../data/vrt_r385__387_Lb_sigMC_Lb_lb11_displ.root", "lb", "", "HLT_jpsiDispl", "HLT_matched", 3.0);
    mkDecTree("../data/run385__run387.root", "../data/vrt_r385__387_Lb_sigMC_Lb_lb11_barrel.root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run388.root", "../data/vrt_r388_Lb_bgrMC_Jp_lb11_notrigsel.root", "lb", "", "", "", 3.0);
    mkDecTree("../data/run389.root", "../data/vrt_r389_Lb_bgrMC_Bs_lb11_notrigsel.root", "lb", "", "", "", 3.0);
    mkDecTree("../data/run390.root", "../data/vrt_r390_Lb_bgrMC_Bp_lb11_notrigsel.root", "lb", "", "", "", 3.0);
    mkDecTree("../data/run391.root", "../data/vrt_r391_Lb_bgrMC_B0_lb11_notrigsel.root", "lb", "", "", "", 3.0);
    mkDecTree("../data/run392.root", "../data/vrt_r392_Lb_offMC_Lb_lb11_notrigsel.root", "lb", "", "", "", 3.0);

    mkDecTree("../data/run393__run397.root", "../data/vrt_r393__397_B0_data_B005_barrel.root", "B0", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run393__run397.root", "../data/vrt_r393__397_B0_data_B005_displ.root", "B0", "", "HLT_jpsiDispl", "HLT_matched", 3.0);
    mkDecTree("../data/run398__run399.root", "../data/vrt_r398__399_B0_sigMC_B0_B005_barrel.root", "B0", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run398__run399.root", "../data/vrt_r398__399_B0_sigMC_B0_B005_displ.root", "B0", "", "HLT_jpsiDispl", "HLT_matched", 3.0);
    mkDecTree("../data/run400.root", "../data/vrt_r400_B0_bgrMC_Jp_B005_notrigsel.root", "B0", "", "", "", 3.0);
    mkDecTree("../data/run401.root", "../data/vrt_r401_B0_bgrMC_Bs_B005_notrigsel.root", "B0", "", "", "", 3.0);
    mkDecTree("../data/run402.root", "../data/vrt_r402_B0_bgrMC_Bp_B005_notrigsel.root", "B0", "", "", "", 3.0);
    mkDecTree("../data/run403.root", "../data/vrt_r403_B0_bgrMC_Lb_B005_notrigsel.root", "B0", "", "", "", 3.0);
    mkDecTree("../data/run404.root", "../data/vrt_r404_B0_offMC_B0_B005_notrigsel.root", "B0", "", "", "", 3.0);
    */

    mkDecTree("../data/run460.root", "../data/vrt_r460_B0_sigMC_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run460.root", "../data/vrt_r460_B0_sigMC_B0_B006_displ.root", "B0", "", "HLT_jpsiDispl", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_displ.root", "B0", "", "HLT_jpsiDispl", "HLT_matched", 3.0);

    mkDecTree("../data/run456.root", "../data/vrt_r456_Lb_sigMC_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDispl", "HLT_matched", 3.0);
    mkDecTree("../data/run456.root", "../data/vrt_r456_Lb_sigMC_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDispl", "HLT_matched", 3.0);

    // some crosscheck samples
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_barrel_phiplus.root", "B0", "phibc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_barrel_phiminus.root", "B0", "phibc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_barrel_etaplus.root", "B0", "etabc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_barrel_etaminus.root", "B0", "etabc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_barrel_runA.root", "B0", "run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_barrel_runB.root", "B0", "run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);

    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel_phiplus.root", "lb", "phibc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel_phiminus.root", "lb", "phibc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel_etaplus.root", "lb", "etabc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel_etaminus.root", "lb", "etabc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel_runA.root", "lb", "run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel_runB.root", "lb", "run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel_lb.root", "lb", "rqha1>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel_lbbar.root", "lb", "rqha1<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
}

void mkRandomTree(string outfile, Long64_t N, double bgrFrac = .5, double nonpromptFrac = .5, UInt_t seed = 4357)
{
    // destination tree
    TFile* tfOut = new TFile(outfile.c_str(),"RECREATE");
    if (tfOut == 0)
    {
	cout << "Problem opening outfile \"" << outfile << "\" - exting";
	return;
    }

    double outmass, outtau, outtauRed, outtauE, outtauTruth;
    TTree * outtree = new TTree("fittree","Tree with very reduced selection");
    outtree->Branch("mass", &outmass, "mass/D");
    outtree->Branch("t",    &outtau,  "t/D");
    outtree->Branch("tRed", &outtauRed,"tRed/D");
    outtree->Branch("tE",   &outtauE, "tE/D");
    outtree->Branch("tTruth",&outtauTruth,  "tTruth/D");

    TRandom3 *ran = new TRandom3(seed);

    bool isB0(true);

    // signal
    const double mass_truth(isB0 ? 5.27950 : 5.6202);
    const double width_truth(0.015);

    const double tau_truth(isB0 ? 1.536e-12 : 1.391e-12);
    const double taureso_truth(0.2e-12);
    const double tauresoE_truth(.04e-12);

    const double reducedSigma(3);

    // background
    const double mass_window(0.5);

    const double bgr_prompt_width(0.2e-12);
    const double bgr_nonprompt_tau1(1.1e-12);
    const double bgr_nonprompt_tau2(0.75e-12);
    const double bgr_nonprompt_frac1(.3);
    const double bgr_taureso(0.2e-12);

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
	    outtree->Fill();
	}
	else
	{   // we do a background event
	    outmass = mass_truth-mass_window+ran->Uniform(2*mass_window);
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

