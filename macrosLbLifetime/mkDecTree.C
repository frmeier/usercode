#include <string>
#include <iostream>
#include <vector>
#include <utility>

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

const double timeFactor(1e12); // choose 1e12 if you want to have all data in ps, choose 1 if s

void mkDecTree(string infile, string cut, string outfile, string addCut, double reducedSigma)
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

    if (addCut.size() > 0)
	if (cut.size() > 0)
	    intree->Draw(">>lst", (cut+"&&"+addCut).c_str());
	else
	    intree->Draw(">>lst", addCut.c_str());
    else
	intree->Draw(">>lst", cut.c_str());
    TEventList *lst = (TEventList*)gDirectory->Get("lst");

    double inmass, inp, intau, intauE, intauTruth, ind3, ind3E;
    intree->SetBranchAddress("mbc",&inmass);
    intree->SetBranchAddress("pbc",&inp);
    intree->SetBranchAddress("ct3dbc",&intau);
    intree->SetBranchAddress("ct3dbcE",&intauE);
    intree->SetBranchAddress("ctbctruth",&intauTruth);
    intree->SetBranchAddress("d3bc",&ind3);
    intree->SetBranchAddress("d3Ebc",&ind3E);

    /*
    int qha1;
    double rrs;
    intree->SetBranchAddress("rqha1", &qha1);
    intree->SetBranchAddress("vrrs", &rrs);
    */

    // destination tree
    TFile* tfOut = new TFile(outfile.c_str(),"RECREATE");
    if (tfOut == 0)
    {
	cout << "Problem opening outfile \"" << outfile << "\" - exting";
	return;
    }

    double outmass, outp, outtau, outtauRed, outtauE, outtauTruth, outtauDiff, outd3sig;
    //double outWeight;
    TTree * outtree = new TTree("fittree","Tree with very reduced selection");
    outtree->Branch("mass", &outmass, "mass/D");
    outtree->Branch("t",    &outtau,  "t/D");
    outtree->Branch("tRed", &outtauRed,"tRed/D");
    outtree->Branch("tE",   &outtauE, "tE/D");
    outtree->Branch("tTruth",&outtauTruth,  "tTruth/D");
    outtree->Branch("tDiff",&outtauDiff,  "tDiff/D");
    outtree->Branch("p",    &outp, "p/D");
    outtree->Branch("d3sig",&outd3sig, "d3sig/D");
    //outtree->Branch("weight",&outWeight, "weight/D");

    //double rlo, effl0, effl0bar;
    //efft->SetBranchAddress("r", &rlo);
    //efft->SetBranchAddress("l0", &effl0);

    cout << "Selection has " << lst->GetN() << " entries" << endl;
    for(Long64_t i = 0; i!=lst->GetN(); i++)
    {
	intree->GetEntry(lst->GetEntry(i));
	outmass = inmass;
	outp = inp;
	outtau = intau*timeFactor;
	outtauRed = intau*timeFactor-reducedSigma*intauE*timeFactor;
	outtauE = intauE*timeFactor;
	outtauTruth = intauTruth*timeFactor;
	outtauDiff = (intau - intauTruth)*timeFactor;
	outd3sig = (ind3E!=0) ? ind3/ind3E : -9999.0;
	/*
	outWeight = 0;
	for (int j = 0; j!= efft->GetEntries(); j++)
	{
	    efft->GetEntry(j);
	    if (rlo > rrs) break;
	    outWeight = qha1>0 ? effl0 : effl0bar;
	}
	if (outWeight > 0) 
	    outWeight = 1.0 / outWeight;
	else
	    outWeight = 0;
	    */
	outtree->Fill();
    }
    tfOut->Write();
    tfOut->Close();
    tfIn->Close();
    delete tfIn;
    delete tfOut;
}

void mkDecTree(string infile, string outfile, string obj = "lb", string addCut = "", string addCut1 = "HLT_jpsiBarrel", string addCut2 = "HLT_matched", double reducedSigma = 3.0)
{
    //TTree *efft = new TTree;
    //if (obj=="lb") efft->ReadFile("efftable_L0L0bar.dat","r/D:l0:l0bar");
    //if (obj=="B0") efft->ReadFile("efftable_Ks.dat","r/D:l0");

    if ( (obj!="lb") && (obj!="B0") )
    {
	cout << "ERROR: Value in obj invalid: " << obj << endl;
	throw;
    }

    // cut selection
    Cuts cut;
    //cut.selectCut("isMCmatch");
    //if (obj == "B0") cut.selectCut("acc05B0", "muSoft", "B004", "HLT_jpsiBarrel", "HLT_matched");
    //if (obj == "lb") cut.selectCut("acc05Lb", "muSoft", "lb11", addCut1, addCut2);
    //if (obj == "B0") cut.selectCut("acc05B0", "muSoft", "B005", addCut1, addCut2);
    //if (obj == "lb") cut.selectCut("acc06Lb", "muSoft", "lb12", addCut1, addCut2);
    //if (obj == "B0") cut.selectCut("acc06B0", "muSoft", "B006", addCut1, addCut2);
    //if (obj == "lb") cut.selectCut("acc06Lb", "muSoft", "lb13exp", addCut1, addCut2);
    if (obj == "lb") cut.selectCut("acc06Lb", "muSoft", "lb14", addCut1, addCut2);
    //if (obj == "lb") cut.selectCut( "lb13exp", addCut1, addCut2);
    if (obj == "B0") cut.selectCut("acc06B0", "muSoft", "B008", addCut1, addCut2);
    //cut.removeOneCut("rptha1");
    //cut.removeOneCut("rptha2");

    mkDecTree(infile, cut.getCut(), outfile, addCut, reducedSigma);
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

    /*
    mkDecTree("../data/run460.root", "../data/vrt_r460_B0_sigMC_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run460.root", "../data/vrt_r460_B0_sigMC_B0_B006_displ.root", "B0", "", "HLT_jpsiDispl", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run461.root", "../data/vrt_r461_B0_data_B0_B006_displ.root", "B0", "", "HLT_jpsiDispl", "HLT_matched", 3.0);

    mkDecTree("../data/run456.root", "../data/vrt_r456_Lb_sigMC_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDispl", "HLT_matched", 3.0);
    mkDecTree("../data/run456.root", "../data/vrt_r456_Lb_sigMC_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run458.root", "../data/vrt_r458_Lb_data_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDispl", "HLT_matched", 3.0);

    mkDecTree("../data/run456.root", "../data/vrt_r456_Lb_sigMC_Lb_lb12_barrel_lb.root", "lb", "rqha1>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run456.root", "../data/vrt_r456_Lb_sigMC_Lb_lb12_barrel_lbbar.root", "lb", "rqha1<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run456.root", "../data/vrt_r456_Lb_sigMC_Lb_lb12_barrel_cow.root", "lb", "isCowboy==1", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run456.root", "../data/vrt_r456_Lb_sigMC_Lb_lb12_barrel_sea.root", "lb", "isCowboy==0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    */

    /*
    mkDecTree("../data/run462.root", "../data/vrt_r462_Jp_bgrMC_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);
    mkDecTree("../data/run463.root", "../data/vrt_r463_Bs_bgrMC_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);
    mkDecTree("../data/run464.root", "../data/vrt_r464_Bp_bgrMC_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);
    mkDecTree("../data/run465.root", "../data/vrt_r465_B0_bgrMC_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);
    mkDecTree("../data/run466.root", "../data/vrt_r466_Lb_offMC_Lb_lb12_barrel.root", "lb", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);

    mkDecTree("../data/run462.root", "../data/vrt_r462_Jp_bgrMC_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    mkDecTree("../data/run463.root", "../data/vrt_r463_Bs_bgrMC_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    mkDecTree("../data/run464.root", "../data/vrt_r464_Bp_bgrMC_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    mkDecTree("../data/run465.root", "../data/vrt_r465_B0_bgrMC_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    mkDecTree("../data/run466.root", "../data/vrt_r466_Lb_offMC_Lb_lb12_displ.root", "lb", "", "HLT_jpsiDisplMCPseudo", "", 3.0);

    mkDecTree("../data/run467.root", "../data/vrt_r467_Jp_bgrMC_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);
    mkDecTree("../data/run468.root", "../data/vrt_r468_Bs_bgrMC_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);
    mkDecTree("../data/run469.root", "../data/vrt_r469_Bp_bgrMC_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);
    mkDecTree("../data/run470.root", "../data/vrt_r470_B0_offMC_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);
    mkDecTree("../data/run471.root", "../data/vrt_r471_Lb_bgrMC_B0_B006_barrel.root", "B0", "", "HLT_jpsiBarrelMCPseudo", "", 3.0);

    mkDecTree("../data/run467.root", "../data/vrt_r467_Jp_bgrMC_B0_B006_displ.root", "B0", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    mkDecTree("../data/run468.root", "../data/vrt_r468_Bs_bgrMC_B0_B006_displ.root", "B0", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    mkDecTree("../data/run469.root", "../data/vrt_r469_Bp_bgrMC_B0_B006_displ.root", "B0", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    mkDecTree("../data/run470.root", "../data/vrt_r470_B0_offMC_B0_B006_displ.root", "B0", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    mkDecTree("../data/run471.root", "../data/vrt_r471_Lb_bgrMC_B0_B006_displ.root", "B0", "", "HLT_jpsiDisplMCPseudo", "", 3.0);
    */

    // some crosscheck samples
    /*
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_phiplus.root", "B0", "phibc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_phiminus.root", "B0", "phibc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_phiLTpihalve.root", "B0", "TMath::Abs(phibc)>.5*TMath::Pi()", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_phiGTpihalve.root", "B0", "TMath::Abs(phibc)<=.5*TMath::Pi()", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_etaplus.root", "B0", "etabc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_etaminus.root", "B0", "etabc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_runA.root", "B0", "run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_runB.root", "B0", "run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_cow.root", "B0", "isCowboy==1", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_sea.root", "B0", "isCowboy==0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_seaA.root", "B0", "isCowboy==0&&run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_seaB.root", "B0", "isCowboy==0&&run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_cowA.root", "B0", "isCowboy==1&&run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_cowB.root", "B0", "isCowboy==1&&run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_PVlo.root", "B0", "nPV<=6", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_PVhi.root", "B0", "nPV>6", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_ptlo.root", "B0", "ptbc<17", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run473.root", "../data/vrt_r473_B0_data_B0_B006_barrel_pthi.root", "B0", "ptbc>=17", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    */

    /*
    double binPt[13] = { 0, 14.214, 15.395, 16.464, 17.531, 18.689, 19.916, 21.298, 23.076, 25.405, 28.578, 33.571, 126.328};
    for (int i=0; i!=12; i++) { mkDecTree(infileB0, outfilestemB0+"_ptbin"+toString(i)+".root", "B0", "ptbc>="+toString(binPt[i])+"&&ptbc<"+toString(binPt[i+1]), "HLT_jpsiBarrel", "HLT_matched", 3.0); };

    double binEta[13] = { -2.5, -1.5037649, -0.9562179, -0.6916759, -0.4587945, -0.229798, 0, 0.229798, 0.4587945, 0.6916759, 0.9562179, 1.5037649, 2.5};
    for (int i=0; i!=12; i++) { mkDecTree(infileB0, outfilestemB0+"_etabin"+toString(i)+".root", "B0", "etabc>="+toString(binEta[i])+"&&etabc<"+toString(binEta[i+1]), "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    */

    /*
    double binPhi[13] = { -3.142, -2.618, -2.094, -1.571, -1.047, -0.524, 0.000, 0.524, 1.047, 1.571, 2.094, 2.618, 3.142};
    for (int i=0; i!=12; i++) { mkDecTree(infileB0, outfilestemB0+"_bcphibin"+toString(i)+".root", "B0", "phibc>="+toString(binPhi[i])+"&&phibc<"+toString(binPhi[i+1]), "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    for (int i=0; i!=12; i++) { mkDecTree(infileB0, outfilestemB0+"_jpphibin"+toString(i)+".root", "B0", "phijp>="+toString(binPhi[i])+"&&phijp<"+toString(binPhi[i+1]), "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    for (int i=0; i!=12; i++) { mkDecTree(infileB0, outfilestemB0+"_rsphibin"+toString(i)+".root", "B0", "phirs>="+toString(binPhi[i])+"&&phirs<"+toString(binPhi[i+1]), "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    */

    /*
    const string infileLb = "../data/run480.root";
    const string outfilestemLb = "../data/vrt_r480_data_lb_acc_barrel";
    const int n_bin_alpha = 4, n_bin_drs = 4, n_bin_drssig = 4;
    double bin_alpha[n_bin_alpha]   = { .005, .012, .020, .1 };
    double bin_drs[n_bin_drs]       = { 0.0, 0.5, 1.0, 1.5};
    double bin_drssig[n_bin_drssig] = { 0.0, 1.0, 2.0, 5.0};
    for (int i=0; i!=n_bin_alpha; i++)
	for(int j=0; j!=n_bin_drs; j++)
	    for(int k=0; k!=n_bin_drssig; k++)
	    {
		const string addcut = "alphars<"+toString(bin_alpha[i])+"&&d3rs>"+toString(bin_drs[j])+"&&d3rs/d3Ers>"+toString(bin_drssig[k]);
		cout << "Additional cut: " << addcut << endl;
		mkDecTree(infileLb, outfilestemLb+"_cutbin_"+toString(i)+"_"+toString(j)+"_"+toString(k)+".root", "lb", addcut, "HLT_jpsiBarrel", "HLT_matched", 3.0);
	    }
    */

    /*
    const string infileLb = "../data/run480.root";
    const string outfilestemLb = "../data/vrt_r480_data_lb_acc_barrel_ser2";
    const int n_bin_ptha1 = 5;
    //double bin_ptha1[n_bin_ptha1]   = { 2.0, 1.9, 1.8, 1.7, 1.6 };
    double bin_ptha1[n_bin_ptha1]   = { 2.4, 2.8, 3.2, 3.6, 4.0 };
    //double bin_ptha1[n_bin_ptha1]   = { 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0 };
    for (int i=0; i!=n_bin_ptha1; i++)
    {
	const string addcut = "rptha1>"+toString(bin_ptha1[i]);
	cout << "Changed cut: " << addcut << endl;
	mkDecTree(infileLb, outfilestemLb+"_ptpr_lball_"+toString(i)+".root", "lb", addcut, "HLT_jpsiBarrel", "HLT_matched", 3.0);
	mkDecTree(infileLb, outfilestemLb+"_ptpr_lb_"+toString(i)+".root", "lb", addcut+"&&rqha1>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
	mkDecTree(infileLb, outfilestemLb+"_ptpr_lbbar_"+toString(i)+".root", "lb", addcut+"&&rqha1<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    }
    */

    /*
    const string infileB0 = "../data/run460_472.root";
    const string outfilestemB0 = "../data/vrt_r460_472_MC_B0_acc_barrel";
    int i=0;
    for (double pt=0.5; pt<=1.8; pt+=.1)
    {
	const string addcut = "(rqha1>0?rptha1:rptha2)>"+toString(pt)+"&&(rqha1<0?rptha1:rptha2)>"+toString(2.3-pt);
	cout << "Changed cut: " << addcut << endl;
	mkDecTree(infileB0, outfilestemB0+"_ptbins_B0_"+toString(i)+".root", "B0", addcut, "HLT_jpsiBarrel", "HLT_matched", 3.0);
	i++;
    }
    */

    /*
    const string fileLb = "../data/run477_478.root";
    const string outfilestemLb = "../data/vrt_r477_478_data_lb_acc_barrel";
    mkDecTree(fileLb, outfilestemLb+".root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_phiplus.root", "lb", "phibc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_phiminus.root", "lb", "phibc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_phiLTpihalve.root", "lb", "TMath::Abs(phibc)>.5*TMath::Pi()", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_phiGTpihalve.root", "lb", "TMath::Abs(phibc)<=.5*TMath::Pi()", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_etaplus.root", "lb", "etabc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_etaminus.root", "lb", "etabc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_runA.root", "lb", "run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_runB.root", "lb", "run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_lb.root", "lb", "rqha1>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_lbbar.root", "lb", "rqha1<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_cow.root", "lb", "isCowboy==1", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_sea.root", "lb", "isCowboy==0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_seaA.root", "lb", "isCowboy==0&&run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_seaB.root", "lb", "isCowboy==0&&run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_cowA.root", "lb", "isCowboy==1&&run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_cowB.root", "lb", "isCowboy==1&&run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_PVlo.root", "lb", "nPV<=6", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_PVhi.root", "lb", "nPV>6", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_ptlo.root", "lb", "ptbc<17", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(fileLb, outfilestemLb+"_pthi.root", "lb", "ptbc>=17", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    */

    /*
    // for misalignmen studies
    // B0
    //for (int i=485; i!=503; i++) { mkDecTree("../data/run"+toString(i)+".root", "../data/vrt_r"+toString(i)+"_B0_MC_misali_B0_B006_barrel_MCmatch.root", "B0", "isMCmatch==1", "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    for (int i=527; i!=532; i++) { mkDecTree("../data/run"+toString(i)+".root", "../data/vrt_r"+toString(i)+"_B0_MC_misali_B0_B006_barrel_MCmatch.root", "B0", "isMCmatch==1", "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    // Lb
    //for (int i=503; i<=520; i++) { mkDecTree("../data/run"+toString(i)+".root", "../data/vrt_r"+toString(i)+"_Lb_MC_misali_Lb_lb12_barrel_MCmatch.root", "lb", "isMCmatch==1", "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    //for (int i=503; i<=520; i++) { mkDecTree("../data/run"+toString(i)+".root", "../data/vrt_r"+toString(i)+"_Lb_MC_misali_Lb_lb12_barrel_MCmatch_Lb.root", "lb", "isMCmatch==1&&rqha1>0", "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    //for (int i=503; i<=520; i++) { mkDecTree("../data/run"+toString(i)+".root", "../data/vrt_r"+toString(i)+"_Lb_MC_misali_Lb_lb12_barrel_MCmatch_Lbbar.root", "lb", "isMCmatch==1&&rqha1<0", "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    const int istart=525, iend=526;
    for (int i=istart; i<=iend; i++) { mkDecTree("../data/run"+toString(i)+".root", "../data/vrt_r"+toString(i)+"_Lb_MC_misali_Lb_lb12_barrel_MCmatch.root", "lb", "isMCmatch==1", "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    for (int i=istart; i<=iend; i++) { mkDecTree("../data/run"+toString(i)+".root", "../data/vrt_r"+toString(i)+"_Lb_MC_misali_Lb_lb12_barrel_MCmatch_Lb.root", "lb", "isMCmatch==1&&rqha1>0", "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    for (int i=istart; i<=iend; i++) { mkDecTree("../data/run"+toString(i)+".root", "../data/vrt_r"+toString(i)+"_Lb_MC_misali_Lb_lb12_barrel_MCmatch_Lbbar.root", "lb", "isMCmatch==1&&rqha1<0", "HLT_jpsiBarrel", "HLT_matched", 3.0); };
    */

    mkDecTree("../data/run552_556.root", "../data/vrt_r552_556_mc_lb_cuts.root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run552_556.root", "../data/vrt_r552_556_mc_lb_sig.root", "lb", "isSig==1", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree("../data/run552_556.root", "../data/vrt_r552_556_mc_lb_match.root", "lb", "isMCmatch==1", "HLT_jpsiBarrel", "HLT_matched", 3.0);
}

// example: infileLb = "../data/run480.root"; outfilestemLb = "../data/vrt_r480_lb_data_lb_lb13exp_barrel";
void mkSomeDecTreesXchecksLb(string infileLb, string outfilestemLb)
{
    mkDecTree(infileLb, outfilestemLb+".root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_phiplus.root", "lb", "phibc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_phiminus.root", "lb", "phibc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_phiLTpihalve.root", "lb", "TMath::Abs(phibc)>.5*TMath::Pi()", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_phiGTpihalve.root", "lb", "TMath::Abs(phibc)<=.5*TMath::Pi()", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_etaplus.root", "lb", "etabc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_etaminus.root", "lb", "etabc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_runA.root", "lb", "run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_runB.root", "lb", "run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_lb.root", "lb", "rqha1>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_lbbar.root", "lb", "rqha1<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_cow.root", "lb", "isCowboy==1", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_sea.root", "lb", "isCowboy==0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_seaA.root", "lb", "isCowboy==0&&run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_seaB.root", "lb", "isCowboy==0&&run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_cowA.root", "lb", "isCowboy==1&&run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_cowB.root", "lb", "isCowboy==1&&run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_PVlo.root", "lb", "nPV<=6", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_PVhi.root", "lb", "nPV>6", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_ptlo.root", "lb", "ptbc<17", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_pthi.root", "lb", "ptbc>=17", "HLT_jpsiBarrel", "HLT_matched", 3.0);
}

void mkDecTreesMisaliLb(string infileLb, string outfilestemLb)
{
    mkDecTree(infileLb, outfilestemLb+".root", "lb", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_lb.root", "lb", "rqha1>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileLb, outfilestemLb+"_lbbar.root", "lb", "rqha1<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
}

void mkSomeDecTreesMisaliLb()
{
    typedef vector< pair<string,string> > list_t;
    const string path = "../data/";
    list_t liste;
    liste.push_back(make_pair("592","bowingEpsilon_minus"));
    liste.push_back(make_pair("570","bowingEpsilon_plus"));
    liste.push_back(make_pair("571","ellipticalEpsilon_minus"));
    //liste.push_back(make_pair("572","ellipticalEpsilon_plus"));
    //liste.push_back(make_pair("573","layerRotEpsilon_minus"));
    //liste.push_back(make_pair("574","layerRotEpsilon_plus"));
    //liste.push_back(make_pair("575","nomisali"));
    liste.push_back(make_pair("576","radialEpsilon_layer1minus2"));
    //liste.push_back(make_pair("577","radialEpsilon_layer1plus2"));
    //liste.push_back(make_pair("578","radialEpsilon_minus"));
    //liste.push_back(make_pair("579","radialEpsilon_plus"));
    liste.push_back(make_pair("590","saggitaEpsilon_minus"));
    liste.push_back(make_pair("591","saggitaEpsilon_plus"));
    liste.push_back(make_pair("582","skewEpsilon_minus2"));
    liste.push_back(make_pair("583","skewEpsilon_plus2"));
    //liste.push_back(make_pair("584","telescopeEpsilon_minus"));
    //liste.push_back(make_pair("585","telescopeEpsilon_plus"));
    //liste.push_back(make_pair("586","twistEpsilon_minus"));
    //liste.push_back(make_pair("587","twistEpsilon_plus"));
    //liste.push_back(make_pair("588","zExpEpsilon_minus"));
    //liste.push_back(make_pair("589","zExpEpsilon_plus"));

    for (list_t::const_iterator it=liste.begin(); it!=liste.end(); it++)
    {
	mkDecTreesMisaliLb(path+"run"+it->first+".root", path+"vrt_r"+it->first+"_"+it->second);
    }
    cout << "List of files created:" << endl;
    for (list_t::const_iterator it=liste.begin(); it!=liste.end(); it++)
    {
	cout << path+"vrt_r"+it->first+"_"+it->second << endl;
    }
}

void mkSomeDecTreesCutvarLb(string infileLb, string outfilestemLb, string cutToVary, string cutLabel, int cutN, double cutLo, double cutHi)
{
    Cuts cut;
    cut.selectCut("acc06Lb", "muSoft", "lb13exp", "HLT_jpsiBarrel", "HLT_matched");
    const double cutStep = (cutHi-cutLo) / (double)cutN;
    cutN++;
    for (int i = 0; i!=cutN; i++)
    {
	const double cutVal = cutLo + i*cutStep;
	string cutstring = cut.getCutChangeOne(cutToVary, cutVal);
	cout << "-----------------------------------" << endl;
	cout << i << ": " << cutstring << endl;
	const string targetFileName = outfilestemLb+cutLabel+"_"+toString(i)+".root";
	mkDecTree(infileLb, cutstring, targetFileName, "", 3.0);
    }
}

// example: infileB0 = "../data/run479.root"; outfilestemB0 = "../data/vrt_r479_B0_data_B0_B007_barrel"
void mkSomeDecTreesXchecksB0(string infileB0, string outfilestemB0)
{
    mkDecTree(infileB0, outfilestemB0+".root", "B0", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_phiplus.root", "B0", "phibc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_phiminus.root", "B0", "phibc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_phiLTpihalve.root", "B0", "TMath::Abs(phibc)>.5*TMath::Pi()", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_phiGTpihalve.root", "B0", "TMath::Abs(phibc)<=.5*TMath::Pi()", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_etaplus.root", "B0", "etabc>0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_etaminus.root", "B0", "etabc<0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_runA.root", "B0", "run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_runB.root", "B0", "run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_cow.root", "B0", "isCowboy==1", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_sea.root", "B0", "isCowboy==0", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_seaA.root", "B0", "isCowboy==0&&run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_seaB.root", "B0", "isCowboy==0&&run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_cowA.root", "B0", "isCowboy==1&&run<174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_cowB.root", "B0", "isCowboy==1&&run>174000", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_PVlo.root", "B0", "nPV<=6", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_PVhi.root", "B0", "nPV>6", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_ptlo.root", "B0", "ptbc<17", "HLT_jpsiBarrel", "HLT_matched", 3.0);
    mkDecTree(infileB0, outfilestemB0+"_pthi.root", "B0", "ptbc>=17", "HLT_jpsiBarrel", "HLT_matched", 3.0);
}

void mkSomeDecTreesMisaliB0(string infileB0, string outfilestemB0)
{
    mkDecTree(infileB0, outfilestemB0+".root", "B0", "", "HLT_jpsiBarrel", "HLT_matched", 3.0);
}


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

    const double tau_truth(isB0 ? 1.536e-12 : 1.507e-12);
    const double taureso_truth(0.2e-12);
    const double tauresoE_truth(.04e-12);

    const double reducedSigma(3);

    // background
    const double mass_lo(isB0 ? 5.16 : 5.4);
    const double mass_hi(isB0 ? 5.75 : 6.0);
    const double mass_window(mass_hi-mass_lo);

    const double bgr_prompt_width(0.2e-12);
    const double bgr_nonprompt_tau1(1.15e-12);
    const double bgr_nonprompt_tau2(0.75e-12);
    const double bgr_nonprompt_frac1(.99);
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

