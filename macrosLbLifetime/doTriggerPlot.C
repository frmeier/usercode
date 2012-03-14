#include "TTree.h"
#include "TH2I.h"
#include <map>
#include <iostream>
#include "utils.h"

void doTriggerPlot(TTree *t)
{
    // First create an ordered list of runs in the tree
    int run;
    typedef std::map<int,int> runmap_t;
    runmap_t runmap;
    t->SetBranchAddress("run",&run);
    for (int i = 0; i!=t->GetEntries(); i++)
    {
	t->GetEntry(i);
	runmap[run]=1;
    }
    {
	int i(0);
	for (runmap_t::const_iterator it=runmap.begin(); it!=runmap.end(); it++)
	    runmap[it->first] = ++i;
    }
    // now lets make the histo
    const int nTriggerBins(15);
    const int nRunBins(runmap.size());
    TH2I *h = new TH2I("htrig","Trigger paths",nRunBins,1,nRunBins+1, nTriggerBins,0,nTriggerBins);
    bool fHLTmatch, fHLTok;
    bool fHLTDMu6p5BarJp, fHLTDMu6p5JpDis, fHLTDMu6p5Jp, fHLTMu5L2Mu2Jpsi, fHLTMu5Tr2Jpsi, fHLTMu5Tr7Jpsi, fHLTDMu10BarJp, fHLTDMu7JpDis;
    bool fHLTSingleMu, fHLTSingleIsoMu, fHLTSingleL1Mu, fHLTSingleL2Mu, fHLTSingleHLTMu;
    t->SetBranchAddress("HLTmatch",&fHLTmatch);
    t->SetBranchAddress("HLTok",&fHLTok);
    t->SetBranchAddress("HLTDMu6p5BarJp", &fHLTDMu6p5BarJp); // Dimuon6p5_Barrel_Jpsi_v1
    t->SetBranchAddress("HLTDMu6p5JpDis", &fHLTDMu6p5JpDis); // Dimuon6p5_Jpsi_Displaced_v1
    t->SetBranchAddress("HLTDMu6p5Jp", &fHLTDMu6p5Jp); // Dimuon6p5_Jpsi_v1
    t->SetBranchAddress("HLTMu5L2Mu2Jpsi", &fHLTMu5L2Mu2Jpsi); // Mu5_L2Mu2_Jpsi_v3
    t->SetBranchAddress("HLTMu5Tr2Jpsi", &fHLTMu5Tr2Jpsi); // Mu5_Track2_Jpsi_v2
    t->SetBranchAddress("HLTMu5Tr7Jpsi", &fHLTMu5Tr7Jpsi); // Mu7_Track7_Jpsi_v3
    t->SetBranchAddress("HLTDMu10BarJp", &fHLTDMu10BarJp); // Dimuon10_Barrel_Jpsi_v1
    t->SetBranchAddress("HLTDMu7JpDis",&fHLTDMu7JpDis);
    t->SetBranchAddress("HLTSingleMu",&fHLTSingleMu); // Single Mu family for efficiencies
    t->SetBranchAddress("HLTSingleIsoMu",&fHLTSingleIsoMu); // Single Mu family for efficiencies
    t->SetBranchAddress("HLTSingleL1Mu",&fHLTSingleL1Mu); // Single Mu family for efficiencies
    t->SetBranchAddress("HLTSingleL2Mu",&fHLTSingleL2Mu); // Single Mu family for efficiencies
    t->SetBranchAddress("HLTSingleHLTMu",&fHLTSingleHLTMu); // Single Mu family for efficiencies
    for (int i = 0; i!=t->GetEntries(); i++)
    {
	t->GetEntry(i);
	if (fHLTmatch) h->Fill(runmap[run],0);
	if (fHLTok) h->Fill(runmap[run],1);
	if (fHLTDMu6p5BarJp) h->Fill(runmap[run],2);
	if (fHLTDMu6p5JpDis) h->Fill(runmap[run],3);
	if (fHLTDMu6p5Jp) h->Fill(runmap[run],4);
	if (fHLTMu5L2Mu2Jpsi) h->Fill(runmap[run],5);
	if (fHLTMu5Tr2Jpsi) h->Fill(runmap[run],6);
	if (fHLTMu5Tr7Jpsi) h->Fill(runmap[run],7);
	if (fHLTDMu10BarJp) h->Fill(runmap[run],8);
	if (fHLTDMu7JpDis) h->Fill(runmap[run],9);
	if (fHLTSingleMu) h->Fill(runmap[run],10);
	if (fHLTSingleIsoMu) h->Fill(runmap[run],11);
	if (fHLTSingleL1Mu) h->Fill(runmap[run],12);
	if (fHLTSingleL2Mu) h->Fill(runmap[run],13);
	if (fHLTSingleHLTMu) h->Fill(runmap[run],14);
    }
    // set bin labels
    {
	int i(0);
	for (runmap_t::const_iterator it=runmap.begin(); it!=runmap.end(); it++)
	{
	    i++;
	    h->GetXaxis()->SetBinLabel(i,toString(it->first).c_str());
	}
    }
    h->GetYaxis()->SetBinLabel(1,"HLTmatch");
    h->GetYaxis()->SetBinLabel(2,"HLTok");
    h->GetYaxis()->SetBinLabel(3,"HLTDMu6p5BarJp");
    h->GetYaxis()->SetBinLabel(4,"HLTDMu6p5JpDis");
    h->GetYaxis()->SetBinLabel(5,"HLTDMu6p5Jp");
    h->GetYaxis()->SetBinLabel(6,"HLTMu5L2Mu2Jpsi");
    h->GetYaxis()->SetBinLabel(7,"HLTMu5Tr2Jpsi");
    h->GetYaxis()->SetBinLabel(8,"HLTMu5Tr7Jpsi");
    h->GetYaxis()->SetBinLabel(9,"HLTDMu10BarJp");
    h->GetYaxis()->SetBinLabel(10,"HLTDMu7JpDis");
    h->GetYaxis()->SetBinLabel(11,"HLTSingleMu");
    h->GetYaxis()->SetBinLabel(12,"HLTSingleIsoMu");
    h->GetYaxis()->SetBinLabel(13,"HLTSingleL1Mu");
    h->GetYaxis()->SetBinLabel(14,"HLTSingleL2Mu");
    h->GetYaxis()->SetBinLabel(15,"HLTSingleHLTMu");
    h->Draw("COLZ");
}

