#include "TROOT.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH2I.h"
#include <map>
#include "utils.h"
#include "TriggerPlot.h"

TriggerPlot::TriggerPlot(TTree *t) : maxRunsPerHisto_(40)
{
    t_ = t;
    initBranchAddresses();
    initRunmap();
}

void TriggerPlot::initRunmap(){
    runmap_.clear();
    // create an ordered list of runs in the tree
    for (int i = 0; i!=t_->GetEntries(); i++)
    {
	t_->GetEntry(i);
	runmap_[fRun]=1;
    }
    {
	int i(0);
	for (runmap_t::const_iterator it=runmap_.begin(); it!=runmap_.end(); it++)
	    runmap_[it->first] = ++i;
    }
}

void TriggerPlot::initBranchAddresses()
{
    t_->SetBranchAddress("run",&fRun);

    t_->SetBranchAddress("HLTmatch",&fHLTmatch);

    t_->SetBranchAddress("HLTqrk",&fHLTqrk);
    t_->SetBranchAddress("HLTDMu3jp",&fHLTDMu3jp);
    t_->SetBranchAddress("HLTDMu6p5BarJp", &fHLTDMu6p5BarJp); // Dimuon6p5_Barrel_Jpsi_v1
    t_->SetBranchAddress("HLTDMu6p5JpDis", &fHLTDMu6p5JpDis); // Dimuon6p5_Jpsi_Displaced_v1
    t_->SetBranchAddress("HLTDMu6p5Jp", &fHLTDMu6p5Jp); // Dimuon6p5_Jpsi_v1
    t_->SetBranchAddress("HLTMu5L2Mu2Jpsi", &fHLTMu5L2Mu2Jpsi); // Mu5_L2Mu2_Jpsi_v3
    t_->SetBranchAddress("HLTMu5Tr2Jpsi", &fHLTMu5Tr2Jpsi); // Mu5_Track2_Jpsi_v2
    t_->SetBranchAddress("HLTMu5Tr7Jpsi", &fHLTMu5Tr7Jpsi); // Mu7_Track7_Jpsi_v3
    t_->SetBranchAddress("HLTDMu10BarJp", &fHLTDMu10BarJp); // Dimuon10_Barrel_Jpsi_v1
    t_->SetBranchAddress("HLTDMu7JpDis",&fHLTDMu7JpDis);

    t_->SetBranchAddress("HLTokJpsi",&fHLTokJpsi);
    t_->SetBranchAddress("HLTokBarrelJpsi",&fHLTokBarrelJpsi);
    t_->SetBranchAddress("HLTokDisplJpsi",&fHLTokDisplJpsi);

    t_->SetBranchAddress("HLTSingleMu",&fHLTSingleMu); // Single Mu family for efficiencies
    t_->SetBranchAddress("HLTSingleIsoMu",&fHLTSingleIsoMu); // Single Mu family for efficiencies
    t_->SetBranchAddress("HLTSingleL1Mu",&fHLTSingleL1Mu); // Single Mu family for efficiencies
    t_->SetBranchAddress("HLTSingleL2Mu",&fHLTSingleL2Mu); // Single Mu family for efficiencies
    t_->SetBranchAddress("HLTSingleHLTMu",&fHLTSingleHLTMu); // Single Mu family for efficiencies
}

void TriggerPlot::initHistoLabels(TH2I *h)
{
    h->GetYaxis()->SetBinLabel(1,"all events");

    h->GetYaxis()->SetBinLabel(2,"HLT muon match");

    h->GetYaxis()->SetBinLabel(3,"HLTDMuX_Qrk_vX");
    h->GetYaxis()->SetBinLabel(4,"HLTDMu3_Jpsi_vX");
    h->GetYaxis()->SetBinLabel(5,"HLTDMu6p5BarJp");
    h->GetYaxis()->SetBinLabel(6,"HLTDMu6p5JpDis");
    h->GetYaxis()->SetBinLabel(7,"HLTDMu6p5Jp");
    h->GetYaxis()->SetBinLabel(8,"HLTMu5L2Mu2Jpsi");
    h->GetYaxis()->SetBinLabel(9,"HLTMu5Tr2Jpsi");
    h->GetYaxis()->SetBinLabel(10,"HLTMu5Tr7Jpsi");
    h->GetYaxis()->SetBinLabel(11,"HLTDMu10BarJp");
    h->GetYaxis()->SetBinLabel(12,"HLTDMu7JpDis");

    h->GetYaxis()->SetBinLabel(13,"HLTSum_Jpsi");
    h->GetYaxis()->SetBinLabel(14,"HLTSum_JpsiBarrel");
    h->GetYaxis()->SetBinLabel(15,"HLTSum_JpsiDispl");
    //h->GetYaxis()->SetBinLabel(,"");

    h->GetYaxis()->SetBinLabel(16,"HLTSingleMu");
    h->GetYaxis()->SetBinLabel(17,"HLTSingleIsoMu");
    h->GetYaxis()->SetBinLabel(18,"HLTSingleL1Mu");
    h->GetYaxis()->SetBinLabel(19,"HLTSingleL2Mu");
    h->GetYaxis()->SetBinLabel(20,"HLTSingleHLTMu");
}

void TriggerPlot::drawHisto(int drawNo)
{
    if (drawNo>getNhistos()) return;
    if (drawNo<0) return;
    // now lets make the histo
    const int nTriggerBins(20);
    // determine boundaries of histo
    unsigned int nRunBins, runIdxMin, runIdxMax;
    if (getNhistos() == 1)
    {
	nRunBins = runmap_.size();
	runIdxMin = 1;
	runIdxMax = runmap_.size();
    }
    else
    {
	nRunBins = maxRunsPerHisto_;
	runIdxMin = drawNo*maxRunsPerHisto_+1;
	runIdxMax = (drawNo+1)*maxRunsPerHisto_;
	if (runIdxMax > runmap_.size()) runIdxMax = runmap_.size();
    }
    // doing the histo
    TH2I *h = new TH2I(("htrig"+toString(drawNo)).c_str(),("Trigger paths "+toString(drawNo+1)+"/"+toString(getNhistos())).c_str(),
	    nRunBins,runIdxMin,runIdxMin+nRunBins, nTriggerBins,0,nTriggerBins);
    for (int i = 0; i!=t_->GetEntries(); i++)
    {
	t_->GetEntry(i);
	const unsigned int curIdx = runmap_[fRun];
	if (curIdx < runIdxMin) continue;
	if (curIdx > runIdxMax) continue;
	// indices here are bin number -1
	h->Fill(curIdx,0); // gives total number of entries in run

	if (fHLTmatch)        h->Fill(curIdx,1); // trigger matcher result

	if (fHLTqrk)          h->Fill(curIdx,2); // individual triggers
	if (fHLTDMu3jp)       h->Fill(curIdx,3);
	if (fHLTDMu6p5BarJp)  h->Fill(curIdx,4);
	if (fHLTDMu6p5JpDis)  h->Fill(curIdx,5);
	if (fHLTDMu6p5Jp)     h->Fill(curIdx,6);
	if (fHLTMu5L2Mu2Jpsi) h->Fill(curIdx,7);
	if (fHLTMu5Tr2Jpsi)   h->Fill(curIdx,8);
	if (fHLTMu5Tr7Jpsi)   h->Fill(curIdx,9);
	if (fHLTDMu10BarJp)   h->Fill(curIdx,10);
	if (fHLTDMu7JpDis)    h->Fill(curIdx,11);

	if (fHLTokJpsi)       h->Fill(curIdx,12); // summary results
	if (fHLTokBarrelJpsi) h->Fill(curIdx,13);
	if (fHLTokDisplJpsi)  h->Fill(curIdx,14);

	if (fHLTSingleMu)     h->Fill(curIdx,15);
	if (fHLTSingleIsoMu)  h->Fill(curIdx,16);
	if (fHLTSingleL1Mu)   h->Fill(curIdx,17);
	if (fHLTSingleL2Mu)   h->Fill(curIdx,18);
	if (fHLTSingleHLTMu)  h->Fill(curIdx,19);
    }
    // set bin labels
    {
	int i(0);
	for (runmap_t::const_iterator it=runmap_.begin(); it!=runmap_.end(); it++)
	{
	    const unsigned int curIdx = it->second;
	    if (curIdx < runIdxMin) continue;
	    if (curIdx > runIdxMax) continue;
	    i++;
	    h->GetXaxis()->SetBinLabel(i,toString(it->first).c_str());
	}
    }
    initHistoLabels(h);
    h->SetStats(false);
    h->GetXaxis()->SetLabelSize(0.03);
    h->GetXaxis()->LabelsOption("v");
    h->GetXaxis()->SetTitle("Run number");
    h->GetXaxis()->SetTitleSize(0.03);
    h->GetXaxis()->SetTitleOffset(2.4);
    h->GetYaxis()->SetTitleOffset(3.0);
    h->GetYaxis()->SetTitle("HLT path");
    h->Draw("COLZtext90");
    gStyle->SetPalette(1);
    gPad->SetLeftMargin(0.22);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.10);
}

