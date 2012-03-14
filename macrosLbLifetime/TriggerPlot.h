#ifndef TRIGGERPLOT_H_GUARD
#define TRIGGERPLOT_H_GUARD

#include "TTree.h"
#include <map>
#include "TMath.h"
#include "TH2I.h"

class TriggerPlot
{
    public:
	typedef std::map<unsigned int,int> runmap_t;
	TriggerPlot(TTree *t);
	int getNhistos() { return TMath::CeilNint( (double)runmap_.size()/(double)maxRunsPerHisto_); };
	void drawHisto(int i);

    private:
	runmap_t runmap_;
	TTree *t_;
	int maxRunsPerHisto_;

	void initRunmap();
	void initBranchAddresses();
	void initHistoLabels(TH2I *h);

	int fRun;
	bool fHLTmatch;
	bool fHLTokJpsi, fHLTokBarrelJpsi, fHLTokDisplJpsi;
	bool fHLTDMu6p5BarJp, fHLTDMu6p5JpDis, fHLTDMu6p5Jp, fHLTMu5L2Mu2Jpsi, fHLTMu5Tr2Jpsi, fHLTMu5Tr7Jpsi, fHLTDMu10BarJp, fHLTDMu7JpDis;
	bool fHLTqrk, fHLTDMu3jp;
	bool fHLTSingleMu, fHLTSingleIsoMu, fHLTSingleL1Mu, fHLTSingleL2Mu, fHLTSingleHLTMu;
};

#endif
