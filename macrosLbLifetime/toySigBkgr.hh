#ifndef GUARD_TOYSIGBKGR
#define GUARD_TOYSIGBKGR

#include "TH1F.h"
#include "TRandom3.h"

class toySigBkgr
{
    public:
	toySigBkgr();

	void drawSig(unsigned int N);
	void drawBkgr(unsigned int N);
	void draw();
	void drawFit();
	void doFit();
	void resetFit();
	void printSN();

    private:
	TH1F* histo;
	TF1* fitfnc, *fitsig, *fitbgr;
	TRandom3* rng;
	double sigMean, sigSigma; // distribution parameters for signal trruth
	unsigned int nSig, nBkgr; // number of diced S and B values
	double calcSig, calcBkgr; // calculated S and B from fit
	double pConst, pSlope, pScale, pMean, pSigma; // parameters obtained from fit
	double fitSigmas; // how many sigmas the fit should cover
	int nBins; // numberof bins for the histo
	double min,max; // range
	double binWidth;
};

#endif
