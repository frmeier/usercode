#include "toySigBkgr.hh"
#include "TF1.h"
#include "TH1F.h"
#include <iostream>

using std::cout;
using std::endl;

toySigBkgr::toySigBkgr() : sigMean(.5), sigSigma(.05), nSig(0), nBkgr(0), fitSigmas(2.0),
    nBins(20),min(0.0),max(1.0)
{
    histo = new TH1F("histo","toySigBkgr",nBins,min,max);
    rng = new TRandom3;
    fitfnc = new TF1("fitfnc", "[0] + [1]*x + gaus(2)", min, max);
    resetFit();
    fitfnc->SetParNames("const","slope","scale","mean","sigma");
    fitsig = new TF1("fitsig", "gaus", min,max);
    fitbgr = new TF1("fitbgr", "[0] + [1]*x", min,max);
    binWidth = (max-min)/nBins;
}

void toySigBkgr::drawSig(unsigned int N)
{
    for (unsigned int i = 0; i!=N; i++)
    {
	histo->Fill(rng->Gaus(sigMean,sigSigma));
    }
    nSig+=N;
}

void toySigBkgr::drawBkgr(unsigned int N)
{
    for (unsigned int i = 0; i!=N; i++)
    {
	histo->Fill(rng->Uniform());
    }
    nBkgr+=N;
}

void toySigBkgr::draw()
{
    histo->Draw();
}

void toySigBkgr::drawFit()
{
}

void toySigBkgr::doFit()
{
    histo->Fit(fitfnc);
    pConst = fitfnc->GetParameter(0);
    pSlope = fitfnc->GetParameter(1);
    pScale = fitfnc->GetParameter(2);
    pMean  = fitfnc->GetParameter(3);
    pSigma = fitfnc->GetParameter(4);
    fitsig->SetParameters(pScale,pMean,pSigma);
    fitbgr->SetParameters(pConst,pSlope);
    const double intLo = pMean-fitSigmas*pSigma;
    const double intUp = pMean+fitSigmas*pSigma;
    calcSig = fitsig->Integral(intLo,intUp)/binWidth;
    calcBkgr = fitbgr->Integral(intLo,intUp)/binWidth;
}

void toySigBkgr::resetFit()
{
    fitfnc->SetParameters(0,0,1,.5,.1);
}

void toySigBkgr::printSN()
{
    cout << "Diced:  S = " << nSig << " B = " << nBkgr << endl;
    const double truthBgr = nBkgr*(2*sigSigma*fitSigmas)*(max-min);
    cout << "Truth:  S = " << nSig << " B = " << truthBgr << endl;
    const double truthSB = nSig/sqrt(truthBgr);
    const double truthSSB = nSig/sqrt(nSig+truthBgr);
    cout << "     S/sqrt(B) = " << truthSB << " S/sqrt(S+B) = " << truthSSB << endl;
    cout << "Observed: S = " << calcSig << " B = " << calcBkgr << endl;
    if (calcSig > 0 && calcBkgr > 0)
    {
	const double calcSB = calcSig/sqrt(calcBkgr);
	const double calcSSB = calcSig/sqrt(calcSig+calcBkgr);
	cout << "     S/sqrt(B) = " << calcSB << " S/sqrt(S+B) = " << calcSSB << endl;
    }
}
