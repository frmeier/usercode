#include "TF1.h"
#include "TH1F.h"
#include <iostream>

void fitgaus(TH1F &h){
    TF1 *f1 = new TF1("f1", "[0] + [1]*x + gaus(2)", 5.0, 6.);
    f1->SetParameters(-10,0,1,5.6,.1);
    f1->SetParNames("const","slope","scale","mean","sigma");
    h.Fit(f1);
    const double pConst = f1->GetParameter(0);
    const double pSlope = f1->GetParameter(1);
    const double pScale = f1->GetParameter(2);
    const double pMean  = f1->GetParameter(3);
    const double pSigma = f1->GetParameter(4);
    const double nSigmas = 2.0;
    const double lowerBound = pMean-nSigmas*pSigma;
    const double upperBound = pMean+nSigmas*pSigma;
    std::cout << "Mean:  " << pMean << std::endl;
    std::cout << "Sigma: " << pSigma << std::endl;
    std::cout << "lower: " << lowerBound << std::endl;
    std::cout << "upper: " << upperBound << std::endl;
    std::cout << "Integral +/- " << nSigmas << " sigmas: "
	<< f1->Integral(lowerBound,upperBound)*h.GetBinWidth(1) << std::endl;;
}

void fitgaus(TH1F* h, double mass, double min, double max){
    TF1 *f1 = new TF1("f1", "[0] + [1]*x + gaus(2)", min, max);
    f1->SetParameters(1,0,1,mass,.1);
    f1->SetParNames("const","slope","scale","mean","sigma");
    h->Fit(f1);
    const double pConst = f1->GetParameter(0);
    const double pSlope = f1->GetParameter(1);
    const double pScale = f1->GetParameter(2);
    const double pMean  = f1->GetParameter(3);
    const double pSigma = f1->GetParameter(4);
    const double nSigmas = 2.5;
    const double lowerBound = pMean-nSigmas*pSigma;
    const double upperBound = pMean+nSigmas*pSigma;
    const double binWidth = h->GetBinWidth(1);
    const double sigPlNoise= f1->Integral(lowerBound,upperBound)/binWidth;
    std::cout << "Mean:  " << pMean << std::endl;
    std::cout << "Sigma: " << pSigma << std::endl;
    std::cout << "lower: " << lowerBound << std::endl;
    std::cout << "upper: " << upperBound << std::endl;
    std::cout << "bin width: " << binWidth << std::endl;
    std::cout << "Integral +/- " << nSigmas << " sigmas: "
	<< sigPlNoise << std::endl;;
    TF1 *f2 = new TF1("f2", "gaus", min, max);
    f2->SetParameters(pScale,pMean,pSigma);
    f2->Draw("same");
    const double sig = f2->Integral(lowerBound,upperBound)/binWidth;
    std::cout << "Integral sig only: " << sig << std::endl;
    std::cout << "S/sqrt(S+N): " << sig/sqrt(sigPlNoise) << std::endl;
    std::cout << "S/sqrt(N): " << sig/sqrt(sigPlNoise-sig) << std::endl;
}

void fitgaus(TH1F* h){
    fitgaus(h,5.6,5,6);
}


