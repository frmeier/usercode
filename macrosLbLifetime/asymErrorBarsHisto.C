// Taken from https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars

// CMS Statistics Committee Recommendations
//
// Asymmetric Error Bars for Poisson Event Counts
//
// While a Poisson distribution of mean μ has a variance equal to μ, an interval [μ-sqrt(μ),μ+sqrt(μ)] may result in undercoverage, especially if μ is not large.
// 
// The Statistics Committee invites you to use asymmetric error bars with correct coverage for event counts with Poisson variates, when the number of entries in a significant portion of the spectrum are small. This may either arise in low-statistics histograms, or in high-statistics ones drawn on a semi-logarithmic scale.
//
// These "correct coverage" error bars, first derived by Garwood, are obtained from the Neyman construction using the central interval convention for N>0, and drawing the upper limit for N=0. The following piece of code exemplifies the calculation and the plotting of these error bars in a histogram, using ROOT (thanks Lorenzo Moneta): 

// To be able to compile the above code, insert the following preprocessor directives:
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
// - Thanks Guillelmo Gomez-Ceballos for these.

void asymErrorBarsHisto()
{
   const double alpha = 1 - 0.6827;
   TH1D * h1 = new TH1D("h1","h1",50,-4,4);
   h1->FillRandom("gaus",100);

   TGraphAsymmErrors * g = new TGraphAsymmErrors(h1);

   for (int i = 0; i < g->GetN(); ++i) {
      int N = g->GetY()[i];
      double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) : ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
      g->SetPointEYlow(i, N-L);
      g->SetPointEYhigh(i, U-N);
   }
   g->Draw("AP");
}
