#include "setTDRStyle_modified.C"

void quickPlotCtauVsVrrs () {
   setTDRStyle();
   c1 = new TCanvas("c1","c1",200,10,700,500);
   const Int_t n = 10;

   Double_t x[n]   = {0.375, 1.125, 2, 3.15, 4.65, 6.65, 9.5, 13.95, 22.6, 34.25};

   Double_t y[n]   = {1.466, 1.549, 1.545, 1.543, 1.541, 1.542, 1.547, 1.544, 1.539, 1.540};
   Double_t exl[n] = {0.375, 0.375, 0.5, 0.65, 0.85, 1.15, 1.7, 2.75, 5.9, 5.75};
   Double_t exh[n] = {0.375, 0.375, 0.5, 0.65, 0.85, 1.15, 1.7, 2.75, 5.9, 5.75};
   Double_t eyl[n] = {.003,.003,.003,.003,.003,.003,.003,.003,.003,.003};
   Double_t eyh[n] = {.003,.003,.003,.003,.003,.003,.003,.003,.003,.003};
   gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
   gr->SetTitle("");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");
   gr->GetXaxis()->SetTitle("r(V^{0}) [cm]");
   gr->GetYaxis()->SetTitle("#tau_{truth} [ps]");
   TLine *l = new TLine(0,1.537,42,1.537);
   l->Draw();
   l->SetLineColor(2);
   c1->SaveAs("Lb_vrrs_bins.pdf");
}

