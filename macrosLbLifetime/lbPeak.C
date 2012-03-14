.L setStage.C+

setStage("mlb","(id1m&4)==4&&(id2m&4)==4&mjp>3.03&&mjp<3.16&&ptpr>ptpi&&alphal0<0.002&&ml0<1.15&&ptl0>4.0&&d3l0>8&&d3lb/d3Elb>1.",20,5,6,false);

.L fitgaus.C+

fitgaus((TH1F*)gDirectory->GetList()->FindObject("hdat033"));
((TH1F*)gDirectory->GetList()->FindObject("hdat033"))->GetXaxis()->SetTitle("m(#Lambda^{0}J/#psi) / GeV/c^{2}");
((TH1F*)gDirectory->GetList()->FindObject("hdat033"))->GetYaxis()->SetTitle("evts per bin of 0.05 GeV/c^{2}");
((TH1F*)gDirectory->GetList()->FindObject("hdat033"))->SetTitle("#Lambda_{b} in 12.1/pb");
((TH1F*)gDirectory->GetList()->FindObject("hdat033"))->GetXaxis()->SetNdivisions(502);

