TCanvas a("Acceptance","Acceptance");
Double_t x[12]={0.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,15.0,20.0,30.0,50.0};
TH1D Acc("Acceptance","Acceptance",11,x);  
TH1D Events("Events","Events",11,x);
TH1D LambdaB_pt("LambdaB_pt", "LambdaB_pt",25,0.,50.);
TFile f("/shome/kaestli/CMSSW_3_9_7/src/reducedTree_MSEL1-nocuts.root");
TTree *t=(TTree* ) f.Get("events");
TString AcceptanceCuts;

void acceptance02(){
  AcceptanceCuts="((mu1eta<1.3 && mu1pT>3.3)||(mu1eta>1.3 && mu1eta<2.2 && mu1p>2.9)||(mu1eta>2.2 && mu1eta<2.4 && mu1pT>0.8))";
  AcceptanceCuts+=" && ((mu2eta<1.3 && mu2pT>3.3)||(mu2eta>1.3 && mu2eta<2.2 && mu2p>2.9)||(mu2eta>2.2 && mu2eta<2.4 && mu2pT>0.8))";
  //  AcceptanceCuts+=" && (p_pT>1 && p_eta<2.3 && pi_eta<2.3 && pi_pT>0.4) ";
  AcceptanceCuts += " && sqrt(p_Vx*p_Vx+p_Vy*p_Vy)>1.0 && sqrt(p_Vx*p_Vx+p_Vy*p_Vy)<35.0 && p_Vz<100.0";
  Double_t n_nocuts=0;
  Double_t n_cuts=0;
  stringstream buffer;
  TString pt_cut;
  Double_t acceptance;
  Int_t events;

  t->Draw("lambdaB_pt>>LambdaB_pt",AcceptanceCuts);
  for(Int_t i=1; i<12; i++){
    buffer.str("");
    buffer << "lambdaB_pt > "<<x[i-1] << " && lambdaB_pt < "<<x[i]; 
    pt_cut = TString(buffer.str());  
    n_nocuts = (Double_t) t->GetEntries(pt_cut);
    n_cuts = (Double_t) t->GetEntries(AcceptanceCuts + "&& "+pt_cut);
    if(n_nocuts>0) acceptance = n_cuts/n_nocuts;
    else acceptance =0.0;
    Acc.SetBinContent(i,acceptance);
    Events.SetBinContent(i,n_cuts/(x[i]-x[i-1]));
  }
  Acc.Draw();
}
