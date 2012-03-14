void overflow()
{
  c1 = new TCanvas("c1","Overflow",200,10,700,700);
c1->Divide(1,2);

  hpx    = new TH1F("hpx","This is the px distribution",50,-4,4);

  gRandom->SetSeed();
  Float_t px,py;
  for (Int_t i = 0; i < 25000; i++) {
     gRandom->Rannor(px,py);
     Float_t random = gRandom->Rndm(1);
     hpx->Fill(px);
     hpx->Fill(px+10,0.01);
     hpx->Fill(px-10,0.01);
  }
  gStyle->SetOptStat(111111);

c1->cd(1);
hpx->Draw();

c1->cd(2);
  PaintOverflow(hpx);
}

void PaintOverflow(TH1 *h)
{
   // This function paint the histogram h with an extra bin for overflows

   char* name  = h->GetName();
   char* title = h->GetTitle();
   Int_t nx    = h->GetNbinsX()+1;
   Double_t x1 = h->GetBinLowEdge(1);
   Double_t bw = h->GetBinWidth(nx);
   Double_t x2 = h->GetBinLowEdge(nx)+bw;

   // Book a temporary histogram having ab extra bin for overflows
   TH1F *htmp = new TH1F(name, title, nx, x1, x2);

   // Fill the new hitogram including the extra bin for overflows
   for (Int_t i=1; i<=nx; i++) {
      htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
   }

   // Fill the underflows
   htmp->Fill(x1-1, h->GetBinContent(0));

   // Restore the number of entries
   htmp->SetEntries(h->GetEntries());

   // Draw the temporary histogram
   htmp->Draw();
   TText *t = new TText(x2-bw/2,h->GetBinContent(nx),"Overflow");
   t->SetTextAngle(90);
   t->SetTextAlign(12);
   t->SetTextSize(0.03);;
   t->Draw();
}

