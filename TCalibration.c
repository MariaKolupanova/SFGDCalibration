Int_t npeaks = 30;

Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+2]; // "height" or "area"
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
#if defined(__PEAKS_C_FIT_AREAS__)
      norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif /* defined(__PEAKS_C_FIT_AREAS__) */
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}



void peaks(){
TFile *file = new TFile("Channels_signal.root");
TH1F *h = (TH1F*)file->Get("Channel_0");
   // Generate n peaks at random
Double_t par[3000];
Int_t p;
TCanvas *c = new TCanvas("c1","c1",10,10,1000,900);
c->Divide(1,2);
c->cd(1);
// c->SetLogy();
TSpectrum *s = new TSpectrum(2*npeaks);

Int_t nfound = s->Search(h,2,"",0.001);
TH1 *hb = s->Background(h,20,"same");
c->Update();
npeaks = 0;

TF1 *fline = new TF1("fline","pol1",0,1000);
// h->Fit("fline","qn");
par[0] = fline->GetParameter(0);
par[1] = fline->GetParameter(1);
Double_t *xpeaks = s->GetPositionX();
for (auto p=0;p<nfound;p++) {
      Double_t xp = xpeaks[p];
      Int_t bin = h->GetXaxis()->FindBin(xp);
      Double_t yp = h->GetBinContent(bin);
      if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
      par[3*p+2] = yp; // "height"
      par[3*p+3] = xp; // "mean"
      par[3*p+4] = 3; // "sigma"
#if defined(__PEAKS_C_FIT_AREAS__)
      par[3*p+2] *= par[3*p+4] * (TMath::Sqrt(TMath::TwoPi())); // "area"
#endif /* defined(__PEAKS_C_FIT_AREAS__) */
      npeaks++;
      }

   // We may have more than the default 25 parameters

TF1 *fit_0 = new TF1("fit_0",fpeaks,par[3]-3*par[4],par[3]+3*par[4],2+3*npeaks);
fit_0->SetParameters(par);
fit_0->SetNpx(1000);
h->Fit("fit_0","r");
c->Update();
Double_t par_error[nfound];
for (auto p=0;p<nfound;p++) {
      TF1 *fit_1 = new TF1("fit_1",fpeaks,par[3*p+3]-2*par[3*p+4],par[3*p+3]+2*par[3*p+4],2+3*nfound);
      fit_1->SetParameters(par);
      fit_1->SetNpx(1000);
      Double_t error = fit_1->GetParError(3*p+3);
      par_error[p] = error;
      cout << par_error[p];
      h->Fit("fit_1","r+");
       c->Update();
      }


c->cd(2);
int n = 5;
Double_t gain[6] = {100.356633,141.298146,179.243937,218.188302,255.135521,294.079886};
Double_t x[5]= {1,2,3,4,5};
Double_t gain_error[6]= {0.035817,0.0020664,0.065720,0.194676,0.320407,1.828337};
Double_t ex[5] = {0,0,0,0,0};
Double_t y[5];
Double_t ey[5];
 for(int i = 1; i < (n+1); i ++){
       y[i-1] = gain[i] - gain[i-1];
       ey[i-1] = gain_error[i]+gain_error[i-1];
//       
//       //cout << ey[i-1]<<" ";
}
TGraphErrors *gr = new TGraphErrors(n,x,y,ex,ey);
      gr->GetXaxis()->SetTitle("Number of difference");
      gr->GetYaxis()->SetTitle("Difference");
      gr->SetMarkerColor(4);
      //gr->SetMarkerStyle(21);
      gr->Draw();
      c->Update();
      gr->Write();
    
}

