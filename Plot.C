#include <TMath.h>
#include <TFile.h>

double fFitMeanQ( double *x, double *p );

void Plot() {

    double pi = TMath::Pi();

    int i;
    const int n = 5;

    TString sFileName = "toyFlow_n-50000_uni.root";
    TFile *fIn = TFile::Open(sFileName, "read");

    TH1D *hPhi = (TH1D*)fIn->Get("hPhi");
    TH1D *hPhiNonuni = (TH1D*)fIn->Get("hPhiNonuni");

    // Q-vector for uniform phi-distribution
    TH1D *hQx[nCoef];
    TH1D *hQy[nCoef];
    TH1D *hQxNonuni[nCoef];
    TH1D *hQyNonuni[nCoef];
    TH1D *hQxRescaled[nCoef];
    TH1D *hQyRescaled[nCoef];
    TH1D *hVnObsUni[nCoef];
    TH1D *hVnObsNonuni[nCoef];
    for (i=0; i<nCoef; i++) {
        hQx[i] = (TH1D*)fIn->Get(Form("hQx%02i", i+1));
        hQy[i] = (TH1D*)fIn->Get(Form("hQy%02i", i+1));
        hQxNonuni[i] = (TH1D*)fIn->Get(Form("hQxNonuni%02i", i+1));
        hQyNonuni[i] = (TH1D*)fIn->Get(Form("hQyNonuni%02i", i+1));
        hQxRescaled[i] = (TH1D*)fIn->Get(Form("hQxRescaled%02i", i+1));
        hQyRescaled[i] = (TH1D*)fIn->Get(Form("hQyRescaled%02i", i+1));
        hVnObsUni[i] = (TH1D*)fIn->Get(Form("hVnObsUni%02i", i+1));
        hVnObsNonuni[i] = (TH1D*)fIn->Get(Form("hVnObsNonuni%02i", i+1));
    }

    TF1 *fitQ = new TF1("fitQ", fFitMeanQ, 0.0, 15.0, 3);

    for(int i = 0; i < nCoef; i++){

        fitQ->SetParameters(1.0, hQ[i]->GetBinCenter( hQ[i]->GetMaximumBin() ), 1.0);
        fitQ->FixParameter(0,1.0);
        hQ[i]->Fit("fitQ","RNO");
        meanQ[i]  = TMath::Max(0.0,fitQ->GetParameter(1));  meanErrorQ[i]  = fitQ->GetParError(1);

        gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
        //gPad->SetLogx(0); gPad->SetLogy(0);
        //mpad->SetGridx(0); mpad->SetGridy(0);

        hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, 0.0, 12.0, 10, 0.0, 0.9);
        //hset( *hfr, "Q", "1/N_{entries} dN/dQ",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
        hfr->Draw();

        hQ[i]->Draw("same,p");
        fitQ->DrawCopy("same,l");

        leg = new TLegend(0.50, 0.68, 0.80, 0.93,"","brNDC");
        leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
        leg->AddEntry(hQ[i],"toy simulation","p");
        leg->AddEntry(fitQ,"fit, resulting parameters:","l");
        leg->AddEntry(fitQ,Form("<Q> = %4.3f +- %4.3f",fitQ->GetParameter(1), TMath::Max(0.001,fitQ->GetParError(1))),"");
        leg->AddEntry(fitQ,Form("#sigma = %4.3f +- %4.3f",fitQ->GetParameter(2), TMath::Max(0.001,fitQ->GetParError(2))),"");
        leg->Draw();
    }
}

double fFitMeanQ( double *x, double *p ) {
    double Q = x[0];
    double norm = p[0]; // should be order 1
    double meanQ = p[1];
    double width = p[2];
    double term1 = norm*2.*Q/width/width;
    double term2 = TMath::Exp( -1. * ( Q*Q + meanQ*meanQ ) / width/width );
    double term3 = TMath::BesselI0( 2.*Q*meanQ/width/width );
    return term1*term2*term3;
}
