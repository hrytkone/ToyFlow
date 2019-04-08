#include <TMath.h>
#include <TFile.h>

#include "ResIter.h"

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr);

void MakeGraphs(TString sInputName = "toyFlow.root", TString sOutputName = "toyFlowGraphs.root") {

    TFile *fIn = TFile::Open(sInputName, "read");
    TFile *fOut = TFile::Open(sOutputName, "recreate");

    const int nCoef = 5;
    int i, n;
    double pi = TMath::Pi();
    double inputFlow[nCoef] = {0.0};

    TH1D *hInputNumbers = (TH1D*)fIn->Get("hInputNumbers");
    double nEvents = hInputNumbers->GetBinContent(1);
    double dNdeta = hInputNumbers->GetBinContent(2);
    double etaRange = hInputNumbers->GetBinContent(3);
    double nMult = hInputNumbers->GetBinContent(4);
    for(int i = 0; i < nCoef; i++) inputFlow[i] = hInputNumbers->GetBinContent(5+i);
    double Tdec = hInputNumbers->GetBinContent(10);
    double vr = hInputNumbers->GetBinContent(11);
    double Teff = hInputNumbers->GetBinContent(12);
    double slope = hInputNumbers->GetBinContent(13);

    TH1D *hInputFlow = new TH1D("hInputFlow", "hInputFlow", nCoef, 0.5, double(nCoef)+0.5);
    hInputFlow->SetLineStyle(1);
    hInputFlow->SetLineColor(1);
    hInputFlow->SetLineWidth(2);
    for(int i = 0; i < nCoef; i++)
        hInputFlow->Fill(double(i+1), inputFlow[i]);

    //=====vn{EP}=====
    // Observed vn
    TH1D *hVnObs[nCoef];
    TH1D *hVnObsCorrected[nCoef];
    for (i = 0; i < nCoef; i++) {
        hVnObs[i] = (TH1D*)fIn->Get(Form("hVnObs%02i", i+1));
        hVnObsCorrected[i] = (TH1D*)fIn->Get(Form("hVnObsCorrected%02i", i+1));
    }

    // Resolution parameter
    TH1D *hRsub[nCoef];
    TH1D *hRsubCorrected[nCoef];
    TH1D *hRtrue[nCoef];
    for (i = 0; i < nCoef; i++) {
        hRsub[i] = (TH1D*)fIn->Get(Form("hRsub%02i", i+1));
        hRsubCorrected[i] = (TH1D*)fIn->Get(Form("hRsubCorrected%02i", i+1));
        hRtrue[i] = (TH1D*)fIn->Get(Form("hRtrue%02i", i+1));
    }

    double vnEP[nCoef], errorVnEP[nCoef];
    double vnEPCorrected[nCoef], errorVnEPCorrected[nCoef];
    double vnEPtrue[nCoef], errorVnEPtrue[nCoef];
    double Rinit;

    double khi0 = 1.0;
    double err = 0.0001;
    double khi;
    double R[nCoef], errorR[nCoef];
    double RCorrected[nCoef], errorRCorrected[nCoef];

    for (i = 0; i < nCoef; i++) {
        n = i+1;
        vnEP[i] = hVnObs[i]->GetMean();
        vnEPCorrected[i] = hVnObsCorrected[i]->GetMean();
        vnEPtrue[i] = hVnObs[i]->GetMean();

        Rinit = hRsub[i]->GetMean();
        khi = RkIter(khi0, Rinit, n, err);
        R[i] = Rk(TMath::Sqrt(2)*khi, n);
        errorR[i] = CalculateRerror(khi, err, n);

        vnEP[i] /= R[i];
        errorVnEP[i] = GetVnError(hVnObs[i]->GetMean(), hVnObs[i]->GetMeanError(), R[i], errorR[i]);

        vnEPtrue[i] /= hRtrue[i]->GetMean();
        errorVnEPtrue[i] = GetVnError(hVnObs[i]->GetMean(), hVnObs[i]->GetMeanError(), hRtrue[i]->GetMean(), hRtrue[i]->GetMeanError());

        Rinit = hRsubCorrected[i]->GetMean();
        khi = RkIter(khi0, Rinit, n, err);
        RCorrected[i] = Rk(TMath::Sqrt(2)*khi, n);
        errorRCorrected[i] = CalculateRerror(khi, err, n);

        vnEPCorrected[i] /= RCorrected[i];
        errorVnEPCorrected[i] = GetVnError(hVnObsCorrected[i]->GetMean(), hVnObsCorrected[i]->GetMeanError(), RCorrected[i], errorRCorrected[i]);
    }

    // pT bins
    TH1D *hPtBin[6];
    for (i=0; i<6; i++) {
        hPtBin[i] = (TH1D*)fIn->Get(Form("hPtBin%02i", i+1));
    }

    double ptBin[6];
    double ptBinError[6];
    for (i=0; i<6; i++) {
        ptBin[i] = hPtBin[i]->GetMean();
        ptBinError[i] = hPtBin[i]->GetMeanError();
    }

    // Make graphs
    TGraphErrors *gVnEP = new TGraphErrors(nCoef);
    TGraphErrors *gVnEPCorrected = new TGraphErrors(nCoef);
    TGraphErrors *gVnEPtrue = new TGraphErrors(nCoef);
    TGraphErrors *gRtrue = new TGraphErrors(nCoef);
    TGraphErrors *gR = new TGraphErrors(nCoef);
    TGraphErrors *gRcorrected = new TGraphErrors(nCoef);
    TGraphErrors *gPtBin = new TGraphErrors(6);

    for(i = 0; i < nCoef; i++){
        gVnEP->SetPoint(i, double(i+1)-0.15, vnEP[i]);
        gVnEP->SetPointError(i, 0.0, errorVnEP[i]);
        gVnEPCorrected->SetPoint(i, double(i+1)+0.15, vnEPCorrected[i]);
        gVnEPCorrected->SetPointError(i, 0.0, errorVnEPCorrected[i]);
        gVnEPtrue->SetPoint(i, double(i+1), vnEPtrue[i]);
        gVnEPtrue->SetPointError(i, 0.0, errorVnEPtrue[i]);
        gRtrue->SetPoint(i, double(i+1)-0.15, hRtrue[i]->GetMean());
        gRtrue->SetPointError(i, 0.0, hRtrue[i]->GetMeanError());
        gR->SetPoint(i, double(i+1), R[i]);
        gR->SetPointError(i, 0.0, errorR[i]);
        gRcorrected->SetPoint(i, double(i+1)+0.15, RCorrected[i]);
        gRcorrected->SetPointError(i, 0.0, errorRCorrected[i]);
    }

    for (i=0; i<6; i++) {
        gPtBin->SetPoint(i, double(i)+0.5, ptBin[i]);
        gPtBin->SetPointError(i, 0.0, ptBinError[i]);
    }

    fOut->cd();
    hInputFlow->Write("hInputFlow");
    gVnEP->Write("gVnEP");
    gVnEPCorrected->Write("gVnEPCorrected");
    gVnEPtrue->Write("gVnEPtrue");
    gR->Write("gR");
    gRtrue->Write("gRtrue");
    gRcorrected->Write("gRcorrected");
    gPtBin->Write("gPtBin");
    fOut->Close();
}

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr) {
    return (vnObs/Rn)*TMath::Sqrt((vnObsErr/vnObs)*(vnObsErr/vnObs) + (RnErr/Rn)*(RnErr/Rn));
}
