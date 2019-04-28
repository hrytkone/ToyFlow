#include <TMath.h>
#include <TFile.h>

#include "ResIter.h"

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr);
double CalculateVn(double QnQnA, double QnAQnB, double w);
double CalculateVnError(double QnQnA, double QnAQnB, double QnQnAerr, double QnAQnBerr, double w, double wErr);

//void MakeGraphs(TString sInputName = "test.root", TString sOutputName = "toyFlowGraphs.root") {
void MakeGraphs(TString sInputName = "toyFlow.root", TString sOutputName = "toyFlowGraphs.root") {
//void MakeGraphs(TString sInputName = "events-10000_dneta-1000.root", TString sOutputName = "toyFlowGraphs.root") {
//void MakeGraphs(TString sInputName = "events-15000_dneta-1000_pt-depend.root", TString sOutputName = "toyFlowGraphs.root") {
//void MakeGraphs(TString sInputName = "events-100_dneta-1000_pt-depend_pt-weight.root", TString sOutputName = "toyFlowGraphs.root") {

    TFile *fIn = TFile::Open(sInputName, "read");
    TFile *fOut = TFile::Open(sOutputName, "recreate");

    const int nCoef = 5;
    int nPtBins = 9;

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

    TH1D *hSqrtSumWeights = (TH1D*)fIn->Get("hSqrtSumWeights");
    TH1D *hSqrtSumWeightsA = (TH1D*)fIn->Get("hSqrtSumWeightsA");
    TH1D *hSqrtSumWeightsB = (TH1D*)fIn->Get("hSqrtSumWeightsB");

    TH1D *hSqrtSumWeightsNonuni = (TH1D*)fIn->Get("hSqrtSumWeightsNonuni");
    TH1D *hSqrtSumWeightsANonuni = (TH1D*)fIn->Get("hSqrtSumWeightsANonuni");
    TH1D *hSqrtSumWeightsBNonuni = (TH1D*)fIn->Get("hSqrtSumWeightsBNonuni");

    TH1D *hInputFlow = new TH1D("hInputFlow", "hInputFlow", nCoef, 0.5, double(nCoef)+0.5);
    hInputFlow->SetLineStyle(1);
    hInputFlow->SetLineColor(1);
    hInputFlow->SetLineWidth(2);
    for(int i = 0; i < nCoef; i++)
        hInputFlow->Fill(double(i+1), inputFlow[i]);

    //=====vn=====
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
    TH1D *hRtrueCorrected[nCoef];
    for (i = 0; i < nCoef; i++) {
        hRsub[i] = (TH1D*)fIn->Get(Form("hRsub%02i", i+1));
        hRsubCorrected[i] = (TH1D*)fIn->Get(Form("hRsubCorrected%02i", i+1));
        hRtrue[i] = (TH1D*)fIn->Get(Form("hRtrue%02i", i+1));
        hRtrueCorrected[i] = (TH1D*)fIn->Get(Form("hRtrueCorrected%02i", i+1));
    }

    double vnEP[nCoef], errorVnEP[nCoef] = {0};
    double vnEPCorrected[nCoef], errorVnEPCorrected[nCoef] = {0};
    double vnEPtrue[nCoef], errorVnEPtrue[nCoef] = {0};
    double vnEPCorrectedTrue[nCoef], errorvnEPCorrectedTrue[nCoef] = {0};
    double Rinit;

    double khi0 = 1.0;
    double err = 0.0001;
    double khi;
    double R[nCoef], errorR[nCoef] = {0};
    double RCorrected[nCoef], errorRCorrected[nCoef] = {0};

    for (i = 0; i < nCoef; i++) {
        n = i+1;
        vnEP[i] = hVnObs[i]->GetMean();
        vnEPtrue[i] = hVnObs[i]->GetMean();

        Rinit = hRsub[i]->GetMean();
        khi = RkIter(khi0, Rinit, n, err);
        R[i] = Rk(TMath::Sqrt(2)*khi, n);
        errorR[i] = CalculateRerror(khi, err, n);

        vnEP[i] /= R[i];
        errorVnEP[i] = GetVnError(hVnObs[i]->GetMean(), hVnObs[i]->GetMeanError(), R[i], errorR[i]);

        vnEPtrue[i] /= hRtrue[i]->GetMean();
        errorVnEPtrue[i] = GetVnError(hVnObs[i]->GetMean(), hVnObs[i]->GetMeanError(), hRtrue[i]->GetMean(), hRtrue[i]->GetMeanError());
    }

    for (i=0; i<nCoef; i++) {
        n = i+1;
        vnEPCorrected[i] = hVnObsCorrected[i]->GetMean();
        vnEPCorrectedTrue[i] = hVnObsCorrected[i]->GetMean();

        Rinit = hRsubCorrected[i]->GetMean();
        khi = RkIter(khi0, Rinit, n, err);
        RCorrected[i] = Rk(TMath::Sqrt(2)*khi, n);
        errorRCorrected[i] = CalculateRerror(khi, err, n);

        vnEPCorrected[i] /= RCorrected[i];
        errorVnEPCorrected[i] = GetVnError(hVnObsCorrected[i]->GetMean(), hVnObsCorrected[i]->GetMeanError(), RCorrected[i], errorRCorrected[i]);

        vnEPCorrectedTrue[i] /=hRtrueCorrected[i]->GetMean();
        errorvnEPCorrectedTrue[i] = GetVnError(hVnObsCorrected[i]->GetMean(), hVnObsCorrected[i]->GetMeanError(), hRtrueCorrected[i]->GetMean(), hRtrueCorrected[i]->GetMeanError());
    }

    //vn{EP}
    TH1D *hQnQnA[nCoef];
    TH1D *hQnAQnB[nCoef];
    TH1D *hQnQnAcorrected[nCoef];
    TH1D *hQnAQnBcorrected[nCoef];
    for (i = 0; i < nCoef; i++) {
        hQnQnA[i] = (TH1D*)fIn->Get(Form("hQnQnA%02i", i+1));
        hQnAQnB[i] = (TH1D*)fIn->Get(Form("hQnAQnB%02i", i+1));
        hQnQnAcorrected[i] = (TH1D*)fIn->Get(Form("hQnQnAcorrected%02i", i+1));
        hQnAQnBcorrected[i] = (TH1D*)fIn->Get(Form("hQnAQnBcorrected%02i", i+1));
    }

    double w = hSqrtSumWeights->GetMean();
    double wError = hSqrtSumWeights->GetMeanError();
    double wCorrected = hSqrtSumWeightsNonuni->GetMean();
    double wCorrectedError = hSqrtSumWeightsNonuni->GetMeanError();

    double vnEPalt[nCoef], errorVnEPalt[nCoef] = {0};
    double vnEPCorrectedAlt[nCoef], errorVnEPCorrectedAlt[nCoef] = {0};
    for (i=0; i<nCoef; i++) {
        vnEPalt[i] = CalculateVn(hQnQnA[i]->GetMean(), hQnAQnB[i]->GetMean(), w);
        errorVnEPalt[i] = CalculateVnError(hQnQnA[i]->GetMean(), hQnAQnB[i]->GetMean(), hQnQnA[i]->GetMeanError(), hQnAQnB[i]->GetMeanError(),  w, wError);
        vnEPCorrectedAlt[i] = CalculateVn(hQnQnAcorrected[i]->GetMean(), hQnAQnBcorrected[i]->GetMean(), wCorrected);
        errorVnEPalt[i] = CalculateVnError(hQnQnAcorrected[i]->GetMean(), hQnAQnBcorrected[i]->GetMean(), hQnQnAcorrected[i]->GetMeanError(), hQnAQnBcorrected[i]->GetMeanError(),  wCorrected, wCorrectedError);
    }

    // pT bins
    TH1D *hPtBin[nPtBins];
    TH1D *hSqrtSumWeightsPtBins[nPtBins];
    for (i=0; i<nPtBins; i++) {
        hPtBin[i] = (TH1D*)fIn->Get(Form("hPtBin%02i", i+1));
        hSqrtSumWeightsPtBins[i] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsPtBins%02i", i+1));
    }

    double ptBin[nPtBins];
    double ptBinError[nPtBins];
    for (i=0; i<nPtBins; i++) {
        ptBin[i] = CalculateVn(hPtBin[i]->GetMean(), hQnAQnB[1]->GetMean(), hSqrtSumWeightsPtBins[i]->GetMean());
        ptBinError[i] = CalculateVnError(hPtBin[i]->GetMean(), hQnAQnB[1]->GetMean(), hPtBin[i]->GetMeanError(), hQnAQnB[1]->GetMeanError(), hSqrtSumWeightsPtBins[i]->GetMean(), hSqrtSumWeightsPtBins[i]->GetMeanError());
    }

    // Make graphs
    TGraphErrors *gVnEP = new TGraphErrors(nCoef);
    TGraphErrors *gVnEPCorrected = new TGraphErrors(nCoef);
    TGraphErrors *gVnEPtrue = new TGraphErrors(nCoef);
    TGraphErrors *gVnEPCorrectedTrue = new TGraphErrors(nCoef);
    TGraphErrors *gRtrue = new TGraphErrors(nCoef);
    TGraphErrors *gRtrueCorrected = new TGraphErrors(nCoef);
    TGraphErrors *gR = new TGraphErrors(nCoef);
    TGraphErrors *gRcorrected = new TGraphErrors(nCoef);
    TGraphErrors *gPtBin = new TGraphErrors(6);

    TGraphErrors *gVnAlt = new TGraphErrors(nCoef);
    TGraphErrors *gVnAltCorrected = new TGraphErrors(nCoef);

    for(i = 0; i < nCoef; i++){
        gVnEP->SetPoint(i, double(i+1)-0.1, vnEP[i]);
        gVnEP->SetPointError(i, 0.0, errorVnEP[i]);
        gVnEPCorrected->SetPoint(i, double(i+1)+0.1, vnEPCorrected[i]);
        gVnEPCorrected->SetPointError(i, 0.0, errorVnEPCorrected[i]);
        gVnEPtrue->SetPoint(i, double(i+1), vnEPtrue[i]);
        gVnEPtrue->SetPointError(i, 0.0, errorVnEPtrue[i]);
        gVnEPCorrectedTrue->SetPoint(i, double(i+1)+0.2, vnEPCorrectedTrue[i]);
        gVnEPCorrectedTrue->SetPointError(i, 0.0, errorvnEPCorrectedTrue[i]);
        gRtrue->SetPoint(i, double(i+1)-0.1, hRtrue[i]->GetMean());
        gRtrue->SetPointError(i, 0.0, hRtrue[i]->GetMeanError());
        gRtrueCorrected->SetPoint(i, double(i+1)-0.2, hRtrueCorrected[i]->GetMean());
        gRtrueCorrected->SetPointError(i, 0.0, hRtrueCorrected[i]->GetMeanError());
        gR->SetPoint(i, double(i+1), R[i]);
        gR->SetPointError(i, 0.0, errorR[i]);
        gRcorrected->SetPoint(i, double(i+1)+0.1, RCorrected[i]);
        gRcorrected->SetPointError(i, 0.0, errorRCorrected[i]);

        gVnAlt->SetPoint(i, double(i+1)-0.2, vnEPalt[i]);
        gVnAlt->SetPointError(i, 0.0, errorVnEPalt[i]);
        gVnAltCorrected->SetPoint(i, double(i+1)-0.3, vnEPCorrectedAlt[i]);
        gVnAltCorrected->SetPointError(i, 0.0, errorVnEPCorrectedAlt[i]);
    }

    double point = 0.0;
    double step = 0.0;
    for (i=0; i<nPtBins; i++) {
        step += 0.2;
        point += step;
        gPtBin->SetPoint(i, point-step/2.0, ptBin[i]);
        gPtBin->SetPointError(i, step/2.0, ptBinError[i]);
    }

    fOut->cd();
    hInputFlow->Write("hInputFlow");
    gVnEP->Write("gVnEP");
    gVnEPCorrected->Write("gVnEPCorrected");
    gVnEPtrue->Write("gVnEPtrue");
    gVnEPCorrectedTrue->Write("gVnEPCorrectedTrue");
    gR->Write("gR");
    gRtrue->Write("gRtrue");
    gRcorrected->Write("gRcorrected");
    gRtrueCorrected->Write("gRtrueCorrected");
    gPtBin->Write("gPtBin");
    gVnAlt->Write("gVnAlt");
    gVnAltCorrected->Write("gVnAltCorrected");
    fOut->Close();

}

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr) {
    return (vnObs/Rn)*TMath::Sqrt((vnObsErr/vnObs)*(vnObsErr/vnObs) + (RnErr/Rn)*(RnErr/Rn));
}

double CalculateVn(double QnQnA, double QnAQnB, double w) {
    return QnQnA/TMath::Sqrt(QnAQnB)/w;
}

double CalculateVnError(double QnQnA, double QnAQnB, double QnQnAerr, double QnAQnBerr, double w, double wErr) {
    double vn = CalculateVn(QnQnA, QnAQnB, w);
    return vn*(QnQnAerr/QnQnA + 0.5*QnAQnBerr/QnAQnB + wErr/w);
}
