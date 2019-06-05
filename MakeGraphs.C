#include <TMath.h>
#include <TFile.h>

#include "src/ResIter.h"
#include "src/JConst.h"

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr);
double CalculateVn(double QnQnA, double QnAQnB, double w);
double CalculateVnError(double QnQnA, double QnAQnB, double QnQnAerr, double QnAQnBerr, double w, double wErr);

void MakeGraphs(TString sInputName = "toyFlow.root", TString sOutputName = "toyFlowGraphs.root", const int iDet=0) {

    TFile *fIn = TFile::Open(sInputName, "read");
    TFile *fOut = TFile::Open(sOutputName, "recreate");

    int nPtBins = 9;

    int i, n;
    double pi = TMath::Pi();
    double inputFlow[nCoef] = {0.0};

    TH1D *hInputNumbers = (TH1D*)fIn->Get("hInputNumbers");
    double nEvents = hInputNumbers->GetBinContent(1);
    double nOfFiles = hInputNumbers->GetBinContent(14);
    for(int i = 0; i < nCoef; i++) {
        inputFlow[i] = hInputNumbers->GetBinContent(2+i);
        inputFlow[i] /= nOfFiles;
    }

    TH1D *hSqrtSumWeights[DET_N];
    TH1D *hSqrtSumWeightsNonuni[DET_N];

    hSqrtSumWeights[iDet] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsD%02i",iDet));
    hSqrtSumWeightsNonuni[iDet] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsNonuniD%02i",iDet));

    TH1D *hInputFlow = new TH1D("hInputFlow", "hInputFlow", nCoef, 0.5, double(nCoef)+0.5);
    hInputFlow->SetLineStyle(1);
    hInputFlow->SetLineColor(1);
    hInputFlow->SetLineWidth(1);
    for(int i = 0; i < nCoef; i++)
        hInputFlow->Fill(double(i+1), inputFlow[i]);
    hInputFlow->Fill(double(nCoef+1), 0.0);

    TH1D *hPhi = (TH1D*)fIn->Get("hPhi");
    TH1D *hPhiNonuni = (TH1D*)fIn->Get("hPhiNonuni");

    //=====vn=====
    // Observed vn
    TH1D *hVnObs[nCoef][DET_N];
    TH1D *hVnObsCorrected[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hVnObs[i][iDet] = (TH1D*)fIn->Get(Form("hVnObsH%02iD%02i", i+1,iDet));
        hVnObsCorrected[i][iDet] = (TH1D*)fIn->Get(Form("hVnObsCorrectedH%02iD%02i", i+1,iDet));
    }

    // Resolution parameter
    TH1D *hRsub[nCoef][DET_N];
    TH1D *hRsubCorrected[nCoef][DET_N];
    TH1D *hRtrue[nCoef][DET_N];
    TH1D *hRtrueCorrected[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hRsub[i][iDet] = (TH1D*)fIn->Get(Form("hRsubH%02iD%02i", i+1,iDet));
        hRsubCorrected[i][iDet] = (TH1D*)fIn->Get(Form("hRsubCorrectedH%02iD%02i", i+1,iDet));
        hRtrue[i][iDet] = (TH1D*)fIn->Get(Form("hRtrueH%02iD%02i", i+1,iDet));
        hRtrueCorrected[i][iDet] = (TH1D*)fIn->Get(Form("hRtrueCorrectedH%02iD%02i", i+1,iDet));
    }

    double vnEP[nCoef], errorVnEP[nCoef] = {0};
    double vnEPCorrected[nCoef], errorVnEPCorrected[nCoef] = {0};
    double vnEPtrue[nCoef], errorVnEPtrue[nCoef] = {0};
    double vnEPCorrectedTrue[nCoef], errorvnEPCorrectedTrue[nCoef] = {0};
    double Rinit;

    double khi0 = 0.5;
    double err = 0.0001;
    double khi;
    double R[nCoef], errorR[nCoef] = {0};
    double RCorrected[nCoef], errorRCorrected[nCoef] = {0};

    cout << "Uniform:\n";
    for (i = 0; i < nCoef; i++) {
        n = i+1;
        vnEP[i] = hVnObs[i][iDet]->GetMean();
        vnEPtrue[i] = hVnObs[i][iDet]->GetMean();

        Rinit = hRsub[i][iDet]->GetMean();
        khi = RIter(khi0, Rinit, err);
        R[i] = R1(TMath::Sqrt(2)*khi);
        errorR[i] = CalculateRerror(khi, err);

        cout << "R=" << R[i] << "  err=" << errorR[i] << "\n";

        vnEP[i] /= R[i];
        errorVnEP[i] = GetVnError(hVnObs[i][iDet]->GetMean(), hVnObs[i][iDet]->GetMeanError(), R[i], errorR[i]);

        vnEPtrue[i] /= hRtrue[i][iDet]->GetMean();
        errorVnEPtrue[i] = GetVnError(hVnObs[i][iDet]->GetMean(), hVnObs[i][iDet]->GetMeanError(), hRtrue[i][iDet]->GetMean(), hRtrue[i][iDet]->GetMeanError());
    }

    cout << "\nNon-uniform:\n";
    for (i=0; i<nCoef; i++) {
        n = i+1;
        vnEPCorrected[i] = hVnObsCorrected[i][iDet]->GetMean();
        vnEPCorrectedTrue[i] = hVnObsCorrected[i][iDet]->GetMean();

        Rinit = hRsubCorrected[i][iDet]->GetMean();
        khi = RIter(khi0, Rinit, err);
        RCorrected[i] = R1(TMath::Sqrt(2)*khi);
        errorRCorrected[i] = CalculateRerror(khi, err);

        cout << "R=" << RCorrected[i] << "  err=" << errorRCorrected[i] << "\n";

        vnEPCorrected[i] /= RCorrected[i];
        if (RCorrected[i]==0) vnEPCorrected[i] = 0.0;
        errorVnEPCorrected[i] = GetVnError(hVnObsCorrected[i][iDet]->GetMean(), hVnObsCorrected[i][iDet]->GetMeanError(), RCorrected[i], errorRCorrected[i]);

        vnEPCorrectedTrue[i] /=hRtrueCorrected[i][iDet]->GetMean();
        errorvnEPCorrectedTrue[i] = GetVnError(hVnObsCorrected[i][iDet]->GetMean(), hVnObsCorrected[i][iDet]->GetMeanError(), hRtrueCorrected[i][iDet]->GetMean(), hRtrueCorrected[i][iDet]->GetMeanError());
    }

    //vn{EP}
    TH1D *hQnQnA[nCoef][DET_N];
    TH1D *hQnAQnB[nCoef][DET_N];
    TH1D *hQnQnAcorrected[nCoef][DET_N];
    TH1D *hQnAQnBcorrected[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hQnQnA[i][iDet] = (TH1D*)fIn->Get(Form("hQnQnAH%02iD%02i", i+1,iDet));
        hQnAQnB[i][iDet] = (TH1D*)fIn->Get(Form("hQnAQnBH%02iD%02i", i+1,iDet));
        hQnQnAcorrected[i][iDet] = (TH1D*)fIn->Get(Form("hQnQnAcorrectedH%02iD%02i", i+1,iDet));
        hQnAQnBcorrected[i][iDet] = (TH1D*)fIn->Get(Form("hQnAQnBcorrectedH%02iD%02i", i+1,iDet));
    }

    double w = hSqrtSumWeights[iDet]->GetMean();
    double wError = hSqrtSumWeights[iDet]->GetMeanError();
    double wCorrected = hSqrtSumWeightsNonuni[iDet]->GetMean();
    double wCorrectedError = hSqrtSumWeightsNonuni[iDet]->GetMeanError();

    double vnEPalt[nCoef], errorVnEPalt[nCoef] = {0};
    double vnEPCorrectedAlt[nCoef], errorVnEPCorrectedAlt[nCoef] = {0};
    for (i=0; i<nCoef; i++) {
        vnEPalt[i] = CalculateVn(hQnQnA[i][iDet]->GetMean(), hQnAQnB[i][iDet]->GetMean(), w);
        errorVnEPalt[i] = CalculateVnError(hQnQnA[i][iDet]->GetMean(), hQnAQnB[i][iDet]->GetMean(), hQnQnA[i][iDet]->GetMeanError(), hQnAQnB[i][iDet]->GetMeanError(),  w, wError);

        vnEPCorrectedAlt[i] = CalculateVn(hQnQnAcorrected[i][iDet]->GetMean(), hQnAQnBcorrected[i][iDet]->GetMean(), wCorrected);
        errorVnEPCorrectedAlt[i] = CalculateVnError(hQnQnAcorrected[i][iDet]->GetMean(), hQnAQnBcorrected[i][iDet]->GetMean(), hQnQnAcorrected[i][iDet]->GetMeanError(), hQnAQnBcorrected[i][iDet]->GetMeanError(),  wCorrected, wCorrectedError);
    }

    // pT bins
    TH1D *hPtBin[nPtBins];
    TH1D *hSqrtSumWeightsPtBins[nPtBins];
    for (i=0; i<nPtBins; i++) {
        hPtBin[i] = (TH1D*)fIn->Get(Form("hPtBinH%02i", i+1));
        hSqrtSumWeightsPtBins[i] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsPtBinsH%02i", i+1));
    }

    double ptBin[nPtBins];
    double ptBinError[nPtBins];
    for (i=0; i<nPtBins; i++) {
        ptBin[i] = CalculateVn(hPtBin[i]->GetMean(), hQnAQnB[1][iDet]->GetMean(), hSqrtSumWeightsPtBins[i]->GetMean());
        ptBinError[i] = CalculateVnError(hPtBin[i]->GetMean(), hQnAQnB[1][iDet]->GetMean(), hPtBin[i]->GetMeanError(), hQnAQnB[1][iDet]->GetMeanError(), hSqrtSumWeightsPtBins[i]->GetMean(), hSqrtSumWeightsPtBins[i]->GetMeanError());
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
        gVnEPtrue->SetPoint(i, double(i+1)-0.25, vnEPtrue[i]);
        gVnEPtrue->SetPointError(i, 0.0, errorVnEPtrue[i]);
        gVnEPCorrectedTrue->SetPoint(i, double(i+1)-0.15, vnEPCorrectedTrue[i]);
        gVnEPCorrectedTrue->SetPointError(i, 0.0, errorvnEPCorrectedTrue[i]);
        gVnEP->SetPoint(i, double(i+1)-0.05, vnEP[i]);
        gVnEP->SetPointError(i, 0.0, errorVnEP[i]);
        gVnEPCorrected->SetPoint(i, double(i+1)+0.05, vnEPCorrected[i]);
        gVnEPCorrected->SetPointError(i, 0.0, errorVnEPCorrected[i]);
        gVnAlt->SetPoint(i, double(i+1)+0.15, vnEPalt[i]);
        gVnAlt->SetPointError(i, 0.0, errorVnEPalt[i]);
        gVnAltCorrected->SetPoint(i, double(i+1)+0.25, vnEPCorrectedAlt[i]);
        gVnAltCorrected->SetPointError(i, 0.0, errorVnEPCorrectedAlt[i]);

        gRtrue->SetPoint(i, double(i+1)-0.15, hRtrue[i][iDet]->GetMean());
        gRtrue->SetPointError(i, 0.0, hRtrue[i][iDet]->GetMeanError());
        gRtrueCorrected->SetPoint(i, double(i+1)-0.05, hRtrueCorrected[i][iDet]->GetMean());
        gRtrueCorrected->SetPointError(i, 0.0, hRtrueCorrected[i][iDet]->GetMeanError());
        gR->SetPoint(i, double(i+1)+0.05, R[i]);
        gR->SetPointError(i, 0.0, errorR[i]);
        gRcorrected->SetPoint(i, double(i+1)+0.15, RCorrected[i]);
        gRcorrected->SetPointError(i, 0.0, errorRCorrected[i]);
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
    hPhi->Write("hPhi");
    hPhiNonuni->Write("hPhiNonuni");
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
