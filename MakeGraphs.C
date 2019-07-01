#include <TMath.h>
#include <TFile.h>

#include "src/ResIter.h"
#include "src/JConst.h"

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr);
double CalculateVn(double QnQnA, double QnAQnB, double w);
double CalculateVnError(double QnQnA, double QnAQnB, double QnQnAerr, double QnAQnBerr, double w, double wErr);
void checkUnderOverFlow( TH1 *h );

void MakeGraphs(TString sInputName = "toyFlow.root", TString sOutputName = "toyFlowGraphs.root", const int iDet=0) {

    TFile *fIn = TFile::Open(sInputName, "read");
    TFile *fOut = TFile::Open(sOutputName, "recreate");

    int nPtBins = 9;

    int i, n;
    double pi = TMath::Pi();
    double inputFlow[nCoef] = {0.0};

    TH1D *hInputNumbers = (TH1D*)fIn->Get("hInputNumbers");
    checkUnderOverFlow(hInputNumbers);
    double nEvents = hInputNumbers->GetBinContent(1);
    double nOfFiles = hInputNumbers->GetBinContent(14);
    for(int i = 0; i < nCoef; i++) {
        inputFlow[i] = hInputNumbers->GetBinContent(2+i);
        inputFlow[i] /= nOfFiles;
    }

    TH1D *hSqrtSumWeights[DET_N];
    TH1D *hSqrtSumWeightsNonuni[DET_N];

    hSqrtSumWeights[iDet] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsD%02i",iDet));
    checkUnderOverFlow(hSqrtSumWeights[iDet]);
    hSqrtSumWeightsNonuni[iDet] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsNonuniD%02i",iDet));
    checkUnderOverFlow(hSqrtSumWeightsNonuni[iDet]);

    TH1D *hInputFlow = new TH1D("hInputFlow", "hInputFlow", nCoef, 0.5, double(nCoef)+0.5);
    hInputFlow->SetLineStyle(1);
    hInputFlow->SetLineColor(1);
    hInputFlow->SetLineWidth(1);
    for(int i = 0; i < nCoef; i++)
        hInputFlow->Fill(double(i+1), inputFlow[i]);
    hInputFlow->Fill(double(nCoef+1), 0.0);

    TH1D *hPhi = (TH1D*)fIn->Get("hPhi");
    checkUnderOverFlow(hPhi);
    TH1D *hPhiNonuni = (TH1D*)fIn->Get("hPhiNonuni");
    checkUnderOverFlow(hPhiNonuni);

    //=====vn=====
    // Observed vn
    TH1D *hVnObs[nCoef][DET_N];
    TH1D *hVnObsNonuni[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hVnObs[i][iDet] = (TH1D*)fIn->Get(Form("hVnObsH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hVnObs[i][iDet]);
        hVnObsNonuni[i][iDet] = (TH1D*)fIn->Get(Form("hVnObsNonuniH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hVnObsNonuni[i][iDet]);
    }

    // Resolution parameter
    TH1D *hRsub[nCoef][DET_N];
    TH1D *hRsubNonuni[nCoef][DET_N];
    TH1D *hRtrue[nCoef][DET_N];
    TH1D *hRtrueNonuni[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hRsub[i][iDet] = (TH1D*)fIn->Get(Form("hRsubH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hRsub[i][iDet]);
        hRsubNonuni[i][iDet] = (TH1D*)fIn->Get(Form("hRsubNonuniH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hRsubNonuni[i][iDet]);
        hRtrue[i][iDet] = (TH1D*)fIn->Get(Form("hRtrueH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hRtrue[i][iDet]);
        hRtrueNonuni[i][iDet] = (TH1D*)fIn->Get(Form("hRtrueNonuniH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hRtrueNonuni[i][iDet]);
    }

    double vn[nCoef], errorVn[nCoef] = {0};
    double vnNonuni[nCoef], errorVnNonuni[nCoef] = {0};
    double vnTrue[nCoef], errorVnTrue[nCoef] = {0};
    double vnTrueNonuni[nCoef], errorVnTrueNonuni[nCoef] = {0};
    double Rinit;

    double khi0 = 0.5;
    double err = 0.0001;
    double khi;
    double R[nCoef], errorR[nCoef] = {0};
    double Rnonuni[nCoef], errorRnonuni[nCoef] = {0};

    cout << "Uniform:\n";
    for (i = 0; i < nCoef; i++) {
        n = i+1;
        vn[i] = hVnObs[i][iDet]->GetMean();
        vnTrue[i] = hVnObs[i][iDet]->GetMean();

        Rinit = TMath::Sqrt(hRsub[i][iDet]->GetMean());
        khi = RIter(khi0, Rinit, err);
        R[i] = R1(TMath::Sqrt(2)*khi); //Because khi sim sqrt(Multi) and full event has twice multiplicity compared to A or B.
        errorR[i] = CalculateRerror(khi, err);

        cout << "R=    " << R[i] <<                       "  err=" << errorR[i] << "\n";
        cout << "Rtrue=" << hRtrue[i][iDet]->GetMean() << "  err=" << hRtrue[i][iDet]->GetMeanError() << "\n";
        cout << "R/Rtrue=" << R[i]/hRtrue[i][iDet]->GetMean() << endl;

        vn[i] /= R[i];
        errorVn[i] = GetVnError(hVnObs[i][iDet]->GetMean(), hVnObs[i][iDet]->GetMeanError(), R[i], errorR[i]);

        vnTrue[i] /= hRtrue[i][iDet]->GetMean();
        errorVnTrue[i] = GetVnError(hVnObs[i][iDet]->GetMean(), hVnObs[i][iDet]->GetMeanError(), hRtrue[i][iDet]->GetMean(), hRtrue[i][iDet]->GetMeanError());
    }

    cout << "\nNon-uniform:\n";
    for (i=0; i<nCoef; i++) {
        n = i+1;
        vnNonuni[i] = hVnObsNonuni[i][iDet]->GetMean();
        vnTrueNonuni[i] = hVnObsNonuni[i][iDet]->GetMean();

        Rinit = hRsubNonuni[i][iDet]->GetMean();
        khi = RIter(khi0, Rinit, err);
        Rnonuni[i] = R1(TMath::Sqrt(2)*khi);
        errorRnonuni[i] = CalculateRerror(khi, err);

        cout << "R=" << Rnonuni[i] << "  err=" << errorRnonuni[i] << "\n";

        vnNonuni[i] /= Rnonuni[i];
        if (Rnonuni[i]==0) vnNonuni[i] = 0.0;
        errorVnNonuni[i] = GetVnError(hVnObsNonuni[i][iDet]->GetMean(), hVnObsNonuni[i][iDet]->GetMeanError(), Rnonuni[i], errorRnonuni[i]);

        vnTrueNonuni[i] /=hRtrueNonuni[i][iDet]->GetMean();
        errorVnTrueNonuni[i] = GetVnError(hVnObsNonuni[i][iDet]->GetMean(), hVnObsNonuni[i][iDet]->GetMeanError(), hRtrueNonuni[i][iDet]->GetMean(), hRtrueNonuni[i][iDet]->GetMeanError());
    }

    //vn{EP} and vn{SP}
    TH1D *hQnQnAEP[nCoef][DET_N];
    TH1D *hQnAQnBEP[nCoef][DET_N];
    TH1D *hQnQnASP[nCoef][DET_N];
    TH1D *hQnAQnBSP[nCoef][DET_N];

    TH1D *hQnQnAEPnonuni[nCoef][DET_N];
    TH1D *hQnAQnBEPnonuni[nCoef][DET_N];
    TH1D *hQnQnASPnonuni[nCoef][DET_N];
    TH1D *hQnAQnBSPnonuni[nCoef][DET_N];
    for (i = 0; i < nCoef; i++) {
        hQnQnAEP[i][iDet] = (TH1D*)fIn->Get(Form("hQnQnAEPH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hQnQnAEP[i][iDet]);
        hQnAQnBEP[i][iDet] = (TH1D*)fIn->Get(Form("hQnAQnBEPH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hQnAQnBEP[i][iDet]);
        hQnQnASP[i][iDet] = (TH1D*)fIn->Get(Form("hQnQnASPH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hQnQnASP[i][iDet]);
        hQnAQnBSP[i][iDet] = (TH1D*)fIn->Get(Form("hQnAQnBSPH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hQnAQnBSP[i][iDet]);

        hQnQnAEPnonuni[i][iDet] = (TH1D*)fIn->Get(Form("hQnQnAEPnonuniH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hQnQnAEPnonuni[i][iDet]);
        hQnAQnBEPnonuni[i][iDet] = (TH1D*)fIn->Get(Form("hQnAQnBEPnonuniH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hQnAQnBEPnonuni[i][iDet]);
        hQnQnASPnonuni[i][iDet] = (TH1D*)fIn->Get(Form("hQnQnASPnonuniH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hQnQnASPnonuni[i][iDet]);
        hQnAQnBSPnonuni[i][iDet] = (TH1D*)fIn->Get(Form("hQnAQnBSPnonuniH%02iD%02i", i+1,iDet));
        checkUnderOverFlow(hQnAQnBSPnonuni[i][iDet]);
    }

    double w = hSqrtSumWeights[iDet]->GetMean();
    double wError = hSqrtSumWeights[iDet]->GetMeanError();
    double wNonuni = hSqrtSumWeightsNonuni[iDet]->GetMean();
    double wNonuniError = hSqrtSumWeightsNonuni[iDet]->GetMeanError();

    double vnEP[nCoef], errorVnEP[nCoef] = {0};
    double vnSP[nCoef], errorVnSP[nCoef] = {0};
    double vnEPnonuni[nCoef], errorVnEPnonuni[nCoef] = {0};
    double vnSPnonuni[nCoef], errorVnSPnonuni[nCoef] = {0};
    for (i=0; i<nCoef; i++) {
        vnEP[i] = CalculateVn(hQnQnAEP[i][iDet]->GetMean(), hQnAQnBEP[i][iDet]->GetMean(), w);
        errorVnEP[i] = CalculateVnError(hQnQnAEP[i][iDet]->GetMean(), hQnAQnBEP[i][iDet]->GetMean(), hQnQnAEP[i][iDet]->GetMeanError(), hQnAQnBEP[i][iDet]->GetMeanError(),  w, wError);
        vnSP[i] = CalculateVn(hQnQnASP[i][iDet]->GetMean(), hQnAQnBSP[i][iDet]->GetMean(), w);
        errorVnSP[i] = CalculateVnError(hQnQnASP[i][iDet]->GetMean(), hQnAQnBSP[i][iDet]->GetMean(), hQnQnASP[i][iDet]->GetMeanError(), hQnAQnBSP[i][iDet]->GetMeanError(),  w, wError);

        vnEPnonuni[i] = CalculateVn(hQnQnAEPnonuni[i][iDet]->GetMean(), hQnAQnBEPnonuni[i][iDet]->GetMean(), wNonuni);
        errorVnEPnonuni[i] = CalculateVnError(hQnQnAEPnonuni[i][iDet]->GetMean(), hQnAQnBEPnonuni[i][iDet]->GetMean(), hQnQnAEPnonuni[i][iDet]->GetMeanError(), hQnAQnBEPnonuni[i][iDet]->GetMeanError(), wNonuni, wNonuniError);
        vnSPnonuni[i] = CalculateVn(hQnQnASPnonuni[i][iDet]->GetMean(), hQnAQnBSPnonuni[i][iDet]->GetMean(), wNonuni);
        errorVnSPnonuni[i] = CalculateVnError(hQnQnASPnonuni[i][iDet]->GetMean(), hQnAQnBSPnonuni[i][iDet]->GetMean(), hQnQnASPnonuni[i][iDet]->GetMeanError(), hQnAQnBSPnonuni[i][iDet]->GetMeanError(), wNonuni, wNonuniError);
    }

    // pT bins
    TH1D *hQnQnAPtBin[nPtBins];
    TH1D *hSqrtSumWeightsPtBins[nPtBins];
    for (i=0; i<nPtBins; i++) {
        hQnQnAPtBin[i] = (TH1D*)fIn->Get(Form("hQnQnAPtBin%02i", i+1));
        checkUnderOverFlow(hQnQnAPtBin[i]);
        hSqrtSumWeightsPtBins[i] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsPtBinsH%02i", i+1));
        checkUnderOverFlow(hSqrtSumWeightsPtBins[i]);
    }

    double ptBin[nPtBins];
    double ptBinError[nPtBins];
    for (i=0; i<nPtBins; i++) {
        ptBin[i] = CalculateVn(hQnQnAPtBin[i]->GetMean(), hQnAQnBEP[1][iDet]->GetMean(), hSqrtSumWeightsPtBins[i]->GetMean());
        ptBinError[i] = CalculateVnError(hQnQnAPtBin[i]->GetMean(), hQnAQnBEP[1][iDet]->GetMean(), hQnQnAPtBin[i]->GetMeanError(), hQnAQnBEP[1][iDet]->GetMeanError(), hSqrtSumWeightsPtBins[i]->GetMean(), hSqrtSumWeightsPtBins[i]->GetMeanError());
    }

    // Make graphs
    TGraphErrors *gVnTrue = new TGraphErrors(nCoef);
    TGraphErrors *gVn = new TGraphErrors(nCoef);
    TGraphErrors *gVnEP = new TGraphErrors(nCoef);
    TGraphErrors *gVnSP = new TGraphErrors(nCoef);
    TGraphErrors *gRtrue = new TGraphErrors(nCoef);
    TGraphErrors *gR = new TGraphErrors(nCoef);

    TGraphErrors *gVnTrueNonuni = new TGraphErrors(nCoef);
    TGraphErrors *gVnNonuni = new TGraphErrors(nCoef);
    TGraphErrors *gVnEPnonuni = new TGraphErrors(nCoef);
    TGraphErrors *gVnSPnonuni = new TGraphErrors(nCoef);
    TGraphErrors *gRtrueNonuni = new TGraphErrors(nCoef);
    TGraphErrors *gRnonuni = new TGraphErrors(nCoef);

    TGraphErrors *gPtBin = new TGraphErrors(6);

    for(i = 0; i < nCoef; i++){
        gVnTrue->SetPoint(i, double(i+1)-0.3, vnTrue[i]);
        gVnTrue->SetPointError(i, 0.0, errorVnTrue[i]);
        gVn->SetPoint(i, double(i+1)-0.1, vnEP[i]);
        gVn->SetPointError(i, 0.0, errorVnEP[i]);
        gVnEP->SetPoint(i, double(i+1)+0.1, vnEP[i]);
        gVnEP->SetPointError(i, 0.0, errorVnEP[i]);
        gVnSP->SetPoint(i, double(i+1)+0.3, vnSP[i]);
        gVnSP->SetPointError(i, 0.0, errorVnSP[i]);
        gRtrue->SetPoint(i, double(i+1)-0.15, hRtrue[i][iDet]->GetMean());
        gRtrue->SetPointError(i, 0.0, hRtrue[i][iDet]->GetMeanError());
        gR->SetPoint(i, double(i+1)+0.05, R[i]);
        gR->SetPointError(i, 0.0, errorR[i]);

        gVnTrueNonuni->SetPoint(i, double(i+1)-0.3, vnTrueNonuni[i]);
        gVnTrueNonuni->SetPointError(i, 0.0, errorVnTrueNonuni[i]);
        gVnNonuni->SetPoint(i, double(i+1)-0.1, vnNonuni[i]);
        gVnNonuni->SetPointError(i, 0.0, errorVnNonuni[i]);
        gVnEPnonuni->SetPoint(i, double(i+1)+0.1, vnEPnonuni[i]);
        gVnEPnonuni->SetPointError(i, 0.0, errorVnEPnonuni[i]);
        gVnSPnonuni->SetPoint(i, double(i+1)+0.3, vnSPnonuni[i]);
        gVnSPnonuni->SetPointError(i, 0.0, errorVnSPnonuni[i]);
        gRtrueNonuni->SetPoint(i, double(i+1)-0.05, hRtrueNonuni[i][iDet]->GetMean());
        gRtrueNonuni->SetPointError(i, 0.0, hRtrueNonuni[i][iDet]->GetMeanError());
        gRnonuni->SetPoint(i, double(i+1)+0.15, Rnonuni[i]);
        gRnonuni->SetPointError(i, 0.0, errorRnonuni[i]);
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
    hInputFlow->Write("hInputFlow");
    gVn->Write("gVn");
    gVnTrue->Write("gVnTrue");
    gR->Write("gR");
    gPtBin->Write("gPtBin");
    gRtrue->Write("gRtrue");
    gVnEP->Write("gVnEP");
    gVnSP->Write("gVnSP");

    hPhiNonuni->Write("hPhiNonuni");
    gVnNonuni->Write("gVnNonuni");
    gVnTrueNonuni->Write("gVnTrueNonuni");
    gRnonuni->Write("gRnonuni");
    gRtrueNonuni->Write("gRtrueNonuni");
    gVnEPnonuni->Write("gVnEPnonuni");
    gVnSPnonuni->Write("gVnSPnonuni");

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
    return vn*TMath::Sqrt(QnQnAerr*QnQnAerr/QnQnA/QnQnA + 0.5*0.5*QnAQnBerr*QnAQnBerr/QnAQnB/QnAQnB + wErr*wErr/w/w);
}

void checkUnderOverFlow( TH1 *h ){
        if(h->GetBinContent(0)>0) cout << h->GetName() << " underflow bin not empty: " << h->GetBinContent(0) << endl;
        if(h->GetBinContent(h->GetXaxis()->GetNbins()+1)>0) cout << h->GetName() << " overflow bin not empty: " << h->GetBinContent(h->GetXaxis()->GetNbins()+1) << endl;
}
