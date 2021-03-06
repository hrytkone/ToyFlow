#include <TMath.h>
#include <TFile.h>

#include "src/ResIter.h"
#include "src/JConst.h"
#include "src/Filipad.h"

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr);
double CalculateVn(double QnQnA, double QnAQnB, double w);
double CalculateVnError(double QnQnA, double QnAQnB, double QnQnAerr, double QnAQnBerr, double w, double wErr);
void checkUnderOverFlow( TH1 *h );
void ErrorExit(TString error, int errnum=1 ){cout<<"ERROR: "<<error<<endl;gSystem->Exit(errnum);}

void MakeGraphs(TString sInputName = "toyFlow.root", TString sOutputName = "toyFlowGraphs.root", const int iDet=0, const int centBin=0) {

    TFile *fIn = TFile::Open(sInputName, "read");
    TFile *fOut = TFile::Open(sOutputName, "recreate");

    int nPtBins = 9;

    int i, n;
    double pi = TMath::Pi();
    double inputFlow[nCoef] = {0.0};

    TH1D *hInputNumbers = (TH1D*)fIn->Get("toyflow/hInputNumbers/hInputNumbers");
    checkUnderOverFlow(hInputNumbers);
    double nEvents = hInputNumbers->GetBinContent(1);
    double nOfFiles = hInputNumbers->GetBinContent(14);
    for(int i = 0; i < nCoef; i++) {
        inputFlow[i] = hInputNumbers->GetBinContent(2+i);
        inputFlow[i] /= nOfFiles;
    }

    TH1D *hSqrtSumWeights[DET_N][CENTBINS_N];
    hSqrtSumWeights[iDet][centBin] = (TH1D*)fIn->Get(Form("toyflow/hSqrtSumWeights/hSqrtSumWeightsD%02iCENT%02i",iDet,centBin));
    checkUnderOverFlow(hSqrtSumWeights[iDet][centBin]);

    TH1D *hInputFlow = new TH1D("hInputFlow", "hInputFlow", nCoef, 0.5, double(nCoef)+0.5);
    hInputFlow->SetLineStyle(1);
    hInputFlow->SetLineColor(1);
    hInputFlow->SetLineWidth(1);
    for(int i = 0; i < nCoef; i++)
        hInputFlow->Fill(double(i+1), inputFlow[i]);
    hInputFlow->Fill(double(nCoef+1), 0.0);

    TH1D *hPhi = (TH1D*)fIn->Get("toyflow/hPhi/hPhi");
    checkUnderOverFlow(hPhi);

    //=====vn=====
    // Observed vn
    TH1D *hVnObs[nCoef][DET_N][CENTBINS_N];
    for (i = 0; i < nCoef; i++) {
        hVnObs[i][iDet][centBin] = (TH1D*)fIn->Get(Form("toyflow/hVnObs/hVnObsH%02iD%02iCENT%02i",i+1,iDet,centBin));
        checkUnderOverFlow(hVnObs[i][iDet][centBin]);
    }

    // Resolution parameter
    TH1D *hRsub[nCoef][DET_N][CENTBINS_N];
    TH1D *hRtrue[nCoef][DET_N][CENTBINS_N];
    for (i = 0; i < nCoef; i++) {
        hRsub[i][iDet][centBin] = (TH1D*)fIn->Get(Form("toyflow/hRsub/hRsubH%02iD%02iCENT%02i",i+1,iDet,centBin));
        checkUnderOverFlow(hRsub[i][iDet][centBin]);
        hRtrue[i][iDet][centBin] = (TH1D*)fIn->Get(Form("toyflow/hRtrue/hRtrueH%02iD%02iCENT%02i",i+1,iDet,centBin));
        checkUnderOverFlow(hRtrue[i][iDet][centBin]);
    }

    double vn[nCoef], errorVn[nCoef] = {0};
    double vnTrue[nCoef], errorVnTrue[nCoef] = {0};

    double vnRatio[nCoef], errorVnRatio[nCoef] = {0};
    double vnTrueRatio[nCoef], errorVnTrueRatio[nCoef] = {0};

    double Rinit;

    double khi0 = 0.5;
    double err = 0.0001;
    double khi;
    double R[nCoef], errorR[nCoef] = {0};
    double Rnonuni[nCoef], errorRnonuni[nCoef] = {0};

    for (i = 0; i < nCoef; i++) {
        n = i+1;
        vn[i] = hVnObs[i][iDet][centBin]->GetMean();
        vnTrue[i] = hVnObs[i][iDet][centBin]->GetMean();

        Rinit = TMath::Sqrt(hRsub[i][iDet][centBin]->GetMean());
        khi = RIter(khi0, Rinit, err);
        R[i] = R1(TMath::Sqrt(2)*khi); //Because khi sim sqrt(Multi) and full event has twice multiplicity compared to A or B.
        errorR[i] = CalculateRerror(khi, err);

        cout << "R=    " << R[i] <<                       "  err=" << errorR[i] << "\n";
        cout << "Rtrue=" << hRtrue[i][iDet][centBin]->GetMean() << "  err=" << hRtrue[i][iDet][centBin]->GetMeanError() << "\n";
        cout << "R/Rtrue=" << R[i]/hRtrue[i][iDet][centBin]->GetMean() << endl;

        vn[i] /= R[i];
        errorVn[i] = GetVnError(hVnObs[i][iDet][centBin]->GetMean(), hVnObs[i][iDet][centBin]->GetMeanError(), R[i], errorR[i]);

        vnTrue[i] /= hRtrue[i][iDet][centBin]->GetMean();
        errorVnTrue[i] = GetVnError(hVnObs[i][iDet][centBin]->GetMean(), hVnObs[i][iDet][centBin]->GetMeanError(), hRtrue[i][iDet][centBin]->GetMean(), hRtrue[i][iDet][centBin]->GetMeanError());

        if(inputFlow[i]==0) {
            vnRatio[i] = 0;
            errorVnRatio[i] = 0;
            vnTrueRatio[i] = 0;
            errorVnTrueRatio[i] = 0;

        } else {
            vnRatio[i] = vn[i]/inputFlow[i];
            errorVnRatio[i] = errorVn[i]/inputFlow[i];
            vnTrueRatio[i] = vnTrue[i]/inputFlow[i];
            errorVnTrueRatio[i] = errorVnTrue[i]/inputFlow[i];
        }
    }

    //vn{EP} and vn{SP}
    TH1D *hQnQnAEP[nCoef][DET_N][CENTBINS_N];
    TH1D *hQnAQnBEP[nCoef][DET_N][CENTBINS_N];
    TH1D *hQnQnASP[nCoef][DET_N][CENTBINS_N];
    TH1D *hQnAQnBSP[nCoef][DET_N][CENTBINS_N];
    for (i = 0; i < nCoef; i++) {
        hQnQnAEP[i][iDet][centBin] = (TH1D*)fIn->Get(Form("toyflow/hQnQnAEP/hQnQnAEPH%02iD%02iCENT%02i",i+1,iDet,centBin));
        checkUnderOverFlow(hQnQnAEP[i][iDet][centBin]);
        hQnAQnBEP[i][iDet][centBin] = (TH1D*)fIn->Get(Form("toyflow/hQnAQnBEP/hQnAQnBEPH%02iD%02iCENT%02i",i+1,iDet,centBin));
        checkUnderOverFlow(hQnAQnBEP[i][iDet][centBin]);
        hQnQnASP[i][iDet][centBin] = (TH1D*)fIn->Get(Form("toyflow/hQnQnASP/hQnQnASPH%02iD%02iCENT%02i",i+1,iDet,centBin));
        checkUnderOverFlow(hQnQnASP[i][iDet][centBin]);
        hQnAQnBSP[i][iDet][centBin] = (TH1D*)fIn->Get(Form("toyflow/hQnAQnBSP/hQnAQnBSPH%02iD%02iCENT%02i",i+1,iDet,centBin));
        checkUnderOverFlow(hQnAQnBSP[i][iDet][centBin]);
    }

    double w = hSqrtSumWeights[iDet][centBin]->GetMean();
    double wError = hSqrtSumWeights[iDet][centBin]->GetMeanError();

    double vnEP[nCoef], errorVnEP[nCoef] = {0};
    double vnSP[nCoef], errorVnSP[nCoef] = {0};

    double vnEPRatio[nCoef], errorVnEPRatio[nCoef] = {0};
    double vnSPRatio[nCoef], errorVnSPRatio[nCoef] = {0};
    for (i=0; i<nCoef; i++) {
        vnEP[i] = CalculateVn(hQnQnAEP[i][iDet][centBin]->GetMean(), hQnAQnBEP[i][iDet][centBin]->GetMean(), w);
        errorVnEP[i] = CalculateVnError(hQnQnAEP[i][iDet][centBin]->GetMean(), hQnAQnBEP[i][iDet][centBin]->GetMean(), hQnQnAEP[i][iDet][centBin]->GetMeanError(), hQnAQnBEP[i][iDet][centBin]->GetMeanError(),  w, wError);
        vnSP[i] = CalculateVn(hQnQnASP[i][iDet][centBin]->GetMean(), hQnAQnBSP[i][iDet][centBin]->GetMean(), w);
        errorVnSP[i] = CalculateVnError(hQnQnASP[i][iDet][centBin]->GetMean(), hQnAQnBSP[i][iDet][centBin]->GetMean(), hQnQnASP[i][iDet][centBin]->GetMeanError(), hQnAQnBSP[i][iDet][centBin]->GetMeanError(),  w, wError);

        if(inputFlow[i]==0) {
            vnEPRatio[i] = 0;
            errorVnEPRatio[i] = 0;
            vnSPRatio[i] = 0;
            errorVnSPRatio[i] = 0;

        } else {
            vnEPRatio[i] = vnEP[i]/inputFlow[i];
            errorVnEPRatio[i] = errorVnEP[i]/inputFlow[i];
            vnSPRatio[i] = vnSP[i]/inputFlow[i];
            errorVnSPRatio[i] = errorVnSP[i]/inputFlow[i];
        }
    }

    // pT bins
    TH1D *hQnQnAPtBinned[nPtBins];
    TH1D *hSqrtSumWeightsPtBinned[nPtBins];
    for (i=0; i<nPtBins; i++) {
        hQnQnAPtBinned[i] = (TH1D*)fIn->Get(Form("toyflow/hQnQnAPtBinned/hQnQnAPtBinnedPtBin%02i",i+1));
        if(hQnQnAPtBinned[i]==0) ErrorExit("Cannot open histo");
        checkUnderOverFlow(hQnQnAPtBinned[i]);
        hSqrtSumWeightsPtBinned[i] = (TH1D*)fIn->Get(Form("toyflow/hSqrtSumWeightsPtBinned/hSqrtSumWeightsPtBinnedPtBin%02i",i+1));
        if(hSqrtSumWeightsPtBinned[i]==0) ErrorExit("Cannot open histo sqrt");
        checkUnderOverFlow(hSqrtSumWeightsPtBinned[i]);
    }

    double ptBin[nPtBins];
    double ptBinError[nPtBins];
    for (i=0; i<nPtBins; i++) {
        ptBin[i] = CalculateVn(hQnQnAPtBinned[i]->GetMean(), hQnAQnBEP[1][iDet][centBin]->GetMean(), hSqrtSumWeightsPtBinned[i]->GetMean());
        ptBinError[i] = CalculateVnError(hQnQnAPtBinned[i]->GetMean(), hQnAQnBEP[1][iDet][centBin]->GetMean(), hQnQnAPtBinned[i]->GetMeanError(), hQnAQnBEP[1][iDet][centBin]->GetMeanError(), hSqrtSumWeightsPtBinned[i]->GetMean(), hSqrtSumWeightsPtBinned[i]->GetMeanError());
    }

    // Make graphs
    TGraphErrors *gVnTrue = new TGraphErrors(nCoef);
    TGraphErrors *gVn = new TGraphErrors(nCoef);
    TGraphErrors *gVnEP = new TGraphErrors(nCoef);
    TGraphErrors *gVnSP = new TGraphErrors(nCoef);
    TGraphErrors *gRtrue = new TGraphErrors(nCoef);
    TGraphErrors *gR = new TGraphErrors(nCoef);

    TGraphErrors *gVnTrueRatio = new TGraphErrors(nCoef);
    TGraphErrors *gVnRatio = new TGraphErrors(nCoef);
    TGraphErrors *gVnEPRatio = new TGraphErrors(nCoef);
    TGraphErrors *gVnSPRatio = new TGraphErrors(nCoef);

    TGraphErrors *gPtBin = new TGraphErrors(6);

    for(i = 0; i < nCoef; i++){
        gVnTrue->SetPoint(i, double(i+1)-0.3, vnTrue[i]);
        gVnTrue->SetPointError(i, 0.0, errorVnTrue[i]);
        gVn->SetPoint(i, double(i+1)-0.1, vn[i]);
        gVn->SetPointError(i, 0.0, errorVn[i]);
        gVnEP->SetPoint(i, double(i+1)+0.1, vnEP[i]);
        gVnEP->SetPointError(i, 0.0, errorVnEP[i]);
        gVnSP->SetPoint(i, double(i+1)+0.3, vnSP[i]);
        gVnSP->SetPointError(i, 0.0, errorVnSP[i]);
        gRtrue->SetPoint(i, double(i+1)-0.15, hRtrue[i][iDet][centBin]->GetMean());
        gRtrue->SetPointError(i, 0.0, hRtrue[i][iDet][centBin]->GetMeanError());
        gR->SetPoint(i, double(i+1)+0.05, R[i]);
        gR->SetPointError(i, 0.0, errorR[i]);

        gVnTrueRatio->SetPoint(i, double(i+1)-0.3, vnTrueRatio[i]);
        gVnTrueRatio->SetPointError(i, 0.0, errorVnTrueRatio[i]);
        gVnRatio->SetPoint(i, double(i+1)-0.1, vnRatio[i]);
        gVnRatio->SetPointError(i, 0.0, errorVnRatio[i]);
        gVnEPRatio->SetPoint(i, double(i+1)+0.1, vnEPRatio[i]);
        gVnEPRatio->SetPointError(i, 0.0, errorVnEPRatio[i]);
        gVnSPRatio->SetPoint(i, double(i+1)+0.3, vnSPRatio[i]);
        gVnSPRatio->SetPointError(i, 0.0, errorVnSPRatio[i]);
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

    gVnRatio->Write("gVnRatio");
    gVnTrueRatio->Write("gVnTrueRatio");
    gVnEPRatio->Write("gVnEPRatio");
    gVnSPRatio->Write("gVnSPRatio");

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
