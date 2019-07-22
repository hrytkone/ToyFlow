#include <TMath.h>
#include <TFile.h>

#include "src/ResIter.h"
#include "src/JConst.h"

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr);
void checkUnderOverFlow( TH1 *h );

void MakeCentralityGraphs(TString sInputName = "toyFlow.root", TString sOutputName = "toyFlowCentralityGraphs.root", const int iDet=3) {

    TFile *fIn = TFile::Open(sInputName, "read");
    TFile *fOut = TFile::Open(sOutputName, "recreate");

    int nPtBins = 9;

    int i, j, n;
    double pi = TMath::Pi();

    double inputFlow[nCoef][CENTBINS_N-4] = {
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0277401, 0.04488324, 0.06521883, 0.08433443, 0.09597485, 0.10087206, 0.09925828},
        {0.02039728, 0.02369955, 0.02670301, 0.02950095, 0.03118808, 0.03120636, 0.02918556},
        {0.01013229, 0.01171893, 0.0131265, 0.01479335, 0.0159713, 0.01644628, 0.01535014},
        {0.00415816, 0.00467961, 0.00528238, 0.006501, 0.0068885, 0.00690379, 0.00575251}
    };

    TH1D *hInputFlow[nCoef];
    for (i=0; i<nCoef; i++) {
        hInputFlow[i] = new TH1D("hInputFlow", "hInputFlow", CENTBINS_N-4, centBins);
        hInputFlow[i]->SetLineStyle(1);
        hInputFlow[i]->SetLineColor(1);
        hInputFlow[i]->SetLineWidth(1);
    }
    for(i=0; i<nCoef; i++) {
        for (j=0; j<CENTBINS_N-4; j++) {
            hInputFlow[i]->Fill(centBins[j], inputFlow[i][j]);
        }
    }

    //=====vn=====
    // Observed vn
    TH1D *hVnObs[nCoef][DET_N][CENTBINS_N];
    TH1D *hVnObsNonuni[nCoef][DET_N][CENTBINS_N];
    for (i=0; i<nCoef; i++) {
        for (j=0; j<CENTBINS_N-4; j++) {
            hVnObs[i][iDet][j] = (TH1D*)fIn->Get(Form("hVnObsH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hVnObs[i][iDet][j]);
            hVnObsNonuni[i][iDet][j] = (TH1D*)fIn->Get(Form("hVnObsNonuniH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hVnObsNonuni[i][iDet][j]);
        }
    }

    // Resolution parameter
    TH1D *hRsub[nCoef][DET_N][CENTBINS_N];
    TH1D *hRtrue[nCoef][DET_N][CENTBINS_N];
    TH1D *hRsubNonuni[nCoef][DET_N][CENTBINS_N];
    TH1D *hRtrueNonuni[nCoef][DET_N][CENTBINS_N];
    for (i=0; i<nCoef; i++) {
        for (j=0; j<CENTBINS_N-4; j++) {
            hRsub[i][iDet][j] = (TH1D*)fIn->Get(Form("hRsubH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hRsub[i][iDet][j]);
            hRtrue[i][iDet][j] = (TH1D*)fIn->Get(Form("hRtrueH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hRtrue[i][iDet][j]);
            hRsubNonuni[i][iDet][j] = (TH1D*)fIn->Get(Form("hRsubNonuniH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hRsubNonuni[i][iDet][j]);
            hRtrueNonuni[i][iDet][j] = (TH1D*)fIn->Get(Form("hRtrueNonuniH%02iD%02iCENT%02i",i+1,iDet,j));
            checkUnderOverFlow(hRtrueNonuni[i][iDet][j]);
        }
    }

    double vn[nCoef][CENTBINS_N], errorVn[nCoef][CENTBINS_N] = {0};
    double vnTrue[nCoef][CENTBINS_N], errorVnTrue[nCoef][CENTBINS_N] = {0};
    double vnNonuni[nCoef][CENTBINS_N], errorVnNonuni[nCoef][CENTBINS_N] = {0};
    double vnTrueNonuni[nCoef][CENTBINS_N], errorVnTrueNonuni[nCoef][CENTBINS_N] = {0};
    double Rinit;

    double khi0 = 0.5;
    double err = 0.0001;
    double khi;
    double R[nCoef][CENTBINS_N], errorR[nCoef][CENTBINS_N] = {0};
    double Rnonuni[nCoef][CENTBINS_N], errorRnonuni[nCoef][CENTBINS_N] = {0};

    cout << "UNIFORM-----------\n";
    for (i=0; i<nCoef; i++) {
        cout << "n=" << i+1 << "----------------\n";
        for (j=0; j<CENTBINS_N-4; j++) {
            cout << "cent=" << j << "\n";
            n = i+1;
            vn[i][j] = hVnObs[i][iDet][j]->GetMean();
            vnTrue[i][j] = hVnObs[i][iDet][j]->GetMean();

            Rinit = TMath::Sqrt(hRsub[i][iDet][j]->GetMean());
            khi = RIter(khi0, Rinit, err);
            R[i][j] = R1(TMath::Sqrt(2)*khi); //Because khi sim sqrt(Multi) and full event has twice multiplicity compared to A or B.
            errorR[i][j] = CalculateRerror(khi, err);

            cout << "R=    " << R[i][j] <<                       "  err=" << errorR[i][j] << "\n";
            cout << "Rtrue=" << hRtrue[i][iDet][j]->GetMean() << "  err=" << hRtrue[i][iDet][j]->GetMeanError() << "\n";
            cout << "R/Rtrue=" << R[i][j]/hRtrue[i][iDet][j]->GetMean() << "\n\n";

            vn[i][j] /= R[i][j];
            errorVn[i][j] = GetVnError(hVnObs[i][iDet][j]->GetMean(), hVnObs[i][iDet][j]->GetMeanError(), R[i][j], errorR[i][j]);

            vnTrue[i][j] /= hRtrue[i][iDet][j]->GetMean();
            errorVnTrue[i][j] = GetVnError(hVnObs[i][iDet][j]->GetMean(), hVnObs[i][iDet][j]->GetMeanError(), hRtrue[i][iDet][j]->GetMean(), hRtrue[i][iDet][j]->GetMeanError());
        }
    }

    cout << "NONUNIFORM--------\n";
    for (i=0; i<nCoef; i++) {
        cout << "n=" << i+1 << "----------------\n";
        for (j=0; j<CENTBINS_N-4; j++) {
            cout << "cent=" << j << "\n";
            n = i+1;
            vnNonuni[i][j] = hVnObsNonuni[i][iDet][j]->GetMean();
            vnTrueNonuni[i][j] = hVnObsNonuni[i][iDet][j]->GetMean();

            Rinit = TMath::Sqrt(hRsubNonuni[i][iDet][j]->GetMean());
            khi = RIter(khi0, Rinit, err);
            Rnonuni[i][j] = R1(TMath::Sqrt(2)*khi); //Because khi sim sqrt(Multi) and full event has twice multiplicity compared to A or B.
            errorRnonuni[i][j] = CalculateRerror(khi, err);

            cout << "R=    " << Rnonuni[i][j] <<                       "  err=" << errorRnonuni[i][j] << "\n";
            cout << "Rtrue=" << hRtrueNonuni[i][iDet][j]->GetMean() << "  err=" << hRtrueNonuni[i][iDet][j]->GetMeanError() << "\n";
            cout << "R/Rtrue=" << Rnonuni[i][j]/hRtrueNonuni[i][iDet][j]->GetMean() << "\n\n";

            vnNonuni[i][j] /= Rnonuni[i][j];
            errorVnNonuni[i][j] = GetVnError(hVnObsNonuni[i][iDet][j]->GetMean(), hVnObsNonuni[i][iDet][j]->GetMeanError(), Rnonuni[i][j], errorRnonuni[i][j]);

            vnTrueNonuni[i][j] /= hRtrueNonuni[i][iDet][j]->GetMean();
            errorVnTrueNonuni[i][j] = GetVnError(hVnObsNonuni[i][iDet][j]->GetMean(), hVnObsNonuni[i][iDet][j]->GetMeanError(), hRtrueNonuni[i][iDet][j]->GetMean(), hRtrueNonuni[i][iDet][j]->GetMeanError());
        }
    }

    // Make graphs
    TGraphErrors *gVnTrue[nCoef];
    TGraphErrors *gVn[nCoef];
    TGraphErrors *gRtrue[nCoef];
    TGraphErrors *gR[nCoef];

    for(i=0; i<nCoef; i++){
        gVnTrue[i] = new TGraphErrors(CENTBINS_N-4);
        gVn[i] = new TGraphErrors(CENTBINS_N-4);
        gRtrue[i] = new TGraphErrors(CENTBINS_N-4);
        gR[i] = new TGraphErrors(CENTBINS_N-4);
        for (j=0; j<CENTBINS_N-4; j++) {
            gVnTrue[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2-1.5, vnTrue[i][j]);
            gVnTrue[i]->SetPointError(j, 0.0, errorVnTrue[i][j]);
            gVn[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2-0.5, vn[i][j]);
            gVn[i]->SetPointError(j, 0.0, errorVn[i][j]);
            gRtrue[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2-1.5, hRtrue[i][iDet][j]->GetMean());
            gRtrue[i]->SetPointError(j, 0.0, hRtrue[i][iDet][j]->GetMeanError());
            gR[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2-0.5, R[i][j]);
            gR[i]->SetPointError(j, 0.0, errorR[i][j]);
        }
    }

    TGraphErrors *gVnTrueNonuni[nCoef];
    TGraphErrors *gVnNonuni[nCoef];
    TGraphErrors *gRtrueNonuni[nCoef];
    TGraphErrors *gRnonuni[nCoef];

    for(i=0; i<nCoef; i++){
        gVnTrueNonuni[i] = new TGraphErrors(CENTBINS_N-4);
        gVnNonuni[i] = new TGraphErrors(CENTBINS_N-4);
        gRtrueNonuni[i] = new TGraphErrors(CENTBINS_N-4);
        gRnonuni[i] = new TGraphErrors(CENTBINS_N-4);
        for (j=0; j<CENTBINS_N-4; j++) {
            gVnTrueNonuni[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2+0.5, vnTrueNonuni[i][j]);
            gVnTrueNonuni[i]->SetPointError(j, 0.0, errorVnTrueNonuni[i][j]);
            gVnNonuni[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2+1.5, vnNonuni[i][j]);
            gVnNonuni[i]->SetPointError(j, 0.0, errorVnNonuni[i][j]);
            gRtrueNonuni[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2+0.5, hRtrueNonuni[i][iDet][j]->GetMean());
            gRtrueNonuni[i]->SetPointError(j, 0.0, hRtrueNonuni[i][iDet][j]->GetMeanError());
            gRnonuni[i]->SetPoint(j, (centBins[j]+centBins[j+1])/2+1.5, Rnonuni[i][j]);
            gRnonuni[i]->SetPointError(j, 0.0, errorRnonuni[i][j]);
        }
    }

    fOut->cd();
    for (i=0; i<nCoef; i++) {
        hInputFlow[i]->Write(Form("hInputFlow%02i", i+1));

        gVnTrue[i]->Write(Form("gVnTrueH%02i", i+1));
        gVn[i]->Write(Form("gVnH%02i", i+1));
        gRtrue[i]->Write(Form("gRtrueH%02i", i+1));
        gR[i]->Write(Form("gRH%02i", i+1));

        gVnTrueNonuni[i]->Write(Form("gVnTrueNonuniH%02i", i+1));
        gVnNonuni[i]->Write(Form("gVnNonuniH%02i", i+1));
        gRtrueNonuni[i]->Write(Form("gRtrueNonuniH%02i", i+1));
        gRnonuni[i]->Write(Form("gRnonuniH%02i", i+1));
    }

    fOut->Close();
}

double GetVnError(double vnObs, double vnObsErr, double Rn, double RnErr) {
    return (vnObs/Rn)*TMath::Sqrt((vnObsErr/vnObs)*(vnObsErr/vnObs) + (RnErr/Rn)*(RnErr/Rn));
}

void checkUnderOverFlow( TH1 *h ){
        if(h->GetBinContent(0)>0) cout << h->GetName() << " underflow bin not empty: " << h->GetBinContent(0) << endl;
        if(h->GetBinContent(h->GetXaxis()->GetNbins()+1)>0) cout << h->GetName() << " overflow bin not empty: " << h->GetBinContent(h->GetXaxis()->GetNbins()+1) << endl;
}
