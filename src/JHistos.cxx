#include "JHistos.h"
#include "TMath.h"

JHistos::JHistos(){

    int i, j, k;
    double pi = TMath::Pi();

    int NBINS = 150;
    double LogBinsX[NBINS+1], LimL=0.1, LimH=100;
    double logBW = (TMath::Log(LimH)-TMath::Log(LimL))/NBINS;
    for (i=0; i<=NBINS; i++) {
        LogBinsX[i] = LimL*exp(i*logBW);
    }

	hPt = new TH1D("hPt","pT - inclusive", NBINS, LogBinsX);
    hPt->Sumw2();
    hPhi = new TH1D("hPhi", "phi - uniform", 129, -pi, pi);
    hPhi->Sumw2();
    hCentrality = new TH1D("hCentrality", "centrality", CENTBINS_N, 0.0, 90.0);
    hCentrality->Sumw2();
    hEta = new TH1D("hEta", "pseudorapidity", 401, -3.8, 5.2);


    for (j=0; j<CENTBINS_N; j++){
        hMultiplicity[j] = new TH1D(Form("hMultiplicityCENT%02i",j), "Total multiplicity", 300, 0.0, 30000.);
        hMultiplicity[j]->Sumw2();
        for (i=0; i<DET_N; i++) {
            hMultiPerDet[i][j] = new TH1D(Form("hMultiPerDetD%02iCENT%02i",i,j), "Detector multiplicity", 300, 0.0, 30000.);
            hMultiPerDet[i][j]->Sumw2();
            hSqrtSumWeights[i][j]= new TH1D(Form("hSqrtSumWeightsD%02iCENT%02i",i,j),"sqrt of sum of weights squares", 480, 0.0, 120.0);
            hSqrtSumWeights[i][j]->Sumw2();
        }
    }

    for (i=0; i<nCoef; i++){
        for (j=0; j<DET_N; j++){
            for (k=0; k<CENTBINS_N; k++){
                hRtrue[i][j][k] = new TH1D(Form("hRtrueH%02iD%02iCENT%02i",i+1,j,k),Form("hRtrueH%02iD%02iCENT%02i",i+1,j,k),404,-1.01,1.01);
                hRtrue[i][j][k]->Sumw2();
                hRsub[i][j][k] = new TH1D(Form("hRsubH%02iD%02iCENT%02i",i+1,j,k),Form("hRsubH%02iD%02iCENT%02i",i+1,j,k),404,-1.01,1.01);
                hRsub[i][j][k]->Sumw2();
                hVnObs[i][j][k] = new TH1D(Form("hVnObsH%02iD%02iCENT%02i",i+1,j,k),Form("hVnObsH%02iD%02iCENT%02i",i+1,j,k),401,-1.5,1.5);
                hVnObs[i][j][k]->Sumw2();

                hQnQnAEP[i][j][k] = new TH1D(Form("hQnQnAEPH%02iD%02iCENT%02i",i+1,j,k),Form("hQnQnAEPH%02iD%02iCENT%02i",i+1,j,k),401,-50.0,50.0);
                hQnQnAEP[i][j][k]->Sumw2();
                hQnAQnBEP[i][j][k] = new TH1D(Form("hQnAQnBEPH%02iD%02iCENT%02i",i+1,j,k),Form("hQnAQnBEPH%02iD%02iCENT%02i",i+1,j,k),404,-1.01,1.01);
                hQnAQnBEP[i][j][k]->Sumw2();

                hQnQnASP[i][j][k] = new TH1D(Form("hQnQnASPH%02iD%02iCENT%02i",i+1,j,k),Form("hQnQnASPH%02iD%02iCENT%02i",i+1,j,k),401,-10.0,110.0);
                hQnQnASP[i][j][k]->Sumw2();
                hQnAQnBSP[i][j][k] = new TH1D(Form("hQnAQnBSPH%02iD%02iCENT%02i",i+1,j,k),Form("hQnAQnBSPH%02iD%02iCENT%02i",i+1,j,k),401,-10.0,110.0);
                hQnAQnBSP[i][j][k]->Sumw2();

                hQvec[i][j][k] = new TH2D(Form("hQvec%02iD%02iCENT%02i",i+1,j,k),Form("hQvec%02iD%02iCENT%02i",i+1,j,k),101,-100.0,100.0, 101, -100.0, 100.0);
                hQvec[i][j][k]->Sumw2();
            }
        }
    }

    for (int i=0; i<nCoef; i++) {
        for (int j=0; j<CENTBINS_N; j++) {
            hRsubAB[i][j] = new TH1D(Form("hRsubAB%dCENT%02d", i+1, j), Form("hRsubAB%dCENT%02d", i+1, j), 100, -1.0, 1.0);
            hRsubAC[i][j] = new TH1D(Form("hRsubAC%dCENT%02d", i+1, j), Form("hRsubAC%dCENT%02d", i+1, j), 100, -1.0, 1.0);
            hRsubBC[i][j] = new TH1D(Form("hRsubBC%dCENT%02d", i+1, j), Form("hRsubBC%dCENT%02d", i+1, j), 100, -1.0, 1.0);
        }
    }

    for (i=0; i<PTBINS_N; i++) {
        hQnQnAPtBin[i] = new TH1D(Form("hQnQnAPtBin%02i", i+1), Form("hQnQnAPtBin%02i", i+1), 482, -11.0, 11.0);
        hQnQnAPtBin[i]->Sumw2();
        hSqrtSumWeightsPtBins[i] = new TH1D(Form("hSqrtSumWeightsPtBinsH%02i", i+1),Form("sqrt of sum of weights squares for pT bin %02i", i+1), 240, 0.0, 60.0);
        hSqrtSumWeightsPtBins[i]->Sumw2();
    }

    //FOR TESTING
    hV2ComplexPart = new TH1D(Form("hV2ComplexPartH%02i",i+1), Form("hV2ComplexPartH%02i",i+1), 401, -5.0, 5.0);
}
