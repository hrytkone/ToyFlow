#include "JHistos.h"
#include "TMath.h"

JHistos::JHistos(){

    int i, j;
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

    hPhiNonuni = new TH1D("hPhiNonuni", "phi - nonuniform", 129, -pi, pi);
    hPhiNonuni->Sumw2();

    hMultiplicity = new TH1D("hMultiplicity", "Multiplicity - uniform", 300, 0.0, 30000.);
    hMultiplicity->Sumw2();
    hMultiplicityNonuni = new TH1D("hMultiplicityNonuni", "Multiplicity - non-uniform", 300, 0.0, 30000);
    hMultiplicityNonuni->Sumw2();


    for (i=0; i<DET_N; i++){
        hSqrtSumWeights[i]= new TH1D(Form("hSqrtSumWeightsD%02i",i),"sqrt of sum of weights squares", 240, 0.0, 60.0);
        hSqrtSumWeights[i]->Sumw2();
        hSqrtSumWeightsNonuni[i] = new TH1D(Form("hSqrtSumWeightsNonuniD%02i",i),"sqrt of sum of weights squares - nonuni", 240, 0.0, 60.0);
        hSqrtSumWeightsNonuni[i]->Sumw2();
    }

    for (i=0; i<nCoef; i++){
        for (j=0; j<DET_N; j++){
            hRtrue[i][j] = new TH1D(Form("hRtrueH%02iD%02i",i+1,j),Form("hRtrueH%02iD%02i",i+1,j),401,-1.0,1.0);
            hRtrue[i][j]->Sumw2();
            hRsub[i][j] = new TH1D(Form("hRsubH%02iD%02i",i+1,j),Form("hRsubH%02iD%02i",i+1,j),401,-1.0,1.0);
            hRsub[i][j]->Sumw2();
            hVnObs[i][j] = new TH1D(Form("hVnObsH%02iD%02i",i+1,j),Form("hVnObsH%02iD%02i",i+1,j),401,-1.5,1.5);
            hVnObs[i][j]->Sumw2();

            hRtrueCorrected[i][j] = new TH1D(Form("hRtrueCorrectedH%02iD%02i",i+1,j),Form("hRtrueCorrectedH%02iD%02i",i+1,j),401,-1.0,1.0);
            hRtrueCorrected[i][j]->Sumw2();
            hRsubCorrected[i][j] = new TH1D(Form("hRsubCorrectedH%02iD%02i",i+1,j),Form("hRsubCorrectedH%02iD%02i",i+1,j),401,-1.0,1.0);
            hRsubCorrected[i][j]->Sumw2();
            hVnObsCorrected[i][j] = new TH1D(Form("hVnObsCorrectedH%02iD%02i",i+1,j),Form("hVnObsCorrectedH%02iD%02i",i+1,j),401,-1.5,1.5);
            hVnObsCorrected[i][j]->Sumw2();

            hQnQnA[i][j] = new TH1D(Form("hQnQnAH%02iD%02i",i+1,j),Form("hQnQnAH%02iD%02i",i+1,j),401,-50.0,50.0);
            hQnQnA[i][j]->Sumw2();
            hQnAQnB[i][j] = new TH1D(Form("hQnAQnBH%02iD%02i",i+1,j),Form("hQnAQnBH%02iD%02i",i+1,j),401,-1.0,1.0);
            hQnAQnB[i][j]->Sumw2();

            hQnQnAcorrected[i][j] = new TH1D(Form("hQnQnAcorrectedH%02iD%02i",i+1,j),Form("hQnQnAcorrectedH%02iD%02i",i+1,j),401,-50.0,50.0);
            hQnQnAcorrected[i][j]->Sumw2();
            hQnAQnBcorrected[i][j] = new TH1D(Form("hQnAQnBcorrectedH%02iD%02i",i+1,j),Form("hQnAQnBcorrectedH%02iD%02i",i+1,j),401,-1.0,1.0);
            hQnAQnBcorrected[i][j]->Sumw2();

        }
    }

    for (i=0; i<PTBINS_N; i++) {
        hPtBin[i] = new TH1D(Form("hPtBinH%02i", i+1), Form("hPtBin%02i", i+1), 401, -10.0, 10.0);
        hPtBin[i]->Sumw2();
        hSqrtSumWeightsPtBins[i] = new TH1D(Form("hSqrtSumWeightsPtBinsH%02i", i+1),Form("sqrt of sum of weights squares for pT bin %02i", i+1), 240, 0.0, 60.0);
        hSqrtSumWeightsPtBins[i]->Sumw2();
    }

    //FOR TESTING
    hV2ComplexPart = new TH1D(Form("hV2ComplexPartH%02i",i+1), Form("hV2ComplexPartH%02i",i+1), 401, -5.0, 5.0);
}
