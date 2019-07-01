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


    for (i=0; i<DET_N; i++) {
        hSqrtSumWeights[i]= new TH1D(Form("hSqrtSumWeightsD%02i",i),"sqrt of sum of weights squares", 480, 0.0, 120.0);
        hSqrtSumWeights[i]->Sumw2();
        hSqrtSumWeightsNonuni[i] = new TH1D(Form("hSqrtSumWeightsNonuniD%02i",i),"sqrt of sum of weights squares - nonuni", 480, 0.0, 120.0);
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

            hRtrueNonuni[i][j] = new TH1D(Form("hRtrueNonuniH%02iD%02i",i+1,j),Form("hRtrueNonuniH%02iD%02i",i+1,j),401,-1.0,1.0);
            hRtrueNonuni[i][j]->Sumw2();
            hRsubNonuni[i][j] = new TH1D(Form("hRsubNonuniH%02iD%02i",i+1,j),Form("hRsubNonuniH%02iD%02i",i+1,j),401,-1.0,1.0);
            hRsubNonuni[i][j]->Sumw2();
            hVnObsNonuni[i][j] = new TH1D(Form("hVnObsNonuniH%02iD%02i",i+1,j),Form("hVnObsNonuniH%02iD%02i",i+1,j),401,-1.5,1.5);
            hVnObsNonuni[i][j]->Sumw2();

            hQnQnAEP[i][j] = new TH1D(Form("hQnQnAEPH%02iD%02i",i+1,j),Form("hQnQnAEPH%02iD%02i",i+1,j),401,-50.0,50.0);
            hQnQnAEP[i][j]->Sumw2();
            hQnAQnBEP[i][j] = new TH1D(Form("hQnAQnBEPH%02iD%02i",i+1,j),Form("hQnAQnBEPH%02iD%02i",i+1,j),401,-1.0,1.0);
            hQnAQnBEP[i][j]->Sumw2();

            hQnQnAEPnonuni[i][j] = new TH1D(Form("hQnQnAEPnonuniH%02iD%02i",i+1,j),Form("hQnQnAEPnonuniH%02iD%02i",i+1,j),401,-50.0,50.0);
            hQnQnAEPnonuni[i][j]->Sumw2();
            hQnAQnBEPnonuni[i][j] = new TH1D(Form("hQnAQnBEPnonuniH%02iD%02i",i+1,j),Form("hQnAQnBEPnonuniH%02iD%02i",i+1,j),401,-1.0,1.0);
            hQnAQnBEPnonuni[i][j]->Sumw2();

            hQnQnASP[i][j] = new TH1D(Form("hQnQnASPH%02iD%02i",i+1,j),Form("hQnQnASPH%02iD%02i",i+1,j),401,-10.0,110.0);
            hQnQnASP[i][j]->Sumw2();
            hQnAQnBSP[i][j] = new TH1D(Form("hQnAQnBSPH%02iD%02i",i+1,j),Form("hQnAQnBSPH%02iD%02i",i+1,j),401,-10.0,110.0);
            hQnAQnBSP[i][j]->Sumw2();

            hQnQnASPnonuni[i][j] = new TH1D(Form("hQnQnASPnonuniH%02iD%02i",i+1,j),Form("hQnQnASPnonuniH%02iD%02i",i+1,j),401,-10.0,110.0);
            hQnQnASPnonuni[i][j]->Sumw2();
            hQnAQnBSPnonuni[i][j] = new TH1D(Form("hQnAQnBSPnonuniH%02iD%02i",i+1,j),Form("hQnAQnBSPnonuniH%02iD%02i",i+1,j),401,-10.0,110.0);
            hQnAQnBSPnonuni[i][j]->Sumw2();
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
