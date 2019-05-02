#include "JHistos.h"
#include "TMath.h"

JHistos::JHistos(){

    int i;
    double pi = TMath::Pi();

    // pT -histograms, taken from JCORRAN version on 11th February 2013
    int NBINS = 150;
    double LogBinsX[NBINS+1], LimL=0.1, LimH=100;
    double logBW = (TMath::Log(LimH)-TMath::Log(LimL))/NBINS;
    for (i=0; i<=NBINS; i++)
        LogBinsX[i] = LimL*exp(i*logBW);

	hPt = new TH1D("hPt","pT - inclusive", NBINS, LogBinsX);
    hPt->Sumw2();
    hMultiplicity = new TH1D("hMultiplicity","Multiplicity - uniform", 125, 0.0, 2501.);
    hMultiplicity->Sumw2();
    hPhi = new TH1D("hPhi", "phi - uniform", 129, -pi, pi);
    hPhi->Sumw2();
    hCentrality = new TH1D("hCentrality", "centrality", 8, 0.0, 70.0);
    hCentrality->Sumw2();

    hMultiplicityNonuni = new TH1D("hMultiplicityNonuni","Multiplicity - nonuniform", 125, 0.0, 2501.);
    hMultiplicityNonuni->Sumw2();
    hPhiNonuni = new TH1D("hPhiNonuni", "phi - nonuniform", 129, -pi, pi);
    hPhiNonuni->Sumw2();

    hSqrtSumWeightsTPC = new TH1D("hSqrtSumWeightsTPC","sqrt of sum of weights squares - TPC", 240, 0.0, 60.0);
    hSqrtSumWeightsTPC->Sumw2();
    hSqrtSumWeightsTPCA = new TH1D("hSqrtSumWeightsTPCA","sqrt of sum of weights squares - TPC-A", 240, 0.0, 60.0);
    hSqrtSumWeightsTPCA->Sumw2();
    hSqrtSumWeightsTPCC = new TH1D("hSqrtSumWeightsTPCC","sqrt of sum of weights squares - TPC-C", 240, 0.0, 60.0);
    hSqrtSumWeightsTPCC->Sumw2();
    hSqrtSumWeightsV0P = new TH1D("hSqrtSumWeightsV0P","sqrt of sum of weights squares - V0+", 240, 0.0, 60.0);
    hSqrtSumWeightsV0P->Sumw2();

    hSqrtSumWeightsTPCNonuni = new TH1D("hSqrtSumWeightsTPCNonuni","sqrt of sum of weights squares - nonuni TPC", 240, 0.0, 60.0);
    hSqrtSumWeightsTPCNonuni->Sumw2();
    hSqrtSumWeightsTPCANonuni = new TH1D("hSqrtSumWeightsTPCANonuni","sqrt of sum of weights squares - nonuni TPC-A", 240, 0.0, 60.0);
    hSqrtSumWeightsTPCANonuni->Sumw2();
    hSqrtSumWeightsTPCCNonuni = new TH1D("hSqrtSumWeightsTPCCNonuni","sqrt of sum of weights squares - nonuni TPC-C", 240, 0.0, 60.0);
    hSqrtSumWeightsTPCCNonuni->Sumw2();
    hSqrtSumWeightsV0PNonuni = new TH1D("hSqrtSumWeightsV0PNonuni","sqrt of sum of weights squares - nonuni V0P", 240, 0.0, 60.0);
    hSqrtSumWeightsV0PNonuni->Sumw2();

    for (i=0; i<5; i++){

        hRtrue[i] = new TH1D(Form("hRtrue%02i",i+1),Form("hRtrue%02i",i+1),401,-1.0,1.0);
        hRtrue[i]->Sumw2();
        hRsub[i] = new TH1D(Form("hRsub%02i",i+1),Form("hRsub%02i",i+1),401,-1.0,1.0);
        hRsub[i]->Sumw2();
        hVnObs[i] = new TH1D(Form("hVnObs%02i",i+1),Form("hVnObs%02i",i+1),401,-1.5,1.5);
        hVnObs[i]->Sumw2();

        hRtrueCorrected[i] = new TH1D(Form("hRtrueCorrected%02i",i+1),Form("hRtrueCorrected%02i",i+1),401,-1.0,1.0);
        hRtrueCorrected[i]->Sumw2();
        hRsubCorrected[i] = new TH1D(Form("hRsubCorrected%02i",i+1),Form("hRsubCorrected%02i",i+1),401,-1.0,1.0);
        hRsubCorrected[i]->Sumw2();
        hVnObsCorrected[i] = new TH1D(Form("hVnObsCorrected%02i",i+1),Form("hVnObsCorrected%02i",i+1),401,-1.5,1.5);
        hVnObsCorrected[i]->Sumw2();

        hQnQnA[i] = new TH1D(Form("hQnQnA%02i",i+1),Form("hQnQnA%02i",i+1),401,-50.0,50.0);
        hQnQnA[i]->Sumw2();
        hQnAQnB[i] = new TH1D(Form("hQnAQnB%02i",i+1),Form("hQnAQnB%02i",i+1),401,-1.0,1.0);
        hQnAQnB[i]->Sumw2();

        hQnQnAcorrected[i] = new TH1D(Form("hQnQnAcorrected%02i",i+1),Form("hQnQnAcorrected%02i",i+1),401,-50.0,50.0);
        hQnQnAcorrected[i]->Sumw2();
        hQnAQnBcorrected[i] = new TH1D(Form("hQnAQnBcorrected%02i",i+1),Form("hQnAQnBcorrected%02i",i+1),401,-1.0,1.0);
        hQnAQnBcorrected[i]->Sumw2();

    }

    for (i=0; i<9; i++) {
        hPtBin[i] = new TH1D(Form("hPtBin%02i", i+1), Form("hPtBin%02i", i+1), 401, -10.0, 10.0);
        hPtBin[i]->Sumw2();
        hSqrtSumWeightsPtBins[i] = new TH1D(Form("hSqrtSumWeightsPtBins%02i", i+1),Form("sqrt of sum of weights squares for pT bin %02i", i+1), 240, 0.0, 60.0);
        hSqrtSumWeightsPtBins[i]->Sumw2();
    }

    //FOR TESTING
    hV2ComplexPart = new TH1D(Form("hV2ComplexPart%02i",i+1), Form("hV2ComplexPart%02i",i+1), 401, -5.0, 5.0);
}
