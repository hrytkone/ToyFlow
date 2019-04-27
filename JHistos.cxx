#include "JHistos.h"
#include "TMath.h"

JHistos::JHistos(){

    int i;

    // pT -histograms, taken from JCORRAN version on 11th February 2013
    int NBINS = 150;
    double LogBinsX[NBINS+1], LimL=0.1, LimH=100;
    double logBW = (TMath::Log(LimH)-TMath::Log(LimL))/NBINS;
    for (i=0; i<=NBINS; i++)
        LogBinsX[i] = LimL*exp(i*logBW);

    // xT -histograms, use also here logarithmic bins
    int binsQ = 400;
    double LogQ2BinsX[binsQ+1];
    double LimQL = 0.1, LimQ2H=2.e2;
    double logQ2BW = (TMath::Log(LimQ2H)-TMath::Log(LimQL))/binsQ;
    for(i=0; i<=binsQ; i++)
        LogQ2BinsX[i] = LimQL*exp(i*logQ2BW);

	hPt = new TH1D("hPt","pT - inclusive", NBINS, LogBinsX);
    hPt->Sumw2();
    hMultiplicity = new TH1D("hMultiplicity","Multiplicity - uniform", 125, 0.0, 2501.);
    hMultiplicity->Sumw2();
    hPhi = new TH1D("hPhi", "phi - uniform", 129, -3.2, 3.2);
    hPhi->Sumw2();

    hMultiplicityNonuni = new TH1D("hMultiplicityNonuni","Multiplicity - nonuniform", 125, 0.0, 2501.);
    hMultiplicityNonuni->Sumw2();
    hPhiNonuni = new TH1D("hPhiNonuni", "phi - nonuniform", 129, -3.2, 3.2);
    hPhiNonuni->Sumw2();

    hSqrtSumWeights = new TH1D("hSqrtSumWeights","sqrt of sum of weights squares", 240, 0.0, 60.0);
    hSqrtSumWeights->Sumw2();
    hSqrtSumWeightsA = new TH1D("hSqrtSumWeightsA","sqrt of sum of weights squares in subevent A", 240, 0.0, 60.0);
    hSqrtSumWeightsA->Sumw2();
    hSqrtSumWeightsB = new TH1D("hSqrtSumWeightsB","sqrt of sum of weights squares in subevent B", 240, 0.0, 60.0);
    hSqrtSumWeightsB->Sumw2();

    hSqrtSumWeightsNonuni = new TH1D("hSqrtSumWeightsNonuni","sqrt of sum of weights squares", 240, 0.0, 60.0);
    hSqrtSumWeightsNonuni->Sumw2();
    hSqrtSumWeightsANonuni = new TH1D("hSqrtSumWeightsANonuni","sqrt of sum of weights squares in subevent A", 240, 0.0, 60.0);
    hSqrtSumWeightsANonuni->Sumw2();
    hSqrtSumWeightsBNonuni = new TH1D("hSqrtSumWeightsBNonuni","sqrt of sum of weights squares in subevent B", 240, 0.0, 60.0);
    hSqrtSumWeightsBNonuni->Sumw2();

    for (i=0; i<5; i++){

        hCosPhi[i] = new TH1D(Form("hCosPhi%02i",i+1),Form("hCosPhi%02i",i+1),401,-1.0,1.0);
        hCosPhi[i]->Sumw2();
        hSinPhi[i] = new TH1D(Form("hSinPhi%02i",i+1),Form("hSinPhi%02i",i+1),401,-1.0,1.0);
        hSinPhi[i]->Sumw2();

        hCosPhi2[i] = new TH1D(Form("hCosPhi2%02i",i+1),Form("hCosPhi2%02i",i+1),401,-1.0,1.0);
        hCosPhi2[i]->Sumw2();
        hSinPhi2[i] = new TH1D(Form("hSinPhi2%02i",i+1),Form("hSinPhi2%02i",i+1),401,-1.0,1.0);
        hSinPhi2[i]->Sumw2();

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

        //FOR TESTING
        hEPcorrealtion[i] = new TH2D(Form("hEPcorrelation%02i", i+1),Form("hEPcorrelation%02i", i+1),50,-2.0,2.0,50,-2.0,2.0);
        hEPcorrealtionCorr[i] = new TH2D(Form("hEPcorrelationCorr%02i", i+1),Form("hEPcorrelationCorr%02i", i+1),50,-2.0,2.0,50,-2.0,2.0);
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
