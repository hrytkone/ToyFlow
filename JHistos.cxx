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

    for (i=0; i<5; i++){
        hRtrue[i] = new TH1D(Form("hRtrue%02i",i+1),Form("hRtrue%02i",i+1),401,-1,1);
        hRtrue[i]->Sumw2();

        hCosPhi[i] = new TH1D(Form("hCosPhi%02i",i+1),Form("hCosPhi%02i",i+1),401,-1.0,1.0);
        hCosPhi[i]->Sumw2();
        hSinPhi[i] = new TH1D(Form("hSinPhi%02i",i+1),Form("hSinPhi%02i",i+1),401,-1.0,1.0);
        hSinPhi[i]->Sumw2();

        hCosPhi2[i] = new TH1D(Form("hCosPhi2%02i",i+1),Form("hCosPhi2%02i",i+1),401,-1.0,1.0);
        hCosPhi2[i]->Sumw2();
        hSinPhi2[i] = new TH1D(Form("hSinPhi2%02i",i+1),Form("hSinPhi2%02i",i+1),401,-1.0,1.0);
        hSinPhi2[i]->Sumw2();

        hVnObs[i] = new TH1D(Form("hVnObs%02i",i+1),Form("hVnObs%02i",i+1),401,-1.5,1.5);
        hVnObs[i]->Sumw2();
        hRsub[i] = new TH1D(Form("hRsub%02i",i+1),Form("hRsub%02i",i+1),401,-1.0,1.0);
        hRsub[i]->Sumw2();

        hVnObsCorrected[i] = new TH1D(Form("hVnObsCorrected%02i",i+1),Form("hVnObsCorrected%02i",i+1),401,-1.5,1.5);
        hVnObsCorrected[i]->Sumw2();
        hRsubCorrected[i] = new TH1D(Form("hRsubCorrected%02i",i+1),Form("hRsubCorrected%02i",i+1),401,-1,1);
        hRsubCorrected[i]->Sumw2();
    }

    for (i=0; i<9; i++) {
        hPtBin[i] = new TH1D(Form("hPtBin%02i", i+1), Form("hPtBin%02i", i+1), 401, -10.0, 10.0);
        hPtBin[i]->Sumw2();
    }
}
