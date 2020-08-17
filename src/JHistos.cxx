#include "JHistos.h"
#include "TMath.h"
#include "TDirectory.h"

JHistos::JHistos(){

    int i, j, k;
    double pi = TMath::Pi();

    int NBINS = 150;
    int N_FOLDERS = 150; //As many as needed.
    double LogBinsX[NBINS+1], LimL=0.1, LimH=100;
    double logBW = (TMath::Log(LimH)-TMath::Log(LimL))/NBINS;
    for (i=0; i<=NBINS; i++) {
        LogBinsX[i] = LimL*exp(i*logBW);
    }

    sHistoBaseNames[nhPt]                     = "hPt";
    sHistoBaseNames[nhPhi]                    = "hPhi";
    sHistoBaseNames[nhCentrality]             = "hCentrality";
    sHistoBaseNames[nhEta]                    = "hEta";
    sHistoBaseNames[nhMultiplicity]           = "hMultiplicity";
    sHistoBaseNames[nhMultiPerDet]            = "hMultiPerDet";
    sHistoBaseNames[nhSqrtSumWeights]         = "hSqrtSumWeights";
    sHistoBaseNames[nhRtrue]                  = "hRtrue";
    sHistoBaseNames[nhRsub]                   = "hRsub";
    sHistoBaseNames[nhVnObs]                  = "hVnObs";
    sHistoBaseNames[nhQnQnAEP]                = "hQnQnAEP";
    sHistoBaseNames[nhQnAQnBEP]               = "hQnAQnBEP";
    sHistoBaseNames[nhQnQnASP]                = "hQnQnASP";
    sHistoBaseNames[nhQnAQnBSP]               = "hQnAQnBSP";
    sHistoBaseNames[nhQvec]                   = "hQvec";
    sHistoBaseNames[nhRsubAB]                 = "hRsubAB";
    sHistoBaseNames[nhRsubAC]                 = "hRsubAC";
    sHistoBaseNames[nhRsubBC]                 = "hRsubBC";
    sHistoBaseNames[nhQnQnAPtBinned]          = "hQnQnAPtBinned";
    sHistoBaseNames[nhSqrtSumWeightsPtBinned] = "hSqrtSumWeightsPtBinned";
    sHistoBaseNames[nhV2ComplexPart]          = "hV2ComplexPart";
    sHistoBaseNames[nhCorrectionParameters]   = "hCorrectionParameters";
    sHistoBaseNames[nhInputNumbers]           = "hInputNumbers";

    TDirectory *fDir = (TDirectory*)gDirectory->mkdir("toyflow");
    TDirectory *fSubDir[N_TOTAL_HISTOS];
    for (i=0; i<N_TOTAL_HISTOS; i++){
        fSubDir[i] = fDir->mkdir(sHistoBaseNames[i]);
    }

    hPt = new TH1D(Form("%s",sHistoBaseNames[nhPt].Data()),"pT - inclusive", NBINS, LogBinsX);
    hPt->Sumw2();
    hPt->SetDirectory(fSubDir[nhPt]);
    hPhi = new TH1D(Form("%s",sHistoBaseNames[nhPhi].Data()), "phi - uniform", 129, -pi, pi);
    hPhi->Sumw2();
    hPhi->SetDirectory(fSubDir[nhPhi]);
    hCentrality = new TH1D(Form("%s",sHistoBaseNames[nhCentrality].Data()), "centrality", CENTBINS_N, 0.0, 90.0);
    hCentrality->Sumw2();
    hCentrality->SetDirectory(fSubDir[nhCentrality]);
    hEta = new TH1D(Form("%s",sHistoBaseNames[nhEta].Data()), "pseudorapidity", 401, -3.8, 5.2);
    hEta->Sumw2();
    hEta->SetDirectory(fSubDir[nhEta]);


    for (j=0; j<CENTBINS_N; j++){
        hMultiplicity[j] = new TH1D(Form("%sCENT%02i",sHistoBaseNames[nhMultiplicity].Data(),j), "Total multiplicity", 300, 0.0, 30000.);
        hMultiplicity[j]->Sumw2();
        hMultiplicity[j]->SetDirectory(fSubDir[nhMultiplicity]);
        for (i=0; i<DET_N; i++) {
            hMultiPerDet[i][j] = new TH1D(Form("%sD%02iCENT%02i",sHistoBaseNames[nhMultiPerDet].Data(),i,j), "Detector multiplicity", 300, 0.0, 30000.);
            hMultiPerDet[i][j]->Sumw2();
            hMultiPerDet[i][j]->SetDirectory(fSubDir[nhMultiPerDet]);
            hSqrtSumWeights[i][j]= new TH1D(Form("%sD%02iCENT%02i",sHistoBaseNames[nhSqrtSumWeights].Data(),i,j),"sqrt of sum of weights squares", 480, 0.0, 120.0);
            hSqrtSumWeights[i][j]->Sumw2();
            hSqrtSumWeights[i][j]->SetDirectory(fSubDir[nhSqrtSumWeights]);
        }
    }
    

    for (i=0; i<nCoef; i++){
        for (j=0; j<DET_N; j++){
            for (k=0; k<CENTBINS_N; k++){
                hRtrue[i][j][k] = new TH1D(Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhRtrue].Data(),i+1,j,k),Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhRtrue].Data(),i+1,j,k),404,-1.01,1.01);
                hRtrue[i][j][k]->Sumw2();
                hRtrue[i][j][k]->SetDirectory(fSubDir[nhRtrue]);
                hRsub[i][j][k] = new TH1D(Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhRsub].Data(),i+1,j,k),Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhRsub].Data(),i+1,j,k),404,-1.01,1.01);
                hRsub[i][j][k]->Sumw2();
                hRsub[i][j][k]->SetDirectory(fSubDir[nhRsub]);
                hVnObs[i][j][k] = new TH1D(Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhVnObs].Data(),i+1,j,k),Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhVnObs].Data(),i+1,j,k),401,-1.5,1.5);
                hVnObs[i][j][k]->Sumw2();
                hVnObs[i][j][k]->SetDirectory(fSubDir[nhVnObs]);

                hQnQnAEP[i][j][k] = new TH1D(Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQnQnAEP].Data(),i+1,j,k),Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQnQnAEP].Data(),i+1,j,k),401,-50.0,50.0);
                hQnQnAEP[i][j][k]->Sumw2();
                hQnQnAEP[i][j][k]->SetDirectory(fSubDir[nhQnQnAEP]);
                hQnAQnBEP[i][j][k] = new TH1D(Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQnAQnBEP].Data(),i+1,j,k),Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQnAQnBEP].Data(),i+1,j,k),404,-1.01,1.01);
                hQnAQnBEP[i][j][k]->Sumw2();
                hQnAQnBEP[i][j][k]->SetDirectory(fSubDir[nhQnAQnBEP]);

                hQnQnASP[i][j][k] = new TH1D(Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQnQnASP].Data(),i+1,j,k),Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQnQnASP].Data(),i+1,j,k),401,-10.0,110.0);
                hQnQnASP[i][j][k]->Sumw2();
                hQnQnASP[i][j][k]->SetDirectory(fSubDir[nhQnQnASP]);
                hQnAQnBSP[i][j][k] = new TH1D(Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQnAQnBSP].Data(),i+1,j,k),Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQnAQnBSP].Data(),i+1,j,k),401,-10.0,110.0);
                hQnAQnBSP[i][j][k]->Sumw2();
                hQnAQnBSP[i][j][k]->SetDirectory(fSubDir[nhQnAQnBSP]);

                hQvec[i][j][k] = new TH2D(Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQvec].Data(),i+1,j,k),Form("%sH%02iD%02iCENT%02i",sHistoBaseNames[nhQvec].Data(),i+1,j,k),101,-100.0,100.0, 101, -100.0, 100.0);
                hQvec[i][j][k]->Sumw2();
                hQvec[i][j][k]->SetDirectory(fSubDir[nhQvec]);
            }
        }
    }

    for (i=0; i<nCoef; i++) {
        for (int j=0; j<CENTBINS_N; j++) {
            hRsubAB[i][j] = new TH1D(Form("%sH%dCENT%02d",sHistoBaseNames[nhRsubAB].Data(), i+1, j), Form("%sH%dCENT%02d",sHistoBaseNames[nhRsubAB].Data(), i+1, j), 100, -1.0, 1.0);
            hRsubAB[i][j]->Sumw2();
            hRsubAB[i][j]->SetDirectory(fSubDir[nhRsubAB]);
            hRsubAC[i][j] = new TH1D(Form("%sH%dCENT%02d",sHistoBaseNames[nhRsubAC].Data(), i+1, j), Form("%sH%dCENT%02d",sHistoBaseNames[nhRsubAC].Data(), i+1, j), 100, -1.0, 1.0);
            hRsubAC[i][j]->Sumw2();
            hRsubAC[i][j]->SetDirectory(fSubDir[nhRsubAC]);
            hRsubBC[i][j] = new TH1D(Form("%sH%dCENT%02d",sHistoBaseNames[nhRsubBC].Data(), i+1, j), Form("%sH%dCENT%02d",sHistoBaseNames[nhRsubBC].Data(), i+1, j), 100, -1.0, 1.0);
            hRsubBC[i][j]->Sumw2();
            hRsubBC[i][j]->SetDirectory(fSubDir[nhRsubBC]);
        }
    }

    for (i=0; i<PTBINS_N; i++) {
        hQnQnAPtBinned[i] = new TH1D(Form("%sPtBin%02i",sHistoBaseNames[nhQnQnAPtBinned].Data(), i+1), Form("%sPtBin%02i",sHistoBaseNames[nhQnQnAPtBinned].Data(), i+1), 482, -11.0, 11.0);
        hQnQnAPtBinned[i]->Sumw2();
        hQnQnAPtBinned[i]->SetDirectory(fSubDir[nhQnQnAPtBinned]);
        hSqrtSumWeightsPtBinned[i] = new TH1D(Form("%sPtBin%02i",sHistoBaseNames[nhSqrtSumWeightsPtBinned].Data(), i+1),Form("sqrt of sum of weights squares for pT bin %02i", i+1), 240, 0.0, 60.0);
        hSqrtSumWeightsPtBinned[i]->Sumw2();
        hSqrtSumWeightsPtBinned[i]->SetDirectory(fSubDir[nhSqrtSumWeightsPtBinned]);
    }

    //FOR TESTING
    hV2ComplexPart = new TH1D(Form("%s",sHistoBaseNames[nhV2ComplexPart].Data()), Form("%s",sHistoBaseNames[nhV2ComplexPart].Data()), 401, -5.0, 5.0);
    hV2ComplexPart->Sumw2();
    hV2ComplexPart->SetDirectory(fSubDir[nhV2ComplexPart]);

    // Save correction parameters
    hCorrectionParameters = new TH2D(Form("%s",sHistoBaseNames[nhCorrectionParameters]), Form("%s",sHistoBaseNames[nhCorrectionParameters]), 5, 0.5, 5.5, 8, 0.5, 8.5);
    hCorrectionParameters->SetDirectory(fSubDir[nhCorrectionParameters]);

    hInputNumbers = new TH1D("hInputNumbers","hInputNumbers",14, 0.5, 14.5);
    hInputNumbers->SetDirectory(fSubDir[nhInputNumbers]);
}
