/*
 *  JHistos.h
 *
 */

#include "TH1D.h"
#include "TH2D.h"

#include "JConst.h"

class JHistos {

public:
	JHistos();
	virtual ~JHistos(){;}

	TH1D *hPt;
    TH1D *hPhi;
    TH1D *hPtNonuni;
    TH1D *hPhiNonuni;
    TH1D *hCentrality;
    TH1D *hEta;

    TH1D *hMultiplicity;
    TH1D *hMultiplicityNonuni;

    TH1D *hSqrtSumWeights[DET_N];
    TH1D *hSqrtSumWeightsNonuni[DET_N];

    // Historgrams for resolutions and vobs
    TH1D *hRtrue[nCoef][DET_N];
    TH1D *hRsub[nCoef][DET_N];
    TH1D *hVnObs[nCoef][DET_N];

    TH1D *hRtrueNonuni[nCoef][DET_N];
    TH1D *hRsubNonuni[nCoef][DET_N];
    TH1D *hVnObsNonuni[nCoef][DET_N];

    // Histograms for EP-method
    TH1D *hQnQnAEP[nCoef][DET_N];
    TH1D *hQnAQnBEP[nCoef][DET_N];

    TH1D *hQnQnAEPnonuni[nCoef][DET_N];
    TH1D *hQnAQnBEPnonuni[nCoef][DET_N];

    // Histograms for SP-method
    TH1D *hQnQnASP[nCoef][DET_N];
    TH1D *hQnAQnBSP[nCoef][DET_N];

    TH1D *hQnQnASPnonuni[nCoef][DET_N];
    TH1D *hQnAQnBSPnonuni[nCoef][DET_N];

    // pT bins for v2
    TH1D *hQnQnAPtBin[PTBINS_N];
    TH1D *hSqrtSumWeightsPtBins[PTBINS_N];

    //FOR TESTING
    TH1D *hV2ComplexPart;

};
