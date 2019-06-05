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

    // Historgrams for calculating vn values
    TH1D *hRtrue[nCoef][DET_N];
    TH1D *hRsub[nCoef][DET_N];
    TH1D *hVnObs[nCoef][DET_N];

    TH1D *hRtrueCorrected[nCoef][DET_N];
    TH1D *hRsubCorrected[nCoef][DET_N];
    TH1D *hVnObsCorrected[nCoef][DET_N];

    // Histograms for the alternative event plane method
    TH1D *hQnQnA[nCoef][DET_N];
    TH1D *hQnAQnB[nCoef][DET_N];

    TH1D *hQnQnAcorrected[nCoef][DET_N];
    TH1D *hQnAQnBcorrected[nCoef][DET_N];

    // pT bins for v2
    TH1D *hPtBin[PTBINS_N];
    TH1D *hSqrtSumWeightsPtBins[PTBINS_N];

    //FOR TESTING
    TH1D *hV2ComplexPart;

};
