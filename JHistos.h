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

    TH1D *hMultiplicity[DET_N];
    TH1D *hMultiplicityNonuni[DET_N];

    TH1D *hSqrtSumWeightsTPC;
    TH1D *hSqrtSumWeightsTPCNonuni;
    TH1D *hSqrtSumWeightsTPCA;
    TH1D *hSqrtSumWeightsTPCANonuni;
    TH1D *hSqrtSumWeightsTPCC;
    TH1D *hSqrtSumWeightsTPCCNonuni;
    TH1D *hSqrtSumWeightsV0P;
    TH1D *hSqrtSumWeightsV0PNonuni;

    // Historgrams for calculating vn values
    TH1D *hRtrue[5];
    TH1D *hRsub[5];
    TH1D *hVnObs[5];

    TH1D *hRtrueCorrected[5];
    TH1D *hRsubCorrected[5];
    TH1D *hVnObsCorrected[5];

    // Histograms for the alternative event plane method
    TH1D *hQnQnA[5];
    TH1D *hQnAQnB[5];

    TH1D *hQnQnAcorrected[5];
    TH1D *hQnAQnBcorrected[5];

    // pT bins for v2
    TH1D *hPtBin[9];
    TH1D *hSqrtSumWeightsPtBins[9];

    //FOR TESTING
    TH1D *hV2ComplexPart;

};
