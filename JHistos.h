/*
 *  JHistos.h
 *
 */

#include "TH1D.h"
#include "TH2D.h"

class JHistos {

public:
	JHistos();
	virtual ~JHistos(){;}

	TH1D *hPt;
    TH1D *hMultiplicity;
    TH1D *hPhi;

    TH1D *hSqrtSumWeights;
    TH1D *hSqrtSumWeightsA;
    TH1D *hSqrtSumWeightsB;

    TH1D *hPtNonuni;
    TH1D *hMultiplicityNonuni;
    TH1D *hPhiNonuni;

    TH1D *hSqrtSumWeightsNonuni;
    TH1D *hSqrtSumWeightsANonuni;
    TH1D *hSqrtSumWeightsBNonuni;

    /**TH1D *hCosPhi[5];
    TH1D *hSinPhi[5];

    TH1D *hCosPhi2[5];
    TH1D *hSinPhi2[5];**/

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

    // FOR TESTING
    TH1D *hPsiA[5];
    TH1D *hPsiB[5];
    TH1D *hPsiAPsiB[5];
    TH1D *hPsiAcorrected[5];
    TH1D *hPsiBcorrected[5];
    TH1D *hPsiAPsiBcorrected[5];

    TH1D *hQ[5];
    TH1D *hQA[5];
    TH1D *hQB[5];
    TH1D *hQcorr[5];
    TH1D *hQAcorr[5];
    TH1D *hQBcorr[5];
};
