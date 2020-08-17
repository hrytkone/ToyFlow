/*
 *  JHistos.h
 *
 */

#include "TH1D.h"
#include "TH2D.h"

#include "JConst.h"

// A recipe for adding a new histogram:
// - Add a new histogram object to this header.
// - Add a new entry to the enumerator, this keeps
//   track of how many histograms we have and gives
//   them a number.
// - Give the histogram a name in sHistoBaseNames list
//   in the JHistos.cxx file
// - Initialize the histogram in the cxx file in a similar
//   manner as previously done. (Remember to use Setw2 and 
//   SetDirectory as well)

class JHistos {

public:
	JHistos();
	virtual ~JHistos(){;}

	TH1D *hPt;
    TH1D *hPhi;
    TH1D *hCentrality;
    TH1D *hEta;
    TH1D *hMultiplicity[CENTBINS_N];
    TH1D *hMultiPerDet[DET_N][CENTBINS_N];

    TH1D *hSqrtSumWeights[DET_N][CENTBINS_N];

    // Historgrams for resolutions and vobs
    TH1D *hRtrue[nCoef][DET_N][CENTBINS_N];
    TH1D *hRsub[nCoef][DET_N][CENTBINS_N];
    TH1D *hVnObs[nCoef][DET_N][CENTBINS_N];

    // Histograms for EP-method
    TH1D *hQnQnAEP[nCoef][DET_N][CENTBINS_N];
    TH1D *hQnAQnBEP[nCoef][DET_N][CENTBINS_N];

    // Histograms for SP-method
    TH1D *hQnQnASP[nCoef][DET_N][CENTBINS_N];
    TH1D *hQnAQnBSP[nCoef][DET_N][CENTBINS_N];
    TH2D *hQvec[nCoef][DET_N][CENTBINS_N];

    // pT bins for v2
    TH1D *hQnQnAPtBinned[PTBINS_N];
    TH1D *hSqrtSumWeightsPtBinned[PTBINS_N];

    // 3-sub event
    TH1D *hRsubAB[nCoef][CENTBINS_N];
    TH1D *hRsubAC[nCoef][CENTBINS_N];
    TH1D *hRsubBC[nCoef][CENTBINS_N];

    //FOR TESTING
    TH1D *hV2ComplexPart;
    TH2D *hCorrectionParameters;
    TH1D *hInputNumbers;


private:
    // This enum is used to find the correct names for histos
    enum HISTOS {
         nhPt
        ,nhPhi
        ,nhCentrality
        ,nhEta
        ,nhMultiplicity
        ,nhMultiPerDet
        ,nhSqrtSumWeights
        ,nhRtrue
        ,nhRsub
        ,nhVnObs
        ,nhQnQnAEP
        ,nhQnAQnBEP
        ,nhQnQnASP
        ,nhQnAQnBSP
        ,nhQvec
        ,nhRsubAB
        ,nhRsubAC
        ,nhRsubBC
        ,nhQnQnAPtBinned
        ,nhSqrtSumWeightsPtBinned
        ,nhV2ComplexPart
        ,nhCorrectionParameters
        ,nhInputNumbers
        ,N_TOTAL_HISTOS // This needs to be last in here. Please add new entries before this one.
    };

    TString sHistoBaseNames[N_TOTAL_HISTOS];
};
