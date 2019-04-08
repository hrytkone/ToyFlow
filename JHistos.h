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

    TH1D *hPtNonuni;
    TH1D *hMultiplicityNonuni;
    TH1D *hPhiNonuni;

    TH1D *hCosPhi[5];
    TH1D *hSinPhi[5];

    TH1D *hCosPhi2[5];
    TH1D *hSinPhi2[5];

    // Historgrams for calculating vn values
    TH1D *hRtrue[5];

    TH1D *hRsub[5];
    TH1D *hVnObs[5];

    TH1D *hRsubCorrected[5];
    TH1D *hVnObsCorrected[5];

    // pT bins for v2
    TH1D *hPtBin[9];

};
