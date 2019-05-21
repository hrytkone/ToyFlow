#include "JInputs.h"

JInputs::JInputs() {
    rand = new TRandom3(0);
    rand->SetSeed(0);

    for (int i=0; i<DET_N; i++)
        gNch[i] = new TGraph(CENTBINS_N-1);
}

void JInputs::Load() {
    double nch;
    int i, j;
    for (i=0; i<CENTBINS_N-1; i++) {
        TGraph gEta(ETADST_N, etadst, etanch[i]);
        TF1 fEtach("etach", [&](double *px, double *pp)->double{ return gEta.Eval(px[0]); }, -3.5, 5.1, 0);

        for (j=0; j<DET_N; j++) {
            nch = fEtach.Integral(cov[j][0], cov[j][1])/(cov[j][1]-cov[j][0]);
            gNch[j]->SetPoint(i, 0.5*(centBins[i]+centBins[i+1]),nch);
        }
    }
}

double JInputs::GetEta(int detId) {
    return rand->Uniform(cov[detId][0], cov[detId][1]);
}