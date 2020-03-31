#include "JInputs.h"

JInputs::JInputs() {;}

JInputs::~JInputs() {
    for (int i=0; i<CENTBINS_N-2; i++) {
        if(hEtaDist[i]!=0x0) delete hEtaDist[i];
        if(hPtDist[i]!=0x0) delete hPtDist[i];
    }
}

// Initializes the eta distribution histogram
// and the total multiplicities for each centrality
// bin.
//
// Initializes also pt distributions for each centrality bin (0-80%)
void JInputs::Load() {

    int i, j;
    for (i=0; i<CENTBINS_N-2; i++) {
        TFile *fPtInput = TFile::Open("pt-spectra.root", "read");
        // Initialize eta histo:
        hEtaDist[i] = new TH1F(Form("hEtaDist%d", i), Form("hEtaDist%d", i), ETADST_N-1, etadst);
        for (j=0; j<ETADST_N; j++) {
            hEtaDist[i]->SetBinContent(j,etanch[i][j]);
        }
        dMulti[i] = hEtaDist[i]->Integral("width");
        cout << "CentBin: " << i << ", Multi: " << dMulti[i] << endl;

        // Initialize pt histo:
        hPtDist[i] = new TH1F(Form("hPtDist%d", i), Form("hPtDist%d", i), PTDST_N-1, ptdst);
        gPtDist = (TGraphErrors*)fPtInput->Get(Form("Table 2/Graph1D_y%d", i+1));
        for (j=0; j<PTDST_N-1; j++) {
            double x, y;
            gPtDist->GetPoint(j, x, y);
            hPtDist[i]->SetBinContent(j, y);
        }
    }
}

bool JInputs::CheckCentBin(int centBin) {
    return centBin<CENTBINS_N-2 && centBin>-1;
}

int JInputs::GetMultiplicity(int centrality) {
    double multi;
    if(CheckCentBin(centrality)) {
        multi = dMulti[centrality];
    } else {
        multi = 0;
    }
    return (int)TMath::Floor(multi);
}

double JInputs::GetEta(double centrality) {
    double eta;
    int iBin = GetCentBin(centrality);

    if(CheckCentBin(iBin)) {
        eta = hEtaDist[iBin]->GetRandom();
    } else {
        eta = -999; //ok?
    }
    return eta;
}

int JInputs::GetCentBin(double centrality) {
    for (int i=0; i<CENTBINS_N-2; i++)
        if (centrality>centBins[i] && centrality<centBins[i+1]) return i;
    return -1;
}

double JInputs::GetCentDependVn(int n, double centrality) {
    int centBin = GetCentBin(centrality);
    return centvn[n-1][centBin];
}

double JInputs::GetPt(double centrality) {
    int iBin = GetCentBin(centrality);
    return hPtDist[iBin]->GetRandom();
}
