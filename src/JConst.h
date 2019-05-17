#ifndef JCONST_H
#define JCONST_H

#define PI TMath::Pi()

const int nCoef = 5;
const int nCorrParam = 8;

// for pT-dependence
const double Tdec = 0.12;
const double vr = 0.6;

#define SECTORS_N 8
#define RINGS_N 5

static double ringEta[RINGS_N+1] = {2.20, 2.81, 3.41, 3.9, 4.47, 5.06};

enum DETECTOR {
    D_TPC,
    D_TPC_A,
    D_TPC_C,
    D_VOP,
    DET_N
};

static double cov[DET_N][2] = {
    {-1.5, 1.5},
    {-1.5, -0.4},
    {0.4, 1.5},
    {2.2, 5.06}
};

#define PTBINS_N 9
static double pTBins[PTBINS_N+1] = {0.0, 0.2, 0.6, 1.2, 2.0, 3.0, 4.2, 5.6, 7.2, 9.0};

#define CENTBINS_N 9
static double centBins[CENTBINS_N] = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};

#endif
