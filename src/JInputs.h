#ifndef JINPUTS_H
#define JINPUTS_H

#include "JConst.h"

#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TFile.h>

using namespace std;

#define ETADST_N 36

const double etadst[ETADST_N] = {
	-3.8, //underflow
	-3.375,-3.125,-2.875,-2.625,-2.375,-2.125,-1.875,-1.625,
    -1.375,-1.125,-0.875,-0.625,-0.375,-0.125, 0.125, 0.375,
     0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 2.375,
     2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125, 4.375,
     4.625, 4.875,
	 5.2 //overflow
};

// Data from arXiv:1612.08966 (10.17182/hepdata.78365)
// The first and the last bins are not data.
// Check sources on this paper: arXiv:1509.07299
const double etanch[CENTBINS_N-1][ETADST_N] = {
	{1640,1643,1670,1718,1787,1835,1912,1968,2001,2021,2017,1995,1970,1943,1929,1929,1943,1970,1995,2017,2021,2001,1968,1912,1835,1787,1718,1670,1643,1563,1474,1370,1324,1281,1244,1240},
	{1360,1364,1391,1424,1474,1507,1569,1644,1679,1682,1672,1646,1621,1597,1583,1583,1597,1621,1646,1672,1682,1679,1644,1569,1507,1474,1424,1391,1364,1292,1218,1132,1093,1062,1032,1030},
	{1000,1038,1061,1080,1114,1136,1178,1229,1253,1256,1247,1229,1210,1191,1181,1181,1191,1210,1229,1247,1256,1253,1229,1178,1136,1114,1080,1061,1038,977,921.3,857.7,829.6,807.4,787,780},
	{700,714,726,738,759,772,797,827,842,844,838,826,811.9,799.2,792.4,792.4,799.2,811.9,826,838,844,842,827,797,772,759,738,726,714,665,625.4,582.6,565.5,551.4,538,530},
	{460,475,482.7,489.7,502.6,510.6,522,539.9,549,549.3,545.5,537.5,527.6,519.3,514.7,514.7,519.3,527.6,537.5,545.5,549.3,549,539.9,522,510.6,502.6,489.7,482.7,475,440,413.6,386.7,375.6,368,359.9,350},
	{300,302,306.3,310.1,317.9,322.3,327.6,335.1,340,340.2,337.7,332.5,326.3,320.7,317.5,317.5,320.7,326.3,332.5,337.7,340.2,340,335.1,327.6,322.3,317.9,310.1,306.3,302,277.5,261.3,244.7,238.4,233.8,229.4,225},
	{177,178,179.9,181.7,186,188.2,189.8,193.5,196.4,196.5,194.8,191.4,187.5,184.3,182.5,182.5,184.3,187.5,191.4,194.8,196.5,196.4,193.5,189.8,188.2,186,181.7,179.9,178,163.2,153.4,143.8,140.3,138.7,136,135},
	{93,94.9,96.1,96.8,98.3,98.8,99.1,101.2,102.7,103.1,102,100.3,98,96.1,95.2,95.2,96.1,98,100.3,102,103.1,102.7,101.2,99.1,98.8,98.3,96.8,96.1,94.9,86.8,81.9,77.3,75.8,75.1,73.8,73},
    {45,45.4,46.0,46.0,46.5,46.8,46.7,47.3,47.9,48.2,47.6,46.7,45.5,44.6,44.1,44.1,44.6,45.5,46.7,47.6,48.2,47.9,47.3,46.7,46.8,46.5,46.0,46.0,45.4,41.4,39.3,37.4,36.8,36.6,35.7,35},
    {18,18.7,18.7,18.6,18.7,18.6,18.6,19.2,19.1,19.2,18.9,18.5,18.0,17.6,17.4,17.4,17.6,18.0,18.5,18.9,19.2,19.1,19.2,18.6,18.6,18.7,18.6,18.7,18.7,17.1,16.27,15.46,15.47,15.6,15.4,15}
};

class JInputs {

public:
    JInputs();
    virtual ~JInputs(){;}

    void Load();
    // double -> int
    int GetMultiplicity(double centrality) { return hEtaDist[centrality]->Integral();}//{ return gNch[detId]->Eval(centrality); }
    double GetEta(int centrality) { return hEtaDist[centrality]->GetRandom();}

private:
    TGraph *gNch[DET_N];
    TH1F *hEtaDist[CENTBINS_N-1];
    TRandom3 *rand;
};

#endif
