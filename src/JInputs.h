#ifndef JINPUTS_H
#define JINPUTS_H

#include "JConst.h"

#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TMath.h>

using namespace std;

#define ETADST_N 36
//#define CENTDST_N 30
#define CENTDST_N 7

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

/**const double centdst[CENTDST_N] = {
    1.0,3.0,5.0,7.0,9.0,11.0,13.0,15.0,17.0,19.0,21.0,23.0,
    25.0,27.0,29.0,31.0,33.0,35.0,37.0,39.0,41.0,43.0,45.0,
    47.0,49.0,51.0,53.0,55.0,57.0,59.0
};**/

// Data from arXiv:1804.02944 (10.17182/hepdata.83737)
// For v2, v3, and v4 every second bin was combined with previous one since v5 has two times bigger bin width
/**const double centvn[nCoef-1][CENTDST_N] = {
    {0.02260806,0.02936416,0.0365414,0.04339233,0.04965241,0.05540649,0.06074488,0.06563346,0.07012421,0.07418512,0.07799644,0.08146661,0.08466233,0.08750616,0.0900406,0.09234899,
    0.0944606,0.09628565,0.09781767,0.09896134,0.09998829,0.10072615,0.1011294,0.1013651,0.10115135,0.10088945,0.10046385,0.09956397,0.09841429,0.09695983},
    {0.0191187,0.02093017,0.02222267,0.02352247,0.02444807,0.02522642,0.02607792,0.02685555,0.02741218,0.027943,0.02864565,0.02898334,0.02958154,0.02989943,0.0303948,0.0307196,
    0.03088819,0.03129652,0.03152673,0.03150934,0.03146594,0.0314677,0.03142536,0.03106437,0.03060843,0.03048753,0.02992579,0.02944655,0.02831613,0.02775181},
    {0.00953763,0.01035041,0.01104624,0.01162182,0.01207197,0.012597,0.01280591,0.01307446,0.01336189,0.01379324,0.01441826,0.01430597,0.01480082,0.01488894,0.01555275,0.01546395,
    0.01579773,0.01568419,0.01659463,0.01631601,0.0167705,0.01604588,0.01673617,0.01647489,0.01620397,0.01566972,0.01487526,0.01594789,0.01557572,0.01468211},
    {0.003721945,0.004261125,0.004491402,0.004781827,0.004765597,0.004676638,0.005054133,0.005370356,0.005561058,0.005749716,0.006059177,0.006337009,0.006657293,0.006606902,
    0.006844629,0.006491805,0.006888308,0.007427039,0.006477826,0.007157546,0.006901006,0.006778134,0.00696109,0.007004513,0.006874222,0.0061236,0.005875487,0.004894928,
    0.005821316,0.006047204}
};**/

// Constant vn in each centrality bin
const double centvn[nCoef][CENTDST_N] = {
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.0277401, 0.04488324, 0.06521883, 0.08433443, 0.09597485, 0.10087206, 0.09925828},
    {0.02039728, 0.02369955, 0.02670301, 0.02950095, 0.03118808, 0.03120636, 0.02918556},
    {0.01013229, 0.01171893, 0.0131265, 0.01479335, 0.0159713, 0.01644628, 0.01535014},
    {0.00415816, 0.00467961, 0.00528238, 0.006501, 0.0068885, 0.00690379, 0.00575251}
};

class JInputs {

public:
    JInputs();
    virtual ~JInputs();

    void Load();
    int GetMultiplicity(int centrality);
    double GetEta(int centrality);
    int GetCentBin(double centrality);
    double GetCentDependVn(int n, double centrality);

private:
    bool CheckCentBin(int centBin);

    TH1F *hEtaDist[CENTBINS_N-1] = {NULL};
    int dMulti[CENTBINS_N-1] = {0};
};

#endif
