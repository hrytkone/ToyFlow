/*
 * Iterointi Newtonin menetelmällä
 */

#include "TMath.h"
#include <boost/math/special_functions/bessel.hpp>

double Rk(double khi, int n) {
    double k = n;
    double bessel = boost::math::cyl_bessel_i((k-1)/2, khi*khi/2) + boost::math::cyl_bessel_i((k+1)/2, khi*khi/2);
    return (TMath::Sqrt(TMath::Pi())/2)*khi*TMath::Exp(-khi*khi/2)*bessel;
}

double func(double *x, double *p) {
    double khi = x[0];
    double k = p[0];
    double Rk = p[1];
    double bessel = boost::math::cyl_bessel_i((k-1)/2, khi*khi/2) + boost::math::cyl_bessel_i((k+1)/2, khi*khi/2);
    return (TMath::Sqrt(TMath::Pi())/2)*khi*TMath::Exp(-khi*khi/2)*bessel - Rk;
}

double RkIter(double x0, double R0, int n, double err) {
    try {
        double x = 0;
        TF1 *fRes = new TF1("fRes", func, 0, 50.0, 2);
        fRes->SetParameters(n, R0);
        while (TMath::Abs(Rk(x, n) - R0) > err) {
            x = x0 - fRes->Eval(x0)/fRes->Derivative(x0);
            x0 = x;
        }
        return x;
    } catch (const std::overflow_error& e) {
        cout << "overflow_error: set khi=0 in case n=" << n << "\n";
        return 0;
    }
}

// Virheen yleisellä etenemisellä R(khi):n lausekkeesta
double CalculateRerror(double khi, double khiErr, double k) {
    double sum1 = (1-2*khi*khi)*(boost::math::cyl_bessel_i((k-1)/2, khi*khi/2) + boost::math::cyl_bessel_i((k+1)/2, khi*khi/2));
    double sum2 = (khi/2)*(boost::math::cyl_bessel_i((k-3)/2, khi*khi/2) + boost::math::cyl_bessel_i((k+1)/2, khi*khi/2) + boost::math::cyl_bessel_i((k-1)/2, khi*khi/2) + boost::math::cyl_bessel_i((k+3)/2, khi*khi/2));
    return TMath::Sqrt(TMath::Pi()/2)*TMath::Exp(-khi*khi)*(sum1 + sum2)*khiErr;
}
