#pragma once

#include <complex>
#include <cassert>
#include <vector>

using namespace std::complex_literals;          // enables 'i' in complex literals: complex x = 2.5i;

typedef std::complex<double> complex;
typedef std::vector<complex> cvector;
typedef std::vector<double> dvector;

const double Pi2        = 6.283185307179586476925286766559;             // 2pi
const double Pi         = 3.1415926535897932384626433832795;            // pi
const double Pi_2       = 1.5707963267948966192313216916398;            // pi/2
const double inv2Pi     = 0.15915494309189533576888376337251;           // 1/(2Pi)
const double invPi_2    = 0.63661977236758134307553505349006;           // 1/(pi/2)
const double sqrtPi2    = 2.506628274631000502415765284811;             // sqrt(2*pi)
const double sqrtPi_2   = 1.2533141373155002512078826424055;            // sqrt(pi/2)
const double C          = 137.0;                    // speed of light in AU
const double LnmToAU    = 45.5513500662;            // wavelength in nm to AU (wavelength in *denominator*)
const double nmToAU     = 1. / 0.052917721092;      // length in nm to AU
const double cmToAU     = 1. / 0.52917721092e-8;    // length in cm to AU
const double Wcm2ToAU   = 1. / 3.51e16;             // 
const double umToAU     = 1.89393939394e4;          // length in um to AU


inline double Dot(double a, double b) {
    return a * b;
}
inline complex Dot(complex a, complex b) {
    return a * std::conj(b);
}

inline dvector Range(double df, double fmax, double fmin = 0) {
    int nf = int((fmax - fmin) / df) + 1;
    dvector f(nf);
    for (int i = 0; i < nf; i++)
        f[i] = fmin + i * df;
    return f;
}
