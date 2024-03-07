#include "zero_range_dtm.h"

namespace SFA {
    ZeroRange::ZeroRange(double Ip) : Ip(Ip) {

    }
    cvec2 ZeroRange::value(const dvec2& p) const {
        const static double Ip2 = 2. * Ip;
        const static double srtK = sqrt(sqrt(2. * Ip));
        const double mag = 1. / (Dot(p, p) + Ip2);

        return (-2.i * (srtK * mag * mag)) * p / sqrtPi;

    }
}

