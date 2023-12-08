#pragma once

#include "sfa/sfa.h"

namespace SFA {
    class ZeroRange : public DTM {
        double Ip;
    public:
        ZeroRange(double Ip);
        cvec2 value(const dvec2& p) const;
    };
}