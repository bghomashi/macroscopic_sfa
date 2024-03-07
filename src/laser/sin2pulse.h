#pragma once

#include "sfa/sfa.h"

namespace SFA {
    class Sin2Pulse : public Pulse {
        double period;
    public:
        typedef std::shared_ptr<Sin2Pulse> Ptr_t;

        double delay;
        double CEP;
        double N;
        double w0;
        double E0;
        dvec2 major_pol;

        Sin2Pulse();
        Sin2Pulse(double cep, double w0, double N, double E0, const dvec2& pol, double delay);
        dvec2 A(double t) const;
        dvec2 E(double t) const;
        dvec2 intA(double t) const;
        double intAA(double t) const;
    };
}