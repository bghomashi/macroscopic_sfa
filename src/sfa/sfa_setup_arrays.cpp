#include "sfa.h"

#include "laser/sin2pulse.h"
#include <iostream>

namespace SFA {
    void SFA::SetupDTM() {
        dtm_el.resize(ps.size());
        for (int ip = 0; ip < ps.size(); ip++) {
            double l = ps[ip];
            dtm_el[ip] = dtm->value({ l, 0 }).x;
            if (l > 0)
                dtm_el[ip] /= l;
        }
    }

    void SFA::SetupFieldArrays() {
        dvector dots(ts.size());
        A = d2vector(ts.size(), { 0,0 });
        E = d2vector(ts.size(), { 0,0 });
        intA = d2vector(ts.size(), { 0,0 });
        intAA = dvector(ts.size(), 0);

        for (auto& p : pulses) {
            for (int i = 0; i < ts.size(); i++) {
                double t = ts[i];
                A[i] += p->A(t);               // compute the total A-field
                E[i] += p->E(t);               // compute the total E-field
                intA[i] += p->intA(t);         // compute the total component-wise integral of A-field
                
                // THIS NEEDS TO BE FIXED
                intAA[i] += p->intAA(t);         // this assumes orthogonal pulses
                //dots[i] = Dot(A[i], A[i]);     // store A.A for integration
                std::cout << A[i].x << " " << A[i].y << std::endl;
            }
            // intAA = TrapzInd<double>(ts, dots);
        }
    }
}