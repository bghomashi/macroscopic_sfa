#include "sfa.h"
#include "maths/maths.h"
#include "maths/vec2.h"
#include "maths/blackman.h"
#include "utility/profiler.h"
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace SFA {
    void SFA::Spectrum() {
#if defined(PROFILING)
        Profile::Push("SFA::Spectrum");
#endif

        int nt = int(ts.size());
        int nf = int(frequencies.size());
        hhg.resize(nf);

        auto window = blackman(nt);



        std::ofstream file("dipole.dat");
        file << std::setprecision(8) << std::scientific;
        for (int itr = 0; itr < nt; itr++) {
            file << ts[itr] << " "  << window[itr] * std::real(dipole[itr].x) << " " << window[itr] * std::real(dipole[itr].y) << std::endl;
        }
        // exit(0);


        double dt = ts[1] - ts[0];
        for (int iff = 0; iff < nf; iff++) {
            double f = frequencies[iff];
            hhg[iff] = cvec2{ 0,0 };
            for (int it = 0; it < nt; it++) {
                double t = ts[it];

                hhg[iff] += std::exp(-1.i * f * t) * window[it] * real(dipole[it]) * dt;
            }
            hhg[iff] *= -f * f;
        }
#if defined(PROFILING)
        Profile::Pop("SFA::Spectrum");
#endif
    }
}