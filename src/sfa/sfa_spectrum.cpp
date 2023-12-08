#include "sfa.h"
#include "maths/vec2.h"
#include "maths/blackman.h"
#include "utility/profiler.h"

namespace SFA {
    void SFA::Spectrum() {
#if defined(PROFILING)
        Profile::Push("SFA::Spectrum");
#endif

        int nt = int(ts.size());
        int nf = int(frequencies.size());
        hhg.resize(nf);

        auto window = blackman(nt);

        double dt = ts[1] - ts[0];
        for (int iff = 0; iff < nf; iff++) {
            double f = frequencies[iff];
            hhg[iff] = cvec2{ 0,0 };
            for (int it = 0; it < nt; it++) {
                double t = ts[it];

                hhg[iff] += std::exp(-1i * f * t) * window[it] * real(dipole[it]) * dt;
            }
            hhg[iff] *= -f * f;
        }
#if defined(PROFILING)
        Profile::Pop("SFA::Spectrum");
#endif
    }
}