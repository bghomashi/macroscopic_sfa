#include "sfa.h"
#include "maths/fast_maths.h"
#include "utility/profiler.h"
#include "utility/logger.h"

namespace SFA {
    void SFA::SetupVectorization() {
        int NTs = ts.size() * (ts.size() + 1) / 2 * sizeof(double);
        LOG_INFO("allocating: " + std::to_string(13*NTs + (ts.size() * sizeof(cvec2))) + " bytes (" + std::to_string(double(13*NTs + (ts.size() * sizeof(cvec2)))/1024./1024./1024.) + " Gb" + ")");
#if defined(PROFILING)
        Profile::Push("SFA::SetupVectorization");
#endif
        dip = (cvec2*)fast::aligned_malloc(ts.size() * sizeof(cvec2), 32);
        realx = (double*)fast::aligned_malloc(NTs, 32);
        imagx = (double*)fast::aligned_malloc(NTs, 32);
        realy = (double*)fast::aligned_malloc(NTs, 32);
        imagy = (double*)fast::aligned_malloc(NTs, 32);

        dtm_iti_real_x = (double*)fast::aligned_malloc(NTs, 32);
        dtm_iti_imag_x = (double*)fast::aligned_malloc(NTs, 32);
        dtm_iti_real_y = (double*)fast::aligned_malloc(NTs, 32);
        dtm_iti_imag_y = (double*)fast::aligned_malloc(NTs, 32);

        dtm_itr_real_x = (double*)fast::aligned_malloc(NTs, 32);
        dtm_itr_imag_x = (double*)fast::aligned_malloc(NTs, 32);
        dtm_itr_real_y = (double*)fast::aligned_malloc(NTs, 32);
        dtm_itr_imag_y = (double*)fast::aligned_malloc(NTs, 32);

        S0_saddle = (double*)fast::aligned_malloc(NTs, 32);

#if defined(PROFILING)
        Profile::Pop("SFA::SetupVectorization");
#endif
    }

    void SFA::FreeVectorization() {
        fast::aligned_free(dip);
        fast::aligned_free(realx);
        fast::aligned_free(imagx);
        fast::aligned_free(realy);
        fast::aligned_free(imagy);
        fast::aligned_free(dtm_iti_real_x);
        fast::aligned_free(dtm_iti_imag_x);
        fast::aligned_free(dtm_iti_real_y);
        fast::aligned_free(dtm_iti_imag_y);
        fast::aligned_free(dtm_itr_real_x);
        fast::aligned_free(dtm_itr_imag_x);
        fast::aligned_free(dtm_itr_real_y);
        fast::aligned_free(dtm_itr_imag_y);
        fast::aligned_free(S0_saddle);

    }
}