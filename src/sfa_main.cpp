
#include "maths/maths.h"
#include "maths/vec2.h"
#include "sfa/sfa.h"
#include "laser/sin2pulse.h"
#include "gas_jet/gas_jet.h"
#include "laser/laser.h"
#include "dtms/zero_range_dtm.h"
#include "utility/save_data.h"
#include "utility/profiler.h"
#include "utility/thread_pool.h"
#include "utility/logger.h"
#include "utility/json.hpp"
#include <iostream>
#include <fstream>
#include "common_params.h"

int SFA_MAIN() {
    auto frequencies = Range(dff, ffmax, ffmin);
    c2vector spectrum(frequencies.size());

    // ----------- set up sfa ------------------------
    LOG_INFO("set up sfa");
    SFA::SFA sfa;
    sfa.ts = Range(dt, tmax);
    std::cout <<sfa.ts.size() << std::endl;

// for (int i = 19900; i < 20000; i++)
// std::cout << i <<" "<<sfa.ts[i] << std::endl;
// exit(0);
    sfa.ps = Range(dp, pmax, pmin);
    sfa.frequencies = Range(dff, ffmax, ffmin);
    sfa.Ip = 0.5;
    sfa.dtm = std::make_shared<SFA::ZeroRange>(sfa.Ip);
    std::vector<SFA::Sin2Pulse::Ptr_t> pulses(peak_I0_wcm2.size());
    
    for (int i = 0; i < peak_I0_wcm2.size(); i++) {
        pulses[i] = std::make_shared<SFA::Sin2Pulse>(0., w0[i], Ncyc[i], 0., major_pol[i], cycles_delay[i]);
        pulses[i]->CEP = cep[i];
        pulses[i]->E0 = sqrt(peak_I0_wcm2[i]*Wcm2ToAU);
    }
    sfa.pulses.insert(sfa.pulses.end(), pulses.begin(), pulses.end());
    sfa.SetupVectorization();
    sfa.SetupDTM();

    LOG_INFO("SetupFieldArrays");
    sfa.SetupFieldArrays();
    LOG_INFO("Execute2D");
    sfa.Execute2D();
    LOG_INFO("Spectrum");
    sfa.Spectrum();

    //sfa.FreeVectorization();


    ThreadPool::Shutdown();

    Store(output_filename, frequencies, sfa.hhg);

    return 0;
}