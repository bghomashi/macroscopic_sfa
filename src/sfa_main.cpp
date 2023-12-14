#if !defined(SFA_ONLY)
#define SFA_MAIN sfa_main
#else
#define SFA_MAIN main
#endif

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

static size_t worker_thread = 16;
// laser parameters
static dvector peak_I0_wcm2 = {1e14};
static dvector waist_um = {30};
static dvector Lnm = {800};
static dvector g0 = {0};
static dvector Ncyc = {10};
static dvector w0 = {LnmToAU / Lnm.back()};
static d2vector major_pol = {{1., 0.}};
static dvector cep = {0};

// gas parameters
static double gas_radius_um = 500;
static double gas_length_um = 6 * waist_um.back();
static size_t gas_cells = 100000;
static double gas_sig_um = 800;
static double gas_density_cm3 = 1;

static double dt = 0.2;
static double tmax = (2. * Pi / w0.back()) * Ncyc.back();
static double dp = 0.00001;
static double pmin = 0, pmax = 3;
static double dff = 0.01;
static double ffmin = 0. * w0.back(), ffmax = 40. * w0.back();

static std::string output_filename = "out.dat";

static bool ReadInput(const std::string& filename);

int SFA_MAIN() {
    Log::set_logger_file("log.txt");
#if defined(PROFILING)
    Profile::Push("total");
#endif
    if (!ReadInput("input.json"))
        LOG_INFO("input.json not found");

    auto frequencies = Range(dff, ffmax, ffmin);
    c2vector spectrum(frequencies.size());

    // ----------- set up sfa ------------------------
    LOG_DEBUG("set up sfa");
    SFA::SFA sfa;
    sfa.ts = Range(dt, tmax);
    sfa.ps = Range(dp, pmax, pmin);
    sfa.frequencies = Range(dff, ffmax, ffmin);
    sfa.Ip = 0.5;
    sfa.dtm = std::make_shared<SFA::ZeroRange>(sfa.Ip);
    std::vector<SFA::Sin2Pulse::Ptr_t> pulses(peak_I0_wcm2.size());
    for (int i = 0; i < peak_I0_wcm2.size(); i++) {
        pulses[i] = std::make_shared<SFA::Sin2Pulse>(0., w0[i], Ncyc[i], 0., major_pol[i]);
        pulses[i]->CEP = cep[i];
        pulses[i]->E0 = sqrt(peak_I0_wcm2[i]*Wcm2ToAU);
    }
    sfa.pulses.insert(sfa.pulses.end(), pulses.begin(), pulses.end());
    sfa.SetupVectorization();
    sfa.SetupDTM();

    LOG_DEBUG("SetupFieldArrays");
    sfa.SetupFieldArrays();
    LOG_DEBUG("Execute2D");
    sfa.Execute2D();
    LOG_DEBUG("Spectrum");
    sfa.Spectrum();

    //sfa.FreeVectorization();


    ThreadPool::Shutdown();

    Store(output_filename, frequencies, sfa.hhg);
#if defined(PROFILING)
    Profile::Pop("total");
#endif

#if defined(PROFILING)
    Profile::Print();
    Profile::PrintTo("profile.txt");
#endif
    return 0;
}


bool ReadInput(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        return false;
        
    LOG_INFO("reading input file");

    nlohmann::json input;
    file >> input;

    worker_thread = input["threads"].get<size_t>();
    dt = input["dt"].get<double>();
    dp = input["dp"].get<double>();
    pmin = input["pmin"].get<double>();
    pmax = input["pmax"].get<double>();
    dff = input["df"].get<double>();
    ffmin = input["fmin"].get<double>();
    ffmax = input["fmax"].get<double>();
    output_filename = input["output_filename"].get<std::string>();

    tmax = 0;
    peak_I0_wcm2.clear();
    waist_um.clear();
    g0.clear();
    Ncyc.clear();
    Lnm.clear();
    w0.clear();
    cep.clear();
    for (auto l : input["lasers"]) {
        peak_I0_wcm2.push_back(l["intensity"].get<double>());
        waist_um.push_back(l["beam_waist_um"].get<double>());
        g0.push_back(l["porras_factor"].get<double>());
        Ncyc.push_back(l["cycles"].get<double>());
        major_pol.push_back({l["polarization"][0].get<double>(), l["polarization"][1].get<double>()});
        if (l.contains("wavelength_nm")) {
            Lnm.push_back(l["wavelength_nm"].get<double>());
            w0.push_back(LnmToAU / Lnm.back());
        } else if (l.contains("frequency")) {
            w0.push_back(l["frequency"].get<double>());
            Lnm.push_back(LnmToAU / w0.back());
        }
        tmax = std::max(tmax, Ncyc.back()*2.*Pi/w0.back());
        if (l.contains("cep"))
            cep.push_back(l["cep"].get<double>());
        else
            cep.push_back(0);
    }

    LOG_DEBUG("gas_jet");

    auto gas_jet = input["gas_jet"];
    // gas parameters
    gas_radius_um = gas_jet["radius_um"].get<double>();
    gas_length_um = gas_jet["length_um"].get<double>();
    gas_cells = gas_jet["cells"].get<size_t>();
    gas_sig_um = gas_jet["sigma_um"].get<double>();
    gas_density_cm3 = gas_jet["density_cm3"].get<double>();



    return true;
}