#include <iostream>
#include <fstream>
#include "utility/json.hpp"
#include "utility/profiler.h"
#include "utility/logger.h"
#include "maths/vec2.h"
#include "common_params.h"

bool sfa_only = false;
size_t worker_thread = 16;
// laser parameters
dvector peak_I0_wcm2 = {1e14};
dvector waist_um = {30};
dvector Lnm = {800};
dvector g0 = {0};
dvector Ncyc = {10};
dvector w0 = {LnmToAU / Lnm.back()};
d2vector major_pol = {{1., 0.}};
dvector cep = {0};
dvector cycles_delay = {0};

// gas parameters
double gas_radius_um = 500;
double gas_length_um = 6 * waist_um.back();
size_t gas_cells = 100000;
double gas_sig_um = 800;
double gas_density_cm3 = 1;

double dt = 0.2;
double tmax = (2. * Pi / w0.back()) * Ncyc.back();
double dp = 0.00001;
double pmin = 0, pmax = 3;
double dff = 0.01;
double ffmin = 0. * w0.back();
double ffmax = 40. * w0.back();

std::string output_filename = "out.dat";

bool ReadInput(const std::string& filename);
int MACRO_MAIN();
int SFA_MAIN();

int main() {
    Log::set_logger_file("log.txt");
#if defined(PROFILING)
    Profile::Push("total");
#endif
    if (!ReadInput("input.json"))
        LOG_INFO("input.json not found");

    LOG_INFO("Go...");

    if (sfa_only)
        SFA_MAIN();
    else 
        MACRO_MAIN();

    #if defined(PROFILING)
        Profile::Pop("total");
        Profile::Print();
        Profile::PrintTo("profile.txt");
    #endif

}

bool ReadInput(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        return false;
        
    LOG_INFO("reading input file");

    nlohmann::json input;
    file >> input;
    if (input.contains("sfa_only"))
        sfa_only = input["sfa_only"].get<bool>();
    else
        sfa_only = false;

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
    major_pol.clear();
    Lnm.clear();
    w0.clear();
    cep.clear();
    cycles_delay.clear();

    for (auto l : input["lasers"]) {
        peak_I0_wcm2.push_back(l["intensity"].get<double>());
        waist_um.push_back(l["beam_waist_um"].get<double>());
        g0.push_back(l["porras_factor"].get<double>());
        Ncyc.push_back(l["cycles"].get<double>());
        major_pol.push_back({l["polarization"][0].get<double>(), 
                             l["polarization"][1].get<double>()});
        


        if (l.contains("cycles_delay")) {
            cycles_delay.push_back(l["cycles_delay"].get<double>());
        } else {
            cycles_delay.push_back(0);
        }
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

    auto gas_jet = input["gas_jet"];
    // gas parameters
    gas_radius_um = gas_jet["radius_um"].get<double>();
    gas_length_um = gas_jet["length_um"].get<double>();
    gas_cells = gas_jet["cells"].get<size_t>();
    gas_sig_um = gas_jet["sigma_um"].get<double>();
    gas_density_cm3 = gas_jet["density_cm3"].get<double>();

    return true;
}
