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
double ffmin = 0. * w0.back(), ffmax = 40. * w0.back();

std::string output_filename = "out.dat";

bool ReadInput(const std::string& filename);

int main() {
    Log::set_logger_file("log.txt");
    LOG_INFO("main()");

    if (!ReadInput("input.json"))
        LOG_INFO("input.json not found");

    std::vector<std::pair<double, double>> detectors = { {0,0} };
    std::vector<c2vector> spectrums(detectors.size());
    auto frequencies = Range(dff, ffmax, ffmin);
    for (auto& s : spectrums) s.resize(frequencies.size(), { 0,0 });

    // ------------- set up lasers --------------------
    LOG_INFO("set up lasers");
    std::vector<Laser> lasers;
    for (int i = 0; i < peak_I0_wcm2.size(); i++) {
        lasers.emplace_back(peak_I0_wcm2[i], waist_um[i], Lnm[i], g0[i]);
    }

    // -------------- set up gas jet ------------------
    LOG_INFO("set up gas jet");
    CylindricalGasJet gas_jet(gas_density_cm3, gas_sig_um, gas_radius_um, gas_length_um, gas_cells);
    gas_jet.SampleCylinder(lasers);

    // -------------- set up threadpool ---------------
    LOG_INFO("set up threadpool");
    ThreadPool::Startup(worker_thread);
    size_t numWorkers = ThreadPool::WorkerCount();
    size_t jobsPerWork = gas_jet.cells.size() / numWorkers;
    size_t remainingJobs = gas_jet.cells.size() % numWorkers;
    std::vector<std::future<void>> futures;
    std::vector<std::vector<c2vector>> interm(numWorkers);
    for (auto& i : interm) i.resize(detectors.size());
    for (auto& i : interm) for (auto& d : i) d.resize(frequencies.size(), {0,0}); 
    size_t jobs_done = 0;

    // ------------- begin calculation ---------------
    LOG_INFO("begin calculation");
    for (size_t iw = 0; iw < numWorkers; iw++) {
        size_t jobs = jobsPerWork + (iw < remainingJobs ? 1 : 0);
            
        futures.push_back(ThreadPool::PushTask([jobs, iw, jobs_done, &interm, &gas_jet, &detectors](size_t)
        {
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
            }
            sfa.pulses.insert(sfa.pulses.end(), pulses.begin(), pulses.end());
            sfa.SetupVectorization();
            sfa.SetupDTM();

            // --------------macroscopic propagation-----------


            for (int ic = 0; ic < jobs; ic++) {                                // for each cell
                auto& cell = gas_jet.cells[ic];
                auto rj = cell.pos;

                LOG_DEBUG("set sfa pulse");
                for (int i = 0; i < cell.intensity.size(); i++) {
                    pulses[i]->CEP = cell.phase[i];
                    pulses[i]->E0 = sqrt(cell.intensity[i]);
                }
                

                LOG_DEBUG("SetupFieldArrays");
                sfa.SetupFieldArrays();
                LOG_DEBUG("Execute2D");
                sfa.Execute2D();
                LOG_DEBUG("Spectrum");
                sfa.Spectrum();

                LOG_DEBUG("macroscopic propagation");
                for (int i = 0; i < interm[iw].size(); i++) {                    // for each detector
                    auto& d = detectors[i];
                    auto& spectrum = interm[iw][i];

                    double td = d.first;
                    double pd = d.second;

                    double arg = rj.x * cos(td) * sin(pd) +
                                rj.y * sin(td) * sin(pd) +
                                rj.z * (cos(td) - 1);

                    for (int j = 0; j < sfa.frequencies.size(); j++) {
                        double w = sfa.frequencies[j];
                        cvec2 radiation = std::exp(-1.i * (w / C) * arg) * sfa.hhg[j];

                        spectrum[j] += cell.density * radiation;
                    }
                }
            }
            //sfa.FreeVectorization();
        }
        ));
        jobs_done += jobs;
    }

    LOG_INFO("waiting for " + std::to_string(ThreadPool::WorkerCount()) + " threads");
    // wait for threads

    for (auto& f: futures)
        f.wait();

    LOG_INFO("combine partial sums");
    // combine partial sums
    for (size_t iw = 0; iw < numWorkers; iw++) {                            // for each worker
        for (int id = 0; id < interm[iw].size(); id++) {                    // for each detector
            for (int j = 0; j < frequencies.size(); j++) {                  // for each frequency
                spectrums[id][j] += interm[iw][id][j];
            }
        }
    }

    ThreadPool::Shutdown();

    Store(output_filename, frequencies, spectrums[0]);


#if defined(PROFILING)
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

    auto gas_jet = input["gas_jet"];
    // gas parameters
    gas_radius_um = gas_jet["radius_um"].get<double>();
    gas_length_um = gas_jet["length_um"].get<double>();
    gas_cells = gas_jet["cells"].get<size_t>();
    gas_sig_um = gas_jet["sigma_um"].get<double>();
    gas_density_cm3 = gas_jet["density_cm3"].get<double>();



    return true;
}