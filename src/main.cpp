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

size_t worker_thread = 16;
// laser parameters
double peak_I0_wcm2 = 1e14;
double waist_um = 30;
double Lnm = 800;
double g0 = 0;
double Ncyc = 10;
double w0 = LnmToAU / Lnm;

// gas parameters
double gas_radius_um = 500;
double gas_length_um = 6 * waist_um;
size_t gas_cells = 1;
double gas_sig_um = 800;
double gas_density_cm3 = 1;

double dt = 0.2;
double tmax = (2. * Pi / w0) * Ncyc;
double dp = 0.00001;
double pmin = 0, pmax = 3;
double dff = 0.01;
double ffmin = 0. * w0, ffmax = 40. * w0;


int main() {
    Log::set_logger_file("log.txt");
    LOG_INFO("main()");

    std::vector<std::pair<double, double>> detectors = { {0,0} };
    std::vector<c2vector> spectrums(detectors.size());
    auto frequencies = Range(dff, ffmax, ffmin);
    for (auto& s : spectrums) s.resize(frequencies.size(), { 0,0 });

    // ------------- set up lasers --------------------
    LOG_INFO("set up lasers");
    std::vector<Laser> lasers;
    lasers.emplace_back(peak_I0_wcm2, waist_um, Lnm, g0);

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
            
        futures.push_back(ThreadPool::PushTask([jobs, iw, jobs_done, &interm, &gas_jet, &detectors](size_t) {
            // ----------- set up sfa ------------------------
            LOG_DEBUG("set up sfa");
            SFA::SFA sfa;
            sfa.ts = Range(dt, tmax);
            sfa.ps = Range(dp, pmax, pmin);
            sfa.frequencies = Range(dff, ffmax, ffmin);
            sfa.Ip = 0.5;
            sfa.dtm = std::make_shared<SFA::ZeroRange>(sfa.Ip);
            std::vector<SFA::Sin2Pulse::Ptr_t> pulses = {
                std::make_shared<SFA::Sin2Pulse>(0., LnmToAU / Lnm, Ncyc, 0., dvec2{1., 0.})
            };
            sfa.pulses.insert(sfa.pulses.begin(), pulses.begin(), pulses.end());
        
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
            sfa.FreeVectorization();
        }));
    }

    LOG_INFO("waiting for " + std::to_string(ThreadPool::WorkerCount()) + " threads");
    // wait for threads
    exit(0);

    for (auto& f: futures)
        f.wait();

    exit(0);

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

    Store("out.dat", frequencies, spectrums[0]);


#if defined(PROFILING)
    Profile::Print();
    Profile::PrintTo("profile.txt");
#endif
}