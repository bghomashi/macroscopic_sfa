#include "maths/maths.h"
#include "maths/vec2.h"
#include "sfa/sfa.h"
#include "laser/sin2pulse.h"
#include "gas_jet/gas_jet.h"
#include "laser/laser.h"
#include "dtms/zero_range_dtm.h"
#include "utility/save_data.h"
#include "utility/profiler.h"

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
size_t gas_cells = 1000;
double gas_sig_um = 800;
double gas_density_cm3 = 1;

double dt = 0.2;
double tmax = (2. * Pi / w0) * Ncyc;
double dp = 0.00001;
double pmin = 0, pmax = 3;
double dff = 0.01;
double ffmin = 0. * w0, ffmax = 40. * w0;


int main() {
    std::vector<std::pair<double, double>> detectors = { {0,0} };
    std::vector<c2vector> spectrums(detectors.size());

    // ----------- set up sfa ------------------------
    SFA::SFA sfa;
    sfa.ts = Range(dt, tmax);
    sfa.ps = Range(dp, pmax, pmin);
    sfa.frequencies = Range(dff, ffmax, ffmin);
    sfa.Ip = 0.5;
    sfa.dtm = std::make_shared<SFA::ZeroRange>(sfa.Ip);
    auto& frequencies = sfa.frequencies;                    // for convenience
    for (auto& s : spectrums)
        s.resize(frequencies.size(), { 0,0 });

    sfa.SetupVectorization();
    sfa.SetupDTM();

    // ------------- set up lasers --------------------
    std::vector<Laser> lasers;
    std::vector<SFA::Sin2Pulse::Ptr_t> pulses;

    lasers.emplace_back(peak_I0_wcm2, waist_um, Lnm, g0);
    pulses.push_back(std::make_shared<SFA::Sin2Pulse>(0., LnmToAU / Lnm, Ncyc, 0., dvec2{1., 0.}));
    sfa.pulses.push_back(pulses.back());                    // add to sfa

    // -------------- set up gas jet ------------------
    CylindricalGasJet gas_jet(gas_density_cm3, gas_sig_um, gas_radius_um, gas_length_um, gas_cells);
    gas_jet.SampleCylinder(lasers);

    for (int ic = 0; ic < gas_jet.cells.size(); ic++) {                                // for each cell
        auto& cell = gas_jet.cells[ic];
        auto rj = cell.pos;
        for (int i = 0; i < lasers.size(); i++) {
            pulses[i]->CEP = cell.phase[i];
            pulses[i]->E0 = sqrt(cell.intensity[i]);
        }

        sfa.SetupFieldArrays();
        sfa.Execute2D();
        sfa.Spectrum();

        for (int i = 0; i < detectors.size(); i++) {                    // for each detector
            auto& d = detectors[i];
            auto& spectrum = spectrums[i];

            double td = d.first;
            double pd = d.second;

            double arg = rj.x * cos(td) * sin(pd) +
                         rj.y * sin(td) * sin(pd) +
                         rj.z * (cos(td) - 1);

            for (int j = 0; j < frequencies.size(); j++) {
                double w = frequencies[j];
                cvec2 radiation = std::exp(-1.i * (w / C) * arg) * sfa.hhg[j];

                spectrum[j] += cell.density * radiation;
            }
        }
    }

    sfa.FreeVectorization();

    Store("out.dat", frequencies, spectrums[0]);


#if defined(PROFILING)
    Profile::Print();
#endif
}