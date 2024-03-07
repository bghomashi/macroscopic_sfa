#include "maths/vec2.h"
#include <string>

extern bool sfa_only;
extern size_t worker_thread;
// laser parameters
extern dvector peak_I0_wcm2;
extern dvector waist_um;
extern dvector Lnm;
extern dvector g0;
extern dvector Ncyc;
extern dvector w0;
extern d2vector major_pol;
extern dvector cep;
extern dvector cycles_delay;

// gas parameters
extern double gas_radius_um;
extern double gas_length_um;
extern size_t gas_cells;
extern double gas_sig_um;
extern double gas_density_cm3;

extern double dt;
extern double tmax;
extern double dp;
extern double pmin;
extern double pmax;
extern double dff;
extern double ffmin;
extern double ffmax;

extern std::string output_filename;