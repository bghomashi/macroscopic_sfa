#include "sfa.h"
#include "maths/fast_maths.h"
#include "utility/profiler.h"

namespace SFA {
	void SFA::Execute2D() {
        int NT = int(ts.size());
        dipole.resize(NT, { 0,0 });

        FillIntermediateArrays();


        LOG_DEBUG("Done");
        double deltat, action;
        cvec2 dtmItr, dtmIti;
        const double dt = (ts[1] - ts[0]);                  // assuming linear spacing

        // temporary storage for vector calculations
        BEGIN_ALIGNED(32) double delt[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double EEx[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double EEy[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp_rx[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp_ry[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp_ix[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp_iy[4] END_ALIGNED(32) = { 0 };

#if defined(PROFILING)
        Profile::Push("SFA::Execute2D");
#endif

        LOG_DEBUG("loop");
        for (int itr = 2; itr < NT; itr++) {
            int l = (itr - 2) / 4;
            int r = (itr - 2) % 4;
            int iti;
            int idx = index(itr, 0);
            for (iti = 1; iti < l * 4; iti += 4) {
                for (int i = 0; i < 4; i++) {
                    delt[i] = ts[itr] - ts[iti + i];
                    EEx[i] = E[iti + i].x;
                    EEy[i] = E[iti + i].y;
                }

                __m256d action  = *(__m256d*)(&S0_saddle[idx + iti]);
                __m256d Ex      = *(__m256d*)EEx,
                        Ey      = *(__m256d*)EEy;
                __m256d dtmTiRealx = _mm256_loadu_pd(&dtm_iti_real_x[idx + iti]),
                        dtmTiImagx = _mm256_loadu_pd(&dtm_iti_imag_x[idx + iti]),
                        dtmTiRealy = _mm256_loadu_pd(&dtm_iti_real_y[idx + iti]),
                        dtmTiImagy = _mm256_loadu_pd(&dtm_iti_imag_y[idx + iti]);
                __m256d dtmTrRealx = _mm256_loadu_pd(&dtm_itr_real_x[idx + iti]),
                        dtmTrImagx = _mm256_loadu_pd(&dtm_itr_imag_x[idx + iti]),
                        dtmTrRealy = _mm256_loadu_pd(&dtm_itr_real_y[idx + iti]),
                        dtmTrImagy = _mm256_loadu_pd(&dtm_itr_imag_y[idx + iti]);

                // ----------------- compute coeff ------------------
                __m256d coeff_real = _mm256_div_pd(fast::sqrtPi2, _mm256_sqrt_pd(*(__m256d*)delt));        // sqrt(2pi)/sqrt(tr-ti)
                __m256d coeff_imag = coeff_real;

                // ----------------- action part --------------------
                __m256d arg = fast::reduceArg(action);
                __m256d exp_real = fast::_mm_cos5(arg);        // real part of exp(-i...)
                __m256d exp_imag = fast::_mm_sin5(arg);        // imag part of exp(-i...)

                coeff_real = _mm256_mul_pd(coeff_real, exp_real);
                coeff_imag = _mm256_mul_pd(coeff_imag, exp_imag);
                // ---------------- E.d -----------------------------
                __m256d ED_real = _mm256_add_pd(_mm256_mul_pd(Ex, dtmTiRealx), _mm256_mul_pd(Ey, dtmTiRealy));       // real(E.D)
                __m256d ED_imag = _mm256_add_pd(_mm256_mul_pd(Ex, dtmTiImagx), _mm256_mul_pd(Ey, dtmTiImagy));       // imag(E.D)

                // ---------------- complex mult to combine ---------
                __m256d factor_real = _mm256_sub_pd(_mm256_mul_pd(coeff_real, ED_real), _mm256_mul_pd(coeff_imag, ED_imag));
                __m256d factor_imag = _mm256_add_pd(_mm256_mul_pd(coeff_imag, ED_real), _mm256_mul_pd(coeff_real, ED_imag));
                // --------------------------------------------------
                _mm256_store_pd(temp_rx,
                    _mm256_sub_pd(
                        _mm256_mul_pd(factor_real, dtmTrRealx),
                        _mm256_mul_pd(factor_imag, dtmTrImagx)
                    )
                );
                _mm256_store_pd(temp_ix,
                    _mm256_add_pd(
                        _mm256_mul_pd(factor_imag, dtmTrRealx),
                        _mm256_mul_pd(factor_real, dtmTrImagx)
                    )
                );
                _mm256_store_pd(temp_ry,
                    _mm256_sub_pd(
                        _mm256_mul_pd(factor_real, dtmTrRealy),
                        _mm256_mul_pd(factor_imag, dtmTrImagy)
                    )
                );
                _mm256_store_pd(temp_iy,
                    _mm256_add_pd(
                        _mm256_mul_pd(factor_imag, dtmTrRealy),
                        _mm256_mul_pd(factor_real, dtmTrImagy)
                    )
                );
                for (int i = 0; i < 4; i++) {
                    realx[idx + iti + i] = temp_rx[i];
                    realy[idx + iti + i] = temp_ry[i];
                    imagx[idx + iti + i] = temp_ix[i];
                    imagy[idx + iti + i] = temp_iy[i];
                }
            }
            for (; iti < itr; iti++) {
                deltat = ts[itr] - ts[iti];

                action = S0_saddle[idx+iti];
                dtmItr = { {dtm_itr_real_x[idx + iti], dtm_itr_imag_x[idx + iti]}, {dtm_itr_real_y[idx + iti], dtm_itr_imag_y[idx + iti]} };
                dtmIti = { {dtm_iti_real_x[idx + iti], dtm_iti_imag_x[idx + iti]}, {dtm_iti_real_y[idx + iti], dtm_iti_imag_y[idx + iti]} };

                auto out = dtmItr * (sqrtPi2 / sqrt(deltat) *
                    std::exp(-1.i * (action + Pi)) *
                    (E[iti].x * dtmIti.x + E[iti].y * dtmIti.y));

                realx[idx] = std::real(out.x);
                imagx[idx] = std::imag(out.x);
                realy[idx] = std::real(out.y);
                imagy[idx] = std::imag(out.y);
            }
        }

        dip[0] = cvec2{ 0,0 };
        dip[1] = cvec2{ 0,0 };

        // reduction - sum across ti
        for (int itr = 2; itr < NT; itr++) {
            int idx = index(itr, 0);
            __m256d ptrX = _mm256_setzero_pd();

            //int iti;
            int l = (itr - 2) / 4;
            int r = (itr - 2) % 4;
            // for (iti = 1; iti < l * 4; iti += 4) {
            cvec2 test = { 0,0 };
            for (int iti = 1; iti < itr; iti++) {
                __m256d ptrXX = _mm256_setr_pd(realx[idx + iti], imagx[idx + iti], realy[idx + iti], imagy[idx + iti]);
                ptrX = _mm256_add_pd(ptrX, ptrXX);
            }
            _mm256_store_pd(temp_rx, ptrX);
            dip[itr].x = complex(temp_rx[0], temp_rx[1]);
            dip[itr].y = complex(temp_rx[2], temp_rx[3]);
        }

        for (int itr = 0; itr < NT; itr++) {
            dip[itr] *= dt;
            dip[itr] += conj(dip[itr]);
            dipole[itr] = { dip[itr].x, dip[itr].y };
        }
#if defined(PROFILING)
        Profile::Pop("SFA::Execute2D");
#endif
	}

}