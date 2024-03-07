#include "sfa.h"
#include "maths/fast_maths.h"
#include "utility/profiler.h"

#include <iostream>

namespace SFA {
    void SFA::FillIntermediateArrays() {
        LOG_DEBUG("FillIntermediateArrays");
#if defined(PROFILING)
        Profile::Push("SFA::FillIntermediateArrays");
#endif
        const double dp = ps[1] - ps[0];
        const double invdp = 1. / dp;

        BEGIN_ALIGNED(32) double delt[4] END_ALIGNED(32) = { 0 } ;
        BEGIN_ALIGNED(32) double PP[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double DintAA[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double trL[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double tiL[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double trDi[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double tiDi[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double trPx[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double trPy[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double trf[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double tiPx[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double tiPy[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double tif[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmtrimagi[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmtrimagj[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmtrreali[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmtrrealj[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmtiimagi[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmtiimagj[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmtireali[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmtirealj[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp_rx[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp_ry[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp_ix[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp_iy[4] END_ALIGNED(32) = { 0 };

        const __m256d vecIp = _mm256_set1_pd(Ip);
        const __m256d vecPs0 = _mm256_set1_pd(ps[0]);
        const __m256d vecInvdp = _mm256_set1_pd(invdp);


std::cout << "FillIntermediateArrays" << std::endl;

        LOG_DEBUG("FillIntermediateArrays " + std::to_string(__LINE__));
        Log::flush();
        for (int itr = 2; itr < ts.size(); itr++) {             // tr
if (itr >= 21090) std::cout << " " << std::to_string(itr) << std::endl;
            int idx = index(itr, 0);
            int l = (itr - 2) / 4;
            int iti = 1;

            for (iti = 1; iti < l * 4; iti += 4) {             // ti

                for (int i = 0; i < 4; i++) {
if (itr >= 21090) std::cout << "iti+i - itr " << std::to_string(iti + i) << "-" << std::to_string(itr) << std::endl;
if (itr >= 21090) std::cout << "iti+i - itr " << std::to_string(ts[iti + i]) << "-" << std::to_string(ts[itr]) << std::endl;

                    delt[i] = ts[iti + i] - ts[itr];
                    dvec2 p = dvec2{ (intA[itr] - intA[iti + i]).x, (intA[itr] - intA[iti + i]).y } / delt[i];
                    PP[i] = p.x * p.x + p.y * p.y;
                    DintAA[i] = intAA[itr] - intAA[iti + i];

                    trPx[i] = p.x + A[itr].x;
                    trPy[i] = p.y + A[itr].y;
                    tiPx[i] = p.x + A[iti].x;
                    tiPy[i] = p.y + A[iti].y;
                    
if (itr >= 21090) std::cout << "delt " << delt[i] << std::endl;
if (itr >= 21090) std::cout << "p " << p.x << " " << p.y << std::endl;
if (itr >= 21090) std::cout << "A[itr] " << A[itr].x << " " << A[itr].y << std::endl;
if (itr >= 21090) std::cout << "A[iti] " << A[iti].x << " " << A[iti].y << std::endl;
                }

if (itr >= 21090) std::cout << " " << std::to_string(__LINE__) << std::endl;
if (itr >= 21090) std::cout << "tr " << trPx[0] << " " << trPy[0]<< std::endl;
if (itr >= 21090) std::cout << "ti " << tiPx[0] << " " << tiPy[0]<< std::endl;
                _mm256_store_pd(trL, _mm256_sqrt_pd(_mm256_add_pd(
                    _mm256_mul_pd(*(__m256d*)trPx, *(__m256d*)trPx),
                    _mm256_mul_pd(*(__m256d*)trPy, *(__m256d*)trPy)
                )));
                _mm256_store_pd(tiL, _mm256_sqrt_pd(_mm256_add_pd(
                    _mm256_mul_pd(*(__m256d*)tiPx, *(__m256d*)tiPx),
                    _mm256_mul_pd(*(__m256d*)tiPy, *(__m256d*)tiPy)
                )));

                _mm256_store_pd(trDi, _mm256_mul_pd(_mm256_sub_pd(*(__m256d*)trL, vecPs0), vecInvdp));
                _mm256_store_pd(tiDi, _mm256_mul_pd(_mm256_sub_pd(*(__m256d*)tiL, vecPs0), vecInvdp));

if (itr >= 21090) std::cout << " " << std::to_string(__LINE__) << std::endl;
                for (int i = 0; i < 4; i++) {
                    trf[i] = (trL[i] - ps[int(trDi[i])]) * invdp;
                    tif[i] = (tiL[i] - ps[int(tiDi[i])]) * invdp;

if (itr >= 21090) std::cout << " " << std::to_string(__LINE__) << std::endl;
                    dtmtrreali[i] = std::real(dtm_el[int(trDi[i])]);
                    dtmtrrealj[i] = std::real(dtm_el[int(trDi[i]) + 1]);
                    dtmtrimagi[i] = -std::imag(dtm_el[int(trDi[i])]);
                    dtmtrimagj[i] = -std::imag(dtm_el[int(trDi[i]) + 1]);

if (itr >= 21090) std::cout << " " << std::to_string(__LINE__) << std::endl;

                    dtmtireali[i] = std::real(dtm_el[int(tiDi[i])]);
                    dtmtirealj[i] = std::real(dtm_el[int(tiDi[i]) + 1]);
                    dtmtiimagi[i] = std::imag(dtm_el[int(tiDi[i])]);
                    dtmtiimagj[i] = std::imag(dtm_el[int(tiDi[i]) + 1]);
                }

if (itr >= 21090) std::cout << " " << std::to_string(__LINE__) << std::endl;
                // dtmtr
                {
                    auto f2 = _mm256_sub_pd(fast::one, *(__m256d*)trf);
                    auto areal = _mm256_mul_pd(*(__m256d*)dtmtrimagj, *(__m256d*)trf);
                    auto aimag = _mm256_mul_pd(*(__m256d*)dtmtrrealj, *(__m256d*)trf);
                    auto breal = _mm256_mul_pd(*(__m256d*)dtmtrimagi, f2);
                    auto bimag = _mm256_mul_pd(*(__m256d*)dtmtrreali, f2);
                    auto imag = _mm256_add_pd(aimag, bimag);
                    auto real = _mm256_add_pd(areal, breal);
                    _mm256_store_pd(temp_ix, _mm256_mul_pd(imag, *(__m256d*)trPx));
                    _mm256_store_pd(temp_iy, _mm256_mul_pd(imag, *(__m256d*)trPy));
                    _mm256_store_pd(temp_rx, _mm256_mul_pd(real, *(__m256d*)trPx));
                    _mm256_store_pd(temp_ry, _mm256_mul_pd(real, *(__m256d*)trPy));

                    for (int i = 0; i < 4; i++) {
                        dtm_itr_real_x[idx + iti + i] = temp_rx[i];
                        dtm_itr_real_y[idx + iti + i] = temp_ry[i];
                        dtm_itr_imag_x[idx + iti + i] = temp_ix[i];
                        dtm_itr_imag_y[idx + iti + i] = temp_ix[i];
                    }
                }
if (itr >= 21090) std::cout << " " << std::to_string(__LINE__) << std::endl;
                {
                    auto f2 = _mm256_sub_pd(fast::one, *(__m256d*)tif);
                    auto areal = _mm256_mul_pd(*(__m256d*)dtmtiimagj, *(__m256d*)tif);
                    auto aimag = _mm256_mul_pd(*(__m256d*)dtmtirealj, *(__m256d*)tif);
                    auto breal = _mm256_mul_pd(*(__m256d*)dtmtiimagi, f2);
                    auto bimag = _mm256_mul_pd(*(__m256d*)dtmtireali, f2);
                    auto imag = _mm256_add_pd(aimag, bimag);
                    auto real = _mm256_add_pd(areal, breal);
                    _mm256_store_pd(temp_ix, _mm256_mul_pd(imag, *(__m256d*)tiPx));
                    _mm256_store_pd(temp_iy, _mm256_mul_pd(imag, *(__m256d*)tiPy));
                    _mm256_store_pd(temp_rx, _mm256_mul_pd(real, *(__m256d*)tiPx));
                    _mm256_store_pd(temp_ry, _mm256_mul_pd(real, *(__m256d*)tiPy));
                    for (int i = 0; i < 4; i++) {
                        dtm_iti_real_x[idx + iti + i] = temp_rx[i];
                        dtm_iti_real_y[idx + iti + i] = temp_ry[i];
                        dtm_iti_imag_x[idx + iti + i] = temp_ix[i];
                        dtm_iti_imag_y[idx + iti + i] = temp_ix[i];
                    }
                }
                
if (itr >= 21090) std::cout << " " << std::to_string(__LINE__) << std::endl;
                // ------------------ S0 
                auto H0 = _mm256_add_pd(_mm256_mul_pd(fast::negHalf, *(__m256d*)PP), vecIp);  // -0.5*(P^2) + Ip
                auto S0 = _mm256_mul_pd(H0, *(__m256d*)delt);
                auto SAA = _mm256_mul_pd(fast::negHalf, *(__m256d*)DintAA);
                auto action = _mm256_add_pd(S0, SAA);
                _mm256_store_pd(temp_rx, _mm256_sub_pd(action, fast::Pi));
                for (int i = 0; i < 4; i++) 
                        S0_saddle[idx + iti + i] = temp_rx[i];

if (itr >= 21090) std::cout << " " << std::to_string(__LINE__) << std::endl;
            }


        LOG_DEBUG("FillIntermediateArrays " + std::to_string(__LINE__));
if (itr >= 21090) std::cout << " " << std::to_string(itr) << std::endl;

        LOG_DEBUG(std::to_string(sizeof(S0_saddle)));
            for (; iti < itr; iti++) {
        LOG_DEBUG(std::to_string(iti));
                double deltat = ts[iti] - ts[itr];
                dvec2 p = dvec2{ (intA[itr] - intA[iti]).x, (intA[itr] - intA[iti]).y } / deltat;
                S0_saddle[idx + iti] = (-0.5 * (p.x * p.x + p.y * p.y) + Ip) * deltat - 0.5 * (intAA[itr] - intAA[iti]) - Pi;

                auto dtmtr = conj(DTM(p + A[itr]));
                auto dtmti = DTM(p + A[iti]);

                dtm_itr_real_x[idx + iti] = std::real(dtmtr.x);
                dtm_itr_imag_x[idx + iti] = std::imag(dtmtr.x);
                dtm_itr_real_y[idx + iti] = std::real(dtmtr.y);
                dtm_itr_imag_y[idx + iti] = std::imag(dtmtr.y);

                dtm_iti_real_x[idx + iti] = std::real(dtmti.x);
                dtm_iti_imag_x[idx + iti] = std::imag(dtmti.x);
                dtm_iti_real_y[idx + iti] = std::real(dtmti.y);
                dtm_iti_imag_y[idx + iti] = std::imag(dtmti.y);
            }

        }
        

        LOG_DEBUG("FillIntermediateArrays " + std::to_string(__LINE__));

        for (int it = 0; it < ts.size(); it++) {             // tr = ti
            int idx = index(it, it);
            auto dtmtr = DTM({ A[it].x, A[it].y });

            S0_saddle[idx] = 0;
            dtm_iti_real_x[idx] = std::real(dtmtr.x);
            dtm_iti_imag_x[idx] = std::imag(dtmtr.x);
            dtm_iti_real_y[idx] = std::real(dtmtr.y);
            dtm_iti_imag_y[idx] = std::imag(dtmtr.y);

            // conj
            dtm_itr_real_x[idx] = std::real(dtmtr.x);
            dtm_itr_imag_x[idx] = -std::imag(dtmtr.x);
            dtm_itr_real_y[idx] = std::real(dtmtr.y);
            dtm_itr_imag_y[idx] = -std::imag(dtmtr.y);
        }
#if defined(PROFILING)
        Profile::Pop("SFA::FillIntermediateArrays");
#endif
    }

}