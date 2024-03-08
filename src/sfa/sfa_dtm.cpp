#include "sfa.h"
#include "utility/logger.h"
#include <algorithm>

namespace SFA {
    cvec2 SFA::DTM(const dvec2& p) {
        static const double dp = ps[1] - ps[0];
        static const double invdp = 1. / dp;

        double l = length(p);
        double di = (l - ps[0]) * invdp;
        int i = int(di);
        i = std::max(0, std::min(i, int(ps.size()-2)));
        int j = i + 1;
        LOG_DEBUG(std::to_string(ps.size()));
        LOG_DEBUG(std::to_string(j));
        if (di == double(i))
            return dtm_el[i] * p;
            


        double f = (l - ps[i]) * invdp;
        complex dtm = dtm_el[j] * f + (1. - f) * dtm_el[i];

        return dtm * p;
    }



    void SFA::ComputeDtmFast(int itr, int iti, 
        __m256d& dtmTrRealx, __m256d& dtmTrRealy, __m256d& dtmTrImagx, __m256d& dtmTrImagy, 
        __m256d& dtmTiRealx, __m256d& dtmTiRealy, __m256d& dtmTiImagx, __m256d& dtmTiImagy) {
        
        BEGIN_ALIGNED(32) double delt[4] END_ALIGNED(32) = { 0 } ;
        BEGIN_ALIGNED(32) double Pxtr[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double Pytr[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double Pxti[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double Pyti[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double trf1[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double tif1[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double trf2[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double tif2[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmreali[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmrealj[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmimagi[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double dtmimagj[4] END_ALIGNED(32) = { 0 };

        for (int i = 0; i < 4; i++) {
            delt[i] = ts[iti + i] - ts[itr];
            dvec2 p = dvec2{ (intA[itr] - intA[iti + i]).x, 
                            (intA[itr] - intA[iti + i]).y } / delt[i];

            Pxtr[i] = p.x + A[itr].x;
            Pytr[i] = p.y + A[itr].y;
            Pxti[i] = p.x + A[iti + i].x;
            Pyti[i] = p.y + A[iti + i].y;
        }

        // dtm(tr)
        {
            // |P(tr)|
            auto Ptr = _mm256_sqrt_pd(_mm256_add_pd(
                _mm256_mul_pd(*(__m256d*)Pxtr, *(__m256d*)Pxtr),
                _mm256_mul_pd(*(__m256d*)Pytr, *(__m256d*)Pytr)
            ));
            
            // Np = (|P| - Ps0) / dp
            auto trNi = _mm256_mul_pd(_mm256_sub_pd(Ptr, vecPs0), vecInvdp);
            
            // dN = (|P(t)| - ps(Np)) / dp
            for (int i = 0; i < 4; i++) {
                trf1[i] = (Ptr[i] - ps[int(trNi[i])]) * invdp;
                trf2[i] = 1. - trf1[i];
            }
            for (int i = 0; i < 4; i++) {
                dtmreali[i] =  std::real(dtm_el[int(trNi[i])]);
                dtmrealj[i] =  std::real(dtm_el[int(trNi[i]) + 1]);
                dtmimagi[i] = -std::imag(dtm_el[int(trNi[i])]);
                dtmimagj[i] = -std::imag(dtm_el[int(trNi[i]) + 1]);
            }
            // lerp = f1 * j + f2 * i; f2 = 1-f1
            auto areal = _mm256_mul_pd(*(__m256d*)dtmimagj, *(__m256d*)trf1);
            auto aimag = _mm256_mul_pd(*(__m256d*)dtmrealj, *(__m256d*)trf1);
            auto breal = _mm256_mul_pd(*(__m256d*)dtmimagi, *(__m256d*)trf2);
            auto bimag = _mm256_mul_pd(*(__m256d*)dtmreali, *(__m256d*)trf2);
            auto imag = _mm256_add_pd(aimag, bimag);
            auto real = _mm256_add_pd(areal, breal);
            dtmTrImagx =  _mm256_mul_pd(imag, *(__m256d*)Pxtr);
            dtmTrImagy = _mm256_mul_pd(imag, *(__m256d*)Pytr);
            dtmTrRealx = _mm256_mul_pd(real, *(__m256d*)Pxtr);
            dtmTrRealy = _mm256_mul_pd(real, *(__m256d*)Pytr);
        }
        // dtm(ti)
        {
            // |P(ti)|
            auto Pti = _mm256_sqrt_pd(_mm256_add_pd(
                _mm256_mul_pd(*(__m256d*)Pxti, *(__m256d*)Pxti),
                _mm256_mul_pd(*(__m256d*)Pyti, *(__m256d*)Pyti)
            ));
            auto tiNi = _mm256_mul_pd(_mm256_sub_pd(Pti, vecPs0), vecInvdp);
            // dN = (|P(t)| - ps(Np)) / dp
            for (int i = 0; i < 4; i++) {
                tif1[i] = (Pti[i] - ps[int(tiNi[i])]) * invdp;
                tif2[i] = 1. - tif1[i];
            }
            for (int i = 0; i < 4; i++) {
                dtmreali[i] = std::real(dtm_el[int(tiNi[i])]);
                dtmrealj[i] = std::real(dtm_el[int(tiNi[i]) + 1]);
                dtmimagi[i] = std::imag(dtm_el[int(tiNi[i])]);
                dtmimagj[i] = std::imag(dtm_el[int(tiNi[i]) + 1]);
            }
            // lerp = f1 * j + f2 * i; f2 = 1-f1
            auto areal = _mm256_mul_pd(*(__m256d*)dtmimagj, *(__m256d*)tif1);
            auto aimag = _mm256_mul_pd(*(__m256d*)dtmrealj, *(__m256d*)tif1);
            auto breal = _mm256_mul_pd(*(__m256d*)dtmimagi, *(__m256d*)tif2);
            auto bimag = _mm256_mul_pd(*(__m256d*)dtmreali, *(__m256d*)tif2);
            auto imag = _mm256_add_pd(aimag, bimag);
            auto real = _mm256_add_pd(areal, breal);
            dtmTiImagx = _mm256_mul_pd(imag, *(__m256d*)Pxti);
            dtmTiImagy = _mm256_mul_pd(imag, *(__m256d*)Pyti);
            dtmTiRealx = _mm256_mul_pd(real, *(__m256d*)Pxti);
            dtmTiRealy = _mm256_mul_pd(real, *(__m256d*)Pyti);
        }
    }


    void SFA::ComputeDtm(int itr, int iti, cvec2& dtmTr, cvec2& dtmTi) {
        double deltat = ts[iti] - ts[itr];
        dvec2 p = dvec2{ (intA[itr] - intA[iti]).x, (intA[itr] - intA[iti]).y } / deltat;
        auto dtmtr = conj(DTM(p + A[itr]));
        auto dtmti = DTM(p + A[iti]);
    }

}