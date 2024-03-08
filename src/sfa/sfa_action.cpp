#include "sfa.h"

namespace SFA {
        
    void SFA::ComputeActionFast(int itr, int iti, double* out) {
        BEGIN_ALIGNED(32) double delt[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double PP[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double DintAA[4] END_ALIGNED(32) = { 0 };
        BEGIN_ALIGNED(32) double temp[4] END_ALIGNED(32) = { 0 };
        
        for (int i = 0; i < 4; i++) {
            delt[i] = ts[iti + i] - ts[itr];
            dvec2 p = dvec2{ (intA[itr] - intA[iti + i]).x, 
                            (intA[itr] - intA[iti + i]).y } / delt[i];
            PP[i] = p.x * p.x + p.y * p.y;
            DintAA[i] = intAA[itr] - intAA[iti + i];

        }
        auto H0 = _mm256_add_pd(_mm256_mul_pd(fast::negHalf, *(__m256d*)PP), vecIp);  // -0.5*(P^2) + Ip
        auto S0 = _mm256_mul_pd(H0, *(__m256d*)delt);
        auto SAA = _mm256_mul_pd(fast::negHalf, *(__m256d*)DintAA);
        auto action = _mm256_add_pd(S0, SAA);
        _mm256_store_pd(out, _mm256_sub_pd(action, fast::Pi));
    }
    double SFA::ComputeAction(int itr, int iti) {
        double deltat = ts[iti] - ts[itr];
        dvec2 p = dvec2{ (intA[itr] - intA[iti]).x, (intA[itr] - intA[iti]).y } / deltat;
        return (-0.5 * (p.x * p.x + p.y * p.y) + Ip) * deltat - 0.5 * (intAA[itr] - intAA[iti]) - Pi;

    }

}