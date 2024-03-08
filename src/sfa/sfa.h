#pragma once

#include <memory>
#include "maths/fast_maths.h"
#include "maths/maths.h"
#include "maths/vec2.h"

namespace SFA {
    class Pulse {
    public:
        typedef std::shared_ptr<Pulse> Ptr_t;

        virtual dvec2 A(double t) const = 0;
        virtual dvec2 E(double t) const = 0;
        virtual dvec2 intA(double t) const = 0;
        virtual double intAA(double t) const = 0;
    };
    class DTM {
    public:
        typedef std::shared_ptr<DTM> Ptr_t;
        virtual cvec2 value(const dvec2& p) const = 0;
    };


    class SFA {
        d2vector A, E, intA;
        dvector intAA;
        cvector dtm_el;

        // 32 byte aligned arrays for vectorization
        cvec2* dip;
        double* realx, * imagx, * realy, * imagy;
        double* dtm_iti_real_x, * dtm_iti_imag_x, * dtm_iti_real_y, * dtm_iti_imag_y;
        double* dtm_itr_real_x, * dtm_itr_imag_x, * dtm_itr_real_y, * dtm_itr_imag_y;
        double* S0_saddle;

        double dp, invdp;
        __m256d vecIp, vecPs0, vecInvdp;

        void ComputeDtmFast(int itr, int iti, 
            __m256d& dtmTrRealx, __m256d& dtmTrRealy, __m256d& dtmTrImagx, __m256d& dtmTrImagy, 
            __m256d& dtmTiRealx, __m256d& dtmTiRealy, __m256d& dtmTiImagx, __m256d& dtmTiImagy);
        void ComputeDtm(int itr, int iti, cvec2& dtmTr, cvec2& dtmTi);
        void ComputeActionFast(int itr, int iti, double* out);
        double ComputeAction(int itr, int iti);
        void FillIntermediateArrays();
        cvec2 DTM(const dvec2& p);

        static inline int index(int row, int col) {
            return (row + 1) * row / 2 + col;
        }
    public:
        // inputs
        std::vector<Pulse::Ptr_t> pulses;
        DTM::Ptr_t dtm;
        dvector ts, ps;
        dvector frequencies;
        double Ip;

        // outputs
        c2vector dipole;
        c2vector hhg;

        void Execute2D();   //
        void Spectrum();    // fourier transform dipole

        void SetupDTM();
        void SetupFieldArrays();

        void SetupVectorization();      // allocate aligned memory
        void FreeVectorization();       // free aligned memory
    };
}
