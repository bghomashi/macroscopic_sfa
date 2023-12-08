#pragma once

#include "maths/maths.h"
#include <immintrin.h>

#if defined(_MSC_VER)
#define BEGIN_ALIGNED(x) __declspec(align(x)) 
#define END_ALIGNED(x)
#elif defined(__GNUC__)
#define BEGIN_ALIGNED(x) 
#define END_ALIGNED(x)  __attribute__((aligned(x)))
#endif

namespace fast {
	constexpr double S0 = 1.;
	constexpr double S1 = -1. / 2. / 3.;
	constexpr double S2 = 1. / 2. / 3. / 4. / 5.;
	constexpr double S3 = -1. / 2. / 3. / 4. / 5. / 6. / 7.;
	constexpr double S4 = 1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9.;
	constexpr double S5 = -1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11.;
	constexpr double S6 = 1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13.;
	constexpr double S7 = -1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13. / 14. / 15.;
	constexpr double S8 = 1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13. / 14. / 15. / 16. / 17.;
	constexpr double S9 = -1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13. / 14. / 15. / 16. / 17. / 18. / 19.;
	constexpr double S10 = 1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13. / 14. / 15. / 16. / 17. / 18. / 19. / 20. / 21.;


	constexpr double C0 = 1.;
	constexpr double C1 = -1. / 2.;
	constexpr double C2 = 1. / 2. / 3. / 4.;
	constexpr double C3 = -1. / 2. / 3. / 4. / 5. / 6.;
	constexpr double C4 = 1. / 2. / 3. / 4. / 5. / 6. / 7. / 8.;
	constexpr double C5 = -1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10.;
	constexpr double C6 = 1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12.;
	constexpr double C7 = -1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13. / 14.;
	constexpr double C8 = 1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13. / 14. / 15. / 16.;
	constexpr double C9 = -1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13. / 14. / 15. / 16. / 17. / 18.;
	constexpr double C10 = 1. / 2. / 3. / 4. / 5. / 6. / 7. / 8. / 9. / 10. / 11. / 12. / 13. / 14. / 15. / 16. / 17. / 18. / 19. / 20.;

	const __m256d one = _mm256_set1_pd(1.);
	const __m256d Half = _mm256_set1_pd(0.5);
	const __m256d negHalf = _mm256_set1_pd(-0.5);
	const __m256d Pi = _mm256_set1_pd(::Pi);
	const __m256d Pi2 = _mm256_set1_pd(::Pi2);
	const __m256d inv2Pi = _mm256_set1_pd(::inv2Pi);
	const __m256d invPi_2 = _mm256_set1_pd(::invPi_2);
	const __m256d Pi_2 = _mm256_set1_pd(::Pi_2);
	const __m256d sqrtPi2 = _mm256_set1_pd(::sqrtPi2);

	const __m256d s0 = _mm256_set1_pd(S0);
	const __m256d s1 = _mm256_set1_pd(S1);
	const __m256d s2 = _mm256_set1_pd(S2);
	const __m256d s3 = _mm256_set1_pd(S3);
	const __m256d s4 = _mm256_set1_pd(S4);
	const __m256d s5 = _mm256_set1_pd(S5);
	const __m256d s6 = _mm256_set1_pd(S6);
	const __m256d s7 = _mm256_set1_pd(S7);
	const __m256d s8 = _mm256_set1_pd(S8);
	const __m256d s9 = _mm256_set1_pd(S9);
	const __m256d s10 = _mm256_set1_pd(S10);

	const __m256d c0 = _mm256_set1_pd(C0);
	const __m256d c1 = _mm256_set1_pd(C1);
	const __m256d c2 = _mm256_set1_pd(C2);
	const __m256d c3 = _mm256_set1_pd(C3);
	const __m256d c4 = _mm256_set1_pd(C4);
	const __m256d c5 = _mm256_set1_pd(C5);
	const __m256d c6 = _mm256_set1_pd(C6);
	const __m256d c7 = _mm256_set1_pd(C7);
	const __m256d c8 = _mm256_set1_pd(C8);
	const __m256d c9 = _mm256_set1_pd(C9);
	const __m256d c10 = _mm256_set1_pd(C10);

	inline __m256d _mm_sin10(const __m256d x) {
		__m256d xx = _mm256_mul_pd(x,x);
		__m256d s;

		s = _mm256_fmadd_pd(s10, xx, s9);
		s = _mm256_fmadd_pd(s, xx, s8);
		s = _mm256_fmadd_pd(s, xx, s7);
		s = _mm256_fmadd_pd(s, xx, s6);
		s = _mm256_fmadd_pd(s, xx, s5);
		s = _mm256_fmadd_pd(s, xx, s4);
		s = _mm256_fmadd_pd(s, xx, s3);
		s = _mm256_fmadd_pd(s, xx, s2);
		s = _mm256_fmadd_pd(s, xx, s1);
		s = _mm256_fmadd_pd(s, xx, s0);
		s = _mm256_mul_pd(s, x);

		return s;
	}
	inline __m256d _mm_sin5(const __m256d x) {
		__m256d xx = _mm256_mul_pd(x, x);
		__m256d s;

		s = _mm256_fmadd_pd(s5, xx, s4);
		s = _mm256_fmadd_pd(s, xx, s3);
		s = _mm256_fmadd_pd(s, xx, s2);
		s = _mm256_fmadd_pd(s, xx, s1);
		s = _mm256_fmadd_pd(s, xx, s0);
		s = _mm256_mul_pd(s, x);

		return s;
	}

	inline __m256d _mm_cos10(const __m256d x) {
		__m256d xx = _mm256_mul_pd(x, x);
		__m256d c;

		c = _mm256_fmadd_pd(c10, xx, c9);
		c = _mm256_fmadd_pd(c, xx, c8);
		c = _mm256_fmadd_pd(c, xx, c7);
		c = _mm256_fmadd_pd(c, xx, c6);
		c = _mm256_fmadd_pd(c, xx, c5);
		c = _mm256_fmadd_pd(c, xx, c4);
		c = _mm256_fmadd_pd(c, xx, c3);
		c = _mm256_fmadd_pd(c, xx, c2);
		c = _mm256_fmadd_pd(c, xx, c1);
		c = _mm256_fmadd_pd(c, xx, c0);

		return c;
	}
	inline __m256d _mm_cos5(const __m256d x) {
		__m256d xx = _mm256_mul_pd(x, x);
		__m256d c;

		c = _mm256_fmadd_pd(c5, xx, c4);
		c = _mm256_fmadd_pd(c, xx, c3);
		c = _mm256_fmadd_pd(c, xx, c2);
		c = _mm256_fmadd_pd(c, xx, c1);
		c = _mm256_fmadd_pd(c, xx, c0);

		return c;
	}

	//inline __m256d quadrant(const __m256d x) {
	//	// x/(pi/4) = [-4, -3]   -> -sin(x+3pi/4)
	//	// x/(pi/4) = [-3, -1]   -> -cos(x+pi/2)
	//	// x/(pi/4) = [-1,  1] 	 ->  sin(x)
	//	// x/(pi/4) = [ 1,  3] 	 ->  cos(x-pi/2)
	//	// x/(pi/4) = [ 3,  4] 	 -> -sin(x-3pi/4)
	//	__m256d arg1 = _mm256_add_pd(x, Pi);									// arg1 = x + pi
	//	__m256d arg = _mm256_mul_pd(arg1, inv2Pi);								// arg = arg1 / 2pi
	//	arg = _mm256_round_pd(arg, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);	// arg = floor(arg)
	//	arg = _mm256_mul_pd(arg, Pi2);											// arg = 2pi * arg
	//	arg = _mm256_sub_pd(arg1, arg);											// arg = arg1 - arg;  [0, 2pi)
	//	arg = _mm256_sub_pd(arg, Pi);											// arg = arg - pi; [-pi,pi)

	//	return arg;
	//}

	inline __m256d reduceArg(const __m256d x) {
		__m256d arg1 = _mm256_add_pd(x, Pi);									// arg1 = x + pi
		__m256d arg = _mm256_mul_pd(arg1, inv2Pi);								// arg = arg1 / 2pi
		arg = _mm256_round_pd(arg, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);	// arg = floor(arg)
		arg = _mm256_mul_pd(arg, Pi2);											// arg = 2pi * arg
		arg = _mm256_sub_pd(arg1, arg);											// arg = arg1 - arg;  [0, 2pi)
		arg = _mm256_sub_pd(arg, Pi);											// arg = arg - pi; [-pi,pi)
		
		return arg;
	}


	
	inline void sin(const double* x, double* out) {
		const __m256d* ptr = (const __m256d*)x;
		__m256d arg = reduceArg(*ptr);
		__m256d sines = _mm_sin10(arg);
		_mm256_store_pd(out, sines);
	}

	inline void cos(const double* x, double* out) {
		const __m256d* ptr = (const __m256d*)x;
		__m256d arg = reduceArg(*ptr);
		__m256d cosines = _mm_cos10(arg);
		_mm256_store_pd(out, cosines);
	}


	inline void* aligned_malloc(size_t size, size_t align) {
		void* result;
#ifdef _MSC_VER 
		result = _aligned_malloc(size, align);
#else 
		if (posix_memalign(&result, align, size)) result = 0;
#endif
		return result;
	}

	inline void aligned_free(void* ptr) {
#ifdef _MSC_VER 
		_aligned_free(ptr);
#else 
		free(ptr);
#endif
	}
}