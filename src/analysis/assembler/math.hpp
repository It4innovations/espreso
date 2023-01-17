
#ifndef SRC_PHYSICS_ASSEMBLER_MATH_HPP_
#define SRC_PHYSICS_ASSEMBLER_MATH_HPP_

#include "math/simd/simd.h"

namespace espreso {

// |cos , -sin| |c[0] , c[2]| | cos , sin|
// |sin ,  cos| |c[3] , c[1]| |-sin , cos|
// |cos * c[0] - sin * c[3] , cos * c[2] - sin * c[1]|  | cos , sin|
// |sin * c[0] + cos * c[3] , sin * c[2] + cos * c[1]|  |-sin , cos|
static inline void rotate2x2(const double &cos, const double &sin, double * __restrict__ m22)
{
	double origin[4] = { m22[0], m22[1], m22[2], m22[3] };
	m22[0] = (cos * origin[0] - sin * origin[2]) * cos - (cos * origin[1] - sin * origin[3]) * sin;
	m22[1] = (cos * origin[0] - sin * origin[2]) * sin + (cos * origin[1] - sin * origin[3]) * cos;
	m22[2] = (sin * origin[0] + cos * origin[2]) * cos - (sin * origin[1] + cos * origin[3]) * sin;
	m22[3] = (sin * origin[0] + cos * origin[2]) * sin + (sin * origin[1] + cos * origin[3]) * cos;
}

static inline void rotate3x3(const double * __restrict__ cos, const double * __restrict__ sin, double * __restrict__ m33)
{
	double t[3][3] {
		{ cos[1] * cos[2]                           , cos[1] * sin[2]                           , -sin[1] },
		{ cos[2] * sin[0] * sin[1] - cos[0] * sin[2], cos[0] * cos[2] + sin[0] * sin[1] * sin[2], cos[1] * sin[0] },
		{ sin[0] * sin[2] + cos[0] * cos[2] * sin[1], cos[0] * sin[1] * sin[2] - cos[2] * sin[0], cos[0] * cos[1] }
	};

	double origin[9] = { m33[0], m33[1], m33[2], m33[3], m33[4], m33[5], m33[6], m33[7], m33[8] };

	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			double _a = t[0][i] * origin[0] + t[1][i] * origin[3] + t[2][i] * origin[6];
			double _b = t[0][i] * origin[1] + t[1][i] * origin[4] + t[2][i] * origin[7];
			double _c = t[0][i] * origin[2] + t[1][i] * origin[5] + t[2][i] * origin[8];
			m33[3 * i + j] = _a * t[0][j] + _b * t[1][j] + _c * t[2][j];
		}
	}
}

template<size_t nodes, size_t dimension>
static void NtoGP(const double * __restrict__ gp_N, const double * __restrict__ nvalues, double * __restrict__ gpvalues)
{
	for (size_t d = 0; d < dimension; ++d) {
		gpvalues[d] = 0;
	}
	for (size_t n = 0; n < nodes; ++n) {
		for (size_t d = 0; d < dimension; ++d) {
			gpvalues[d] += gp_N[n] * nvalues[dimension * n + d];
		}
	}
}

template<size_t nodes, size_t dimension>
static void NtoGPSimd(const double * __restrict__ gp_N, const double * __restrict__ nvalues, double * __restrict__ gpvalues)
{

	SIMD res[dimension];

	for (size_t d = 0; d < dimension; ++d) {
		res[d] = zeros();
	}

	for (size_t n = 0; n < nodes; ++n) {
		SIMD gp = load(&gp_N[n * SIMD::size]);

		for (size_t d = 0; d < dimension; ++d) {
			SIMD nv = load(&nvalues[(dimension * n + d) * SIMD::size]);
			res[d] = res[d] + gp * nv;
		}
	}

	for (size_t d = 0; d < dimension; ++d) {
		store(&gpvalues[d * SIMD::size], res[d]);
	}
}

template<size_t N>
static inline void M12M2N(const double * __restrict__ m12, const double * __restrict__ m2N, double * __restrict__ result)
{
	for (size_t n = 0; n < N; ++n) {
		result[0 + n] = m12[0] * m2N[n] + m12[1] * m2N[N + n];
	}
}

template<size_t N>
static inline void M13M3N(const double * __restrict__ m13, const double * __restrict__ m3N, double * __restrict__ result)
{
	for (size_t n = 0; n < N; ++n) {
		result[n] = m13[0] * m3N[n] + m13[1] * m3N[N + n] + m13[2] * m3N[2 * N + n];
	}
}

template<size_t N>
static inline void M22M2N(const double * __restrict__ m22, const double * __restrict__ m2N, double * __restrict__ result)
{
	for (size_t n = 0; n < N; ++n) {
		result[0 + n] = m22[0] * m2N[n] + m22[1] * m2N[N + n];
		result[N + n] = m22[2] * m2N[n] + m22[3] * m2N[N + n];
	}
}

template<size_t N>
ALWAYS_INLINE static void M22M2NSimd(const double * __restrict__ m22, const double * __restrict__ m2N, double * __restrict__ result)
{
	SIMD m22Simd0 = load(&m22[0 * SIMD::size]);
	SIMD m22Simd1 = load(&m22[1 * SIMD::size]);
	SIMD m22Simd2 = load(&m22[2 * SIMD::size]);
	SIMD m22Simd3 = load(&m22[3 * SIMD::size]);

	#pragma unroll(N)
	for (size_t n = 0; n < N; ++n) {
		SIMD m2N1 = load(&m2N[(0 + n) * SIMD::size]);
		SIMD m2N2 = load(&m2N[(N + n) * SIMD::size]);

		store(&result[(0 + n) * SIMD::size], m22Simd0 * m2N1 + m22Simd1 * m2N2);
		store(&result[(N + n) * SIMD::size], m22Simd2 * m2N1 + m22Simd3 * m2N2);
	}
}

template<size_t N>
static inline void M33M3N(const double * __restrict__ m33, const double * __restrict__ m3N, double * __restrict__ result)
{
	for (size_t n = 0; n < N; ++n) {
		result[0 * N + n] = m33[0] * m3N[n] + m33[1] * m3N[N + n] + m33[2] * m3N[2 * N + n];
		result[1 * N + n] = m33[3] * m3N[n] + m33[4] * m3N[N + n] + m33[5] * m3N[2 * N + n];
		result[2 * N + n] = m33[6] * m3N[n] + m33[7] * m3N[N + n] + m33[8] * m3N[2 * N + n];
	}
}

template<size_t N>
ALWAYS_INLINE static void M33M3NSimd(const double * __restrict__ m33, const double * __restrict__ m3N, double * __restrict__ result)
{
	SIMD m33Simd0 = load(&m33[0 * SIMD::size]);
	SIMD m33Simd1 = load(&m33[1 * SIMD::size]);
	SIMD m33Simd2 = load(&m33[2 * SIMD::size]);
	SIMD m33Simd3 = load(&m33[3 * SIMD::size]);
	SIMD m33Simd4 = load(&m33[4 * SIMD::size]);
	SIMD m33Simd5 = load(&m33[5 * SIMD::size]);
	SIMD m33Simd6 = load(&m33[6 * SIMD::size]);
	SIMD m33Simd7 = load(&m33[7 * SIMD::size]);
	SIMD m33Simd8 = load(&m33[8 * SIMD::size]);

	#pragma unroll(N)
	for (size_t n = 0; n < N; ++n) {
		SIMD m3N1 = load(&m3N[(0 * N + n) * SIMD::size]);
		SIMD m3N2 = load(&m3N[(1 * N + n) * SIMD::size]);
		SIMD m3N3 = load(&m3N[(2 * N + n) * SIMD::size]);

		store(&result[(0 * N + n) * SIMD::size], m33Simd0 * m3N1 + m33Simd1 * m3N2 + m33Simd2 * m3N3);
		store(&result[(1 * N + n) * SIMD::size], m33Simd3 * m3N1 + m33Simd4 * m3N2 + m33Simd5 * m3N3);
		store(&result[(2 * N + n) * SIMD::size], m33Simd6 * m3N1 + m33Simd7 * m3N2 + m33Simd8 * m3N3);
	}
}

template<size_t N>
static inline void MN1M1N(const double &sumscale, const double * __restrict__ mN1, const double * __restrict__ m1N, double * __restrict__ mNN)
{
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			mNN[n * N + m] = sumscale * mN1[n] * m1N[m];
		}
	}
}

template<size_t N>
static inline void M1NMN2(const double &sumscale, const double * __restrict__ m1N, const double * __restrict__ mN2, double * __restrict__ m12)
{
	m12[2 * 0 + 0] = 0;
	m12[2 * 0 + 1] = 0;
	for (size_t n = 0; n < N; ++n) {
		m12[2 * 0 + 0] += m1N[N * 0 + n] * mN2[2 * n + 0];
		m12[2 * 0 + 1] += m1N[N * 0 + n] * mN2[2 * n + 1];
	}
	m12[2 * 0 + 0] *= sumscale;
	m12[2 * 0 + 1] *= sumscale;
}

template<size_t N>
static inline void M1NMN3(const double &sumscale, const double * __restrict__ m1N, const double * __restrict__ mN3, double * __restrict__ m13)
{
	m13[3 * 0 + 0] = 0;
	m13[3 * 0 + 1] = 0;
	m13[3 * 0 + 2] = 0;
	for (size_t n = 0; n < N; ++n) {
		m13[3 * 0 + 0] += m1N[N * 0 + n] * mN3[3 * n + 0];
		m13[3 * 0 + 1] += m1N[N * 0 + n] * mN3[3 * n + 1];
		m13[3 * 0 + 2] += m1N[N * 0 + n] * mN3[3 * n + 2];
	}
	m13[3 * 0 + 0] *= sumscale;
	m13[3 * 0 + 1] *= sumscale;
	m13[3 * 0 + 2] *= sumscale;
}

template<size_t N>
static inline void M2NMN3(const double &sumscale, const double * __restrict__ m2N, const double * __restrict__ mN3, double * __restrict__ m23)
{
	m23[3 * 0 + 0] = 0;
	m23[3 * 0 + 1] = 0;
	m23[3 * 0 + 2] = 0;
	m23[3 * 1 + 0] = 0;
	m23[3 * 1 + 1] = 0;
	m23[3 * 1 + 2] = 0;
	for (size_t n = 0; n < N; ++n) {
		m23[3 * 0 + 0] += m2N[N * 0 + n] * mN3[3 * n + 0];
		m23[3 * 0 + 1] += m2N[N * 0 + n] * mN3[3 * n + 1];
		m23[3 * 0 + 2] += m2N[N * 0 + n] * mN3[3 * n + 2];
		m23[3 * 1 + 0] += m2N[N * 1 + n] * mN3[3 * n + 0];
		m23[3 * 1 + 1] += m2N[N * 1 + n] * mN3[3 * n + 1];
		m23[3 * 1 + 2] += m2N[N * 1 + n] * mN3[3 * n + 2];
	}
	m23[3 * 0 + 0] *= sumscale;
	m23[3 * 0 + 1] *= sumscale;
	m23[3 * 0 + 2] *= sumscale;
	m23[3 * 1 + 0] *= sumscale;
	m23[3 * 1 + 1] *= sumscale;
	m23[3 * 1 + 2] *= sumscale;
}

template<size_t N>
static inline void ADDM2NMN1(const double &sumscale, const double * __restrict__ m2N, const double * __restrict__ mN1, double * __restrict__ m21)
{
	double x[2] = { 0, 0 };
	for (size_t n = 0; n < N; ++n) {
		x[0] += m2N[0 * N + n] * mN1[n];
		x[1] += m2N[1 * N + n] * mN1[n];
	}
	m21[0] += sumscale * x[0];
	m21[1] += sumscale * x[1];
}

template<size_t N>
static inline void ADDM3NMN1(const double &sumscale, const double * __restrict__ m3N, const double * __restrict__ mN1, double * __restrict__ m31)
{
	double x[3] = { 0, 0, 0 };
	for (size_t n = 0; n < N; ++n) {
		x[0] += m3N[0 * N + n] * mN1[n];
		x[1] += m3N[1 * N + n] * mN1[n];
		x[2] += m3N[2 * N + n] * mN1[n];
	}
	m31[0] += sumscale * x[0];
	m31[1] += sumscale * x[1];
	m31[2] += sumscale * x[2];
}

template<size_t N>
static inline void ADDM22M2NMN1(const double &sumscale, const double * __restrict__ m22, const double * __restrict__ m2N, const double * __restrict__ mN1, double * __restrict__ m21)
{
	double x[2] = { 0, 0 };
	for (size_t n = 0; n < N; ++n) {
		x[0] += m2N[0 * N + n] * mN1[n];
		x[1] += m2N[1 * N + n] * mN1[n];
	}
	m21[0] += sumscale * (m22[0 * 2 + 0] * x[0] + m22[0 * 2 + 1] * x[1]);
	m21[1] += sumscale * (m22[1 * 2 + 0] * x[0] + m22[1 * 2 + 1] * x[1]);
}

template<size_t N>
static inline void ADDM33M3NMN1(const double &sumscale, const double * __restrict__ m33, const double * __restrict__ m3N, const double * __restrict__ mN1, double * __restrict__ m31)
{
	double x[3] = { 0, 0, 0 };
	for (size_t n = 0; n < N; ++n) {
		x[0] += m3N[0 * N + n] * mN1[n];
		x[1] += m3N[1 * N + n] * mN1[n];
		x[2] += m3N[2 * N + n] * mN1[n];
	}
	m31[0] += sumscale * (m33[0 * 3 + 0] * x[0] + m33[0 * 3 + 1] * x[1] + m33[0 * 3 + 2] * x[2]);
	m31[1] += sumscale * (m33[1 * 3 + 0] * x[0] + m33[1 * 3 + 1] * x[1] + m33[1 * 3 + 2] * x[2]);
	m31[2] += sumscale * (m33[2 * 3 + 0] * x[0] + m33[2 * 3 + 1] * x[1] + m33[2 * 3 + 2] * x[2]);
}

template<size_t N>
static inline void ADDMN1M1N(const double &sumscale, const double * __restrict__ mN1, const double * __restrict__ m1N, double * __restrict__ mNN)
{
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (mN1[n] * m1N[m]);
		}
	}
}

template<size_t N>
static inline void ADDMN1M1N(const double &sumscale, const double * __restrict__ m1N, double * __restrict__ mNN)
{
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (m1N[n] * m1N[m]);
		}
	}
}

template<size_t N>
static inline void ADDMN2M2N(const double &sumscale, const double * __restrict__ m2N, double * __restrict__ mNN)
{
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * m2N[n] * m2N[m];
			mNN[n * N + m] += sumscale * m2N[N + n] * m2N[N + m];
		}
	}
}

template<size_t N>
ALWAYS_INLINE static void ADDMN2M2NSimd(const SIMD &sumscale, const double * __restrict__ m2N, double * __restrict__ mNN)
{
	for (size_t n = 0; n < N; ++n) {
		SIMD m2N1 = load(&m2N[n * SIMD::size]);
		SIMD m2N2 = load(&m2N[(N + n) * SIMD::size]);
		for (size_t m = 0; m < N; ++m) {
			SIMD m2N3 = load(&m2N[m * SIMD::size]);
			SIMD m2N4 = load(&m2N[(N + m) * SIMD::size]);
			SIMD res  = load(&mNN[(n * N + m) * SIMD::size]);
			res = res + sumscale * (m2N1 * m2N3 + m2N2 * m2N4);
			store(&mNN[(n * N + m) * SIMD::size], res);
		}
	}
}

template<size_t N>
static inline void ADDMN3M3N(const double &sumscale, const double * __restrict__ m3N, double * __restrict__ mNN)
{
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (m3N[n] * m3N[m] + m3N[N + n] * m3N[N + m] + m3N[2 * N + n] * m3N[2 * N + m]);
		}
	}
}

template<size_t N>
static inline void ADDMN2M22M2N(const double &sumscale, const double * __restrict__ m22, const double * __restrict__ m2N, double * __restrict__ mNN)
{
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = m2N[n] * m22[0] + m2N[N + n] * m22[2];
			double b = m2N[n] * m22[1] + m2N[N + n] * m22[3];
			mNN[n * N + m] += sumscale * (a * m2N[m] + b * m2N[N + m]);
		}
	}
}

template<size_t N>
static inline void ADDMN3M33M3N(const double &sumscale, const double * __restrict__ m33, const double * __restrict__ m3N, double * __restrict__ mNN)
{
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = m3N[n] * m33[0] + m3N[N + n] * m33[3] + m3N[2 * N + n] * m33[6];
			double b = m3N[n] * m33[1] + m3N[N + n] * m33[4] + m3N[2 * N + n] * m33[7];
			double c = m3N[n] * m33[2] + m3N[N + n] * m33[5] + m3N[2 * N + n] * m33[8];
			mNN[n * N + m] += sumscale * (a * m3N[m] + b * m3N[N + m] + c * m3N[2 * N + m]);
		}
	}
}

// B * C * Bt
//
// C = 3x3
// B = dX  0 dY
//      0 dY dX
template<size_t N>
static inline void ADDDXDY33DXDYN(const double &sumscale, const double * __restrict__ m33, const double * __restrict__ dX, double * __restrict__ mNN)
{
	const double * __restrict__ dY = dX + N;
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = dX[n] * m33[0] + dY[n] * m33[6];
			double b = dX[n] * m33[1] + dY[n] * m33[7];
			double c = dX[n] * m33[2] + dY[n] * m33[8];
			mNN[n * 2 * N + m    ] += sumscale * (a * dX[m] + c * dY[m]);
			mNN[n * 2 * N + m + N] += sumscale * (b * dY[m] + c * dX[m]);
		}
	}
	mNN += 2 * N * N;
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = dY[n] * m33[3] + dX[n] * m33[6];
			double b = dY[n] * m33[4] + dX[n] * m33[7];
			double c = dY[n] * m33[5] + dX[n] * m33[8];
			mNN[n * 2 * N + m    ] += sumscale * (a * dX[m] + c * dY[m]);
			mNN[n * 2 * N + m + N] += sumscale * (b * dY[m] + c * dX[m]);
		}
	}
}

// B * C * Bt
//
// C = 4x4
// B = dX  0  C dY
//      0 dY  0 dX
template<size_t N>
static inline void ADDDXDYCOO44DXDYN(const double &sumscale, const double * __restrict__ m44, const double * __restrict__ dX, const double * __restrict__ coo, double * __restrict__ mNN)
{
	const double * __restrict__ dY = dX + N;
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = dX[n] * m44[0] + coo[n] * m44[ 8] + dY[n] * m44[12];
			double b = dX[n] * m44[1] + coo[n] * m44[ 9] + dY[n] * m44[13];
			double c = dX[n] * m44[2] + coo[n] * m44[10] + dY[n] * m44[14];
			double d = dX[n] * m44[3] + coo[n] * m44[11] + dY[n] * m44[15];
			mNN[n * 2 * N + m    ] += sumscale * (a * dX[m] + c * coo[m] + d * dY[m]);
			mNN[n * 2 * N + m + N] += sumscale * (b * dY[m] + d * dX[m]);
		}
	}
	mNN += 2 * N * N;
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = dY[n] * m44[4] + dX[n] * m44[12];
			double b = dY[n] * m44[5] + dX[n] * m44[13];
			double c = dY[n] * m44[6] + dX[n] * m44[14];
			double d = dY[n] * m44[7] + dX[n] * m44[15];
			mNN[n * 2 * N + m    ] += sumscale * (a * dX[m] + c * coo[m] + d * dY[m]);
			mNN[n * 2 * N + m + N] += sumscale * (b * dY[m] + d * dX[m]);
		}
	}
}

// B * C * Bt
//
// C = 6x6
// B = dX  0  0 dY  0 dZ
//      0 dY  0 dX dZ  0
//      0  0 dZ  0 dY dX
//     - - - - - - - - -
//      0  6 12 18 24 30
//      a  b  c  d  e  f
template<size_t N>
static inline void ADDDXDYDZ66DXDYDZN(const double &sumscale, const double * __restrict__ m66, const double * __restrict__ dX, double * __restrict__ mNN)
{
	const double * __restrict__ dY = dX + N;
	const double * __restrict__ dZ = dX + 2 * N;
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = dX[n] * m66[ 0] + dY[n] * m66[18] + dZ[n] * m66[30];
			double b = dX[n] * m66[ 1] + dY[n] * m66[19] + dZ[n] * m66[31];
			double c = dX[n] * m66[ 2] + dY[n] * m66[20] + dZ[n] * m66[32];
			double d = dX[n] * m66[ 3] + dY[n] * m66[21] + dZ[n] * m66[33];
			double e = dX[n] * m66[ 4] + dY[n] * m66[22] + dZ[n] * m66[34];
			double f = dX[n] * m66[ 5] + dY[n] * m66[23] + dZ[n] * m66[35];
			mNN[n * 3 * N + m + 0 * N] += sumscale * (a * dX[m] + d * dY[m] + f * dZ[m]);
			mNN[n * 3 * N + m + 1 * N] += sumscale * (b * dY[m] + d * dX[m] + e * dZ[m]);
			mNN[n * 3 * N + m + 2 * N] += sumscale * (c * dZ[m] + e * dY[m] + f * dX[m]);
		}
	}
	mNN += 3 * N * N;
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = dY[n] * m66[ 6] + dX[n] * m66[18] + dZ[n] * m66[24];
			double b = dY[n] * m66[ 7] + dX[n] * m66[19] + dZ[n] * m66[25];
			double c = dY[n] * m66[ 8] + dX[n] * m66[20] + dZ[n] * m66[26];
			double d = dY[n] * m66[ 9] + dX[n] * m66[21] + dZ[n] * m66[27];
			double e = dY[n] * m66[10] + dX[n] * m66[22] + dZ[n] * m66[28];
			double f = dY[n] * m66[11] + dX[n] * m66[23] + dZ[n] * m66[29];
			mNN[n * 3 * N + m + 0 * N] += sumscale * (a * dX[m] + d * dY[m] + f * dZ[m]);
			mNN[n * 3 * N + m + 1 * N] += sumscale * (b * dY[m] + d * dX[m] + e * dZ[m]);
			mNN[n * 3 * N + m + 2 * N] += sumscale * (c * dZ[m] + e * dY[m] + f * dX[m]);
		}
	}
	mNN += 3 * N * N;
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < N; ++m) {
			double a = dZ[n] * m66[12] + dY[n] * m66[24] + dX[n] * m66[30];
			double b = dZ[n] * m66[13] + dY[n] * m66[25] + dX[n] * m66[31];
			double c = dZ[n] * m66[14] + dY[n] * m66[26] + dX[n] * m66[32];
			double d = dZ[n] * m66[15] + dY[n] * m66[27] + dX[n] * m66[33];
			double e = dZ[n] * m66[16] + dY[n] * m66[28] + dX[n] * m66[34];
			double f = dZ[n] * m66[17] + dY[n] * m66[29] + dX[n] * m66[35];
			mNN[n * 3 * N + m + 0 * N] += sumscale * (a * dX[m] + d * dY[m] + f * dZ[m]);
			mNN[n * 3 * N + m + 1 * N] += sumscale * (b * dY[m] + d * dX[m] + e * dZ[m]);
			mNN[n * 3 * N + m + 2 * N] += sumscale * (c * dZ[m] + e * dY[m] + f * dX[m]);
		}
	}
}


template<size_t N>
static inline void ADDM12M2N(const double &sumscale, const double * __restrict__ m12, const double * __restrict__ m2N, double * __restrict__ result)
{
	for (size_t n = 0; n < N; ++n) {
		result[0 + n] += sumscale * (m12[0] * m2N[n] + m12[1] * m2N[N + n]);
	}
}

template<size_t N>
static inline void ADDM13M3N(const double &sumscale, const double * __restrict__ m13, const double * __restrict__ m3N, double * __restrict__ result)
{
	for (size_t n = 0; n < N; ++n) {
		result[n] += sumscale * (m13[0] * m3N[n] + m13[1] * m3N[N + n] + m13[2] * m3N[2 * N + n]);
	}
}

}

#endif /* SRC_PHYSICS_ASSEMBLER_MATH_HPP_ */
