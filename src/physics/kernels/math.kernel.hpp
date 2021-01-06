
#ifndef SRC_PHYSICS_KERNELS_MATH_KERNEL_HPP_
#define SRC_PHYSICS_KERNELS_MATH_KERNEL_HPP_

#include "math/math.h"
#include "math/matrix.dense.h"
#include <vector>

namespace espreso {

// |cos , -sin| |c[0] , c[2]| | cos , sin|
// |sin ,  cos| |c[3] , c[1]| |-sin , cos|
// |cos * c[0] - sin * c[3] , cos * c[2] - sin * c[1]|  | cos , sin|
// |sin * c[0] + cos * c[3] , sin * c[2] + cos * c[1]|  |-sin , cos|
static inline void rotate2Dfull(const double &cos, const double &sin, const double * __restrict__ c, double * __restrict__ k)
{
	k[0] = (cos * c[0] - sin * c[3]) * cos - (cos * c[2] - sin * c[1]) * sin;
	k[1] = (cos * c[0] - sin * c[3]) * sin + (cos * c[2] - sin * c[1]) * cos;
	k[2] = (sin * c[0] + cos * c[3]) * cos - (sin * c[2] + cos * c[1]) * sin;
	k[3] = (sin * c[0] + cos * c[3]) * sin + (sin * c[2] + cos * c[1]) * cos;
}

static inline void rotate2Dsym(const double &cos, const double &sin, const double * __restrict__ c, double * __restrict__ k)
{
	k[0] = (cos * c[0] - sin * c[2]) * cos - (cos * c[2] - sin * c[1]) * sin;
	k[1] = (cos * c[0] - sin * c[2]) * sin + (cos * c[2] - sin * c[1]) * cos;
	k[2] = (sin * c[0] + cos * c[2]) * cos - (sin * c[2] + cos * c[1]) * sin;
	k[3] = (sin * c[0] + cos * c[2]) * sin + (sin * c[2] + cos * c[1]) * cos;
}

static inline void rotate2Ddiag(const double &cos, const double &sin, const double * __restrict__ c, double * __restrict__ k)
{
	k[0] = (cos * c[0]) * cos - (- sin * c[1]) * sin;
	k[1] = (cos * c[0]) * sin + (- sin * c[1]) * cos;
	k[2] = (sin * c[0]) * cos - (+ cos * c[1]) * sin;
	k[3] = (sin * c[0]) * sin + (+ cos * c[1]) * cos;
}

static inline void rotate3Dfull(const double * __restrict__ cos, const double * __restrict__ sin, const double * __restrict__ c, double * __restrict__ k)
{
	double t[3][3] {
		{ cos[1] * cos[2]                           , cos[1] * sin[2]                           , -sin[1] },
		{ cos[2] * sin[0] * sin[1] - cos[0] * sin[2], cos[0] * cos[2] + sin[0] * sin[1] * sin[2], cos[1] * sin[0] },
		{ sin[0] * sin[2] + cos[0] * cos[2] * sin[1], cos[0] * sin[1] * sin[2] - cos[2] * sin[0], cos[0] * cos[1] }
	};

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			double _a = t[0][i] * c[0] + t[1][i] * c[6] + t[2][i] * c[8];
			double _b = t[0][i] * c[3] + t[1][i] * c[1] + t[2][i] * c[7];
			double _c = t[0][i] * c[5] + t[1][i] * c[4] + t[2][i] * c[2];
			k[3 * i + j] = _a * t[0][j] + _b * t[1][j] + _c * t[2][j];
		}
	}
}

static inline void rotate3Dsym(const double * __restrict__ cos, const double * __restrict__ sin, const double * __restrict__ c, double * __restrict__ k)
{
	double t[3][3] {
		{ cos[1] * cos[2]                           , cos[1] * sin[2]                           , -sin[1] },
		{ cos[2] * sin[0] * sin[1] - cos[0] * sin[2], cos[0] * cos[2] + sin[0] * sin[1] * sin[2], cos[1] * sin[0] },
		{ sin[0] * sin[2] + cos[0] * cos[2] * sin[1], cos[0] * sin[1] * sin[2] - cos[2] * sin[0], cos[0] * cos[1] }
	};

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			double _a = t[0][i] * c[0] + t[1][i] * c[3] + t[2][i] * c[5];
			double _b = t[0][i] * c[3] + t[1][i] * c[1] + t[2][i] * c[4];
			double _c = t[0][i] * c[5] + t[1][i] * c[4] + t[2][i] * c[2];
			k[3 * i + j] = _a * t[0][j] + _b * t[1][j] + _c * t[2][j];
		}
	}
}

static inline void rotate3Ddiag(const double * __restrict__ cos, const double * __restrict__ sin, const double * __restrict__ c, double * __restrict__ k)
{
	double t[3][3] {
		{ cos[1] * cos[2]                           , cos[1] * sin[2]                           ,         -sin[1] },
		{ cos[2] * sin[0] * sin[1] - cos[0] * sin[2], cos[0] * cos[2] + sin[0] * sin[1] * sin[2], cos[1] * sin[0] },
		{ sin[0] * sin[2] + cos[0] * cos[2] * sin[1], cos[0] * sin[1] * sin[2] - cos[2] * sin[0], cos[0] * cos[1] }
	};

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			double _a = t[0][i] * c[0];
			double _b = t[1][i] * c[1];
			double _c = t[2][i] * c[2];
			k[3 * i + j] = _a * t[0][j] + _b * t[1][j] + _c * t[2][j];
		}
	}
}

template<int nodes, int dimension>
static void NtoGP(const double * __restrict__ gp_N, const double * __restrict__ nvalues, double * __restrict__ gpvalues)
{
	for (int d = 0; d < dimension; ++d) {
		gpvalues[d] = 0;
	}
	for (int n = 0; n < nodes; ++n) {
		for (int d = 0; d < dimension; ++d) {
			gpvalues[d] += gp_N[n] * nvalues[dimension * n + d];
		}
	}
}

template<int N>
static inline void M12M2N(const double * __restrict__ m12, const double * __restrict__ m2N, double * __restrict__ result)
{
	for (int n = 0; n < N; ++n) {
		result[0 + n] = m12[0] * m2N[n] + m12[1] * m2N[N + n];
	}
}

template<int N>
static inline void M13M3N(const double * __restrict__ m13, const double * __restrict__ m3N, double * __restrict__ result)
{
	for (int n = 0; n < N; ++n) {
		result[n] = m13[0] * m3N[n] + m13[1] * m3N[N + n] + m13[2] * m3N[2 * N + n];
	}
}

template<int N>
static inline void M22M2N(const double * __restrict__ m22, const double * __restrict__ m2N, double * __restrict__ result)
{
	for (int n = 0; n < N; ++n) {
		result[0 + n] = m22[0] * m2N[n] + m22[1] * m2N[N + n];
		result[N + n] = m22[2] * m2N[n] + m22[3] * m2N[N + n];
	}
}

template<int N>
static inline void M33M3N(const double * __restrict__ m33, const double * __restrict__ m3N, double * __restrict__ result)
{
	for (int n = 0; n < N; ++n) {
		result[0 * N + n] = m33[0] * m3N[n] + m33[1] * m3N[N + n] + m33[2] * m3N[2 * N + n];
		result[1 * N + n] = m33[3] * m3N[n] + m33[4] * m3N[N + n] + m33[5] * m3N[2 * N + n];
		result[2 * N + n] = m33[6] * m3N[n] + m33[7] * m3N[N + n] + m33[8] * m3N[2 * N + n];
	}
}

template<int N>
static inline void MN1M1N(const double &sumscale, const double * __restrict__ mN1, const double * __restrict__ m1N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * mN1[n] * m1N[m];
		}
	}
}

template<int N>
static inline void KMN2M2N(const double &sumscale, const double * __restrict__ m2N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (m2N[n] * m2N[m] + m2N[N + n] * m2N[N + m]);
		}
	}
}

template<int N>
static inline void KMN3M3N(const double &sumscale, const double * __restrict__ m3N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			mNN[n * N + m] += sumscale * (m3N[n] * m3N[m] + m3N[N + n] * m3N[N + m] + m3N[2 * N + n] * m3N[2 * N + m]);
		}
	}
}

template<int N>
static inline void KMN2M22M2N(const double &sumscale, const double * __restrict__ m22, const double * __restrict__ m2N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			double a = m2N[n] * m22[0] + m2N[N + n] * m22[2];
			double b = m2N[n] * m22[1] + m2N[N + n] * m22[3];
			mNN[n * N + m] += sumscale * (a * m2N[m] + b * m2N[N + m]);
		}
	}
}

template<int N>
static inline void KMN3M33M3N(const double &sumscale, const double * __restrict__ m33, const double * __restrict__ m3N, double * __restrict__ mNN)
{
	for (int n = 0; n < N; ++n) {
		for (int m = 0; m < N; ++m) {
			double a = m3N[n] * m33[0] + m3N[N + n] * m33[3] + m3N[2 * N + n] * m33[6];
			double b = m3N[n] * m33[1] + m3N[N + n] * m33[4] + m3N[2 * N + n] * m33[7];
			double c = m3N[n] * m33[2] + m3N[N + n] * m33[5] + m3N[2 * N + n] * m33[8];
			mNN[n * N + m] += sumscale * (a * m3N[m] + b * m3N[N + m] + c * m3N[2 * N + m]);
		}
	}
}

}

#endif /* SRC_PHYSICS_KERNELS_MATH_KERNEL_HPP_ */
