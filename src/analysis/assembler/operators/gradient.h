
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_

#include "analysis/assembler/operator.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct TemperatureGradientBase: ActionOperator {
	TemperatureGradientBase(size_t interval, NamedData *gradient)
	: gradient(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin)
	{
		isconst = false;
		action = Action::SOLUTION;
	}

	double* gradient;

	void move(int n)
	{
		gradient += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct TemperatureGradient;

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct TemperatureGradient<nodes, gps, 2, edim, etype, Physics>: TemperatureGradientBase, Physics {
	using TemperatureGradientBase::TemperatureGradientBase;

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = gradient;
		SIMD g0 = zeros(), g1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
			}
		}
		SIMD scale = load1(1. / gps);
		g0 = g0 * scale;
		g1 = g1 * scale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			out[2 * s + 0] = g0[s];
			out[2 * s + 1] = g1[s];
		}
		move(2 * SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		double * __restrict__ out = gradient;
		SIMD g0 = zeros(), g1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
			}
		}
		SIMD scale = load1(1. / gps);
		g0 = g0 * scale;
		g1 = g1 * scale;
		for (size_t s = 0; s < size; ++s) {
			out[2 * s + 0] = g0[s];
			out[2 * s + 1] = g1[s];
		}
		move(2 * size);
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct TemperatureGradient<nodes, gps, 3, edim, etype, Physics>: TemperatureGradientBase, Physics {
	using TemperatureGradientBase::TemperatureGradientBase;

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = gradient;
		SIMD g0 = zeros(), g1 = zeros(), g2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
				g2 = g2 + element.dND[gp][n][2] * element.temp[n];
			}
		}
		SIMD scale = load1(1. / gps);
		g0 = g0 * scale;
		g1 = g1 * scale;
		g2 = g2 * scale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			out[3 * s + 0] = g0[s];
			out[3 * s + 1] = g1[s];
			out[3 * s + 2] = g2[s];
		}
		move(3 * SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		double * __restrict__ out = gradient;
		SIMD g0 = zeros(), g1 = zeros(), g2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
				g2 = g2 + element.dND[gp][n][2] * element.temp[n];
			}
		}
		SIMD scale = load1(1. / gps);
		g0 = g0 * scale;
		g1 = g1 * scale;
		g2 = g2 * scale;
		for (size_t s = 0; s < size; ++s) {
			out[3 * s + 0] = g0[s];
			out[3 * s + 1] = g1[s];
			out[3 * s + 2] = g2[s];
		}
		move(3 * size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_ */
