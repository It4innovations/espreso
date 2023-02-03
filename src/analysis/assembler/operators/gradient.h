
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_

#include "analysis/assembler/operator.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct TemperatureGradientBase: ActionOperator {
	TemperatureGradientBase(int interval, NamedData *gradient)
	: gradient(gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin)
	{
		action = Action::SOLUTION;
	}

	double* gradient;

	void move(int n)
	{
		gradient += n * info::mesh->dimension;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct TemperatureGradient;

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct TemperatureGradient<nodes, gps, 2, edim, etype, Physics>: TemperatureGradientBase, Physics {
	using TemperatureGradientBase::TemperatureGradientBase;

	constexpr static double scale = 1. / gps;

	void sisd(typename Physics::Element &element)
	{
		gradient[0] = gradient[1] = 0;
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				gradient[0] += element.dND[gp][n][0] * element.temp[n];
				gradient[1] += element.dND[gp][n][1] * element.temp[n];
			}
		}
		gradient[0] *= scale;
		gradient[1] *= scale;
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		SIMD g0 = zeros(), g1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
			}
		}
		SIMD sscale; sscale.fill(scale);
		g0 = g0 * sscale;
		g1 = g1 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			gradient[2 * s + 0] = g0[s];
			gradient[2 * s + 1] = g1[s];
		}
		move(SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		SIMD g0 = zeros(), g1 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
			}
		}
		SIMD sscale; sscale.fill(scale);
		g0 = g0 * sscale;
		g1 = g1 * sscale;
		for (size_t s = 0; s < size; ++s) {
			gradient[2 * s + 0] = g0[s];
			gradient[2 * s + 1] = g1[s];
		}
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct TemperatureGradient<nodes, gps, 3, edim, etype, Physics>: TemperatureGradientBase, Physics {
	using TemperatureGradientBase::TemperatureGradientBase;

	constexpr static double scale = 1. / gps;

	void sisd(typename Physics::Element &element)
	{
		gradient[0] = gradient[1] = gradient[2] = 0;
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				gradient[0] += element.dND[gp][n][0] * element.temp[n];
				gradient[1] += element.dND[gp][n][1] * element.temp[n];
				gradient[2] += element.dND[gp][n][2] * element.temp[n];
			}
		}
		gradient[0] *= scale;
		gradient[1] *= scale;
		gradient[2] *= scale;
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		SIMD g0 = zeros(), g1 = zeros(), g2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
				g2 = g2 + element.dND[gp][n][2] * element.temp[n];
			}
		}
		SIMD sscale; sscale.fill(scale);
		g0 = g0 * sscale;
		g1 = g1 * sscale;
		g2 = g2 * sscale;
		for (size_t s = 0; s < SIMD::size; ++s) {
			gradient[3 * s + 0] = g0[s];
			gradient[3 * s + 1] = g1[s];
			gradient[3 * s + 2] = g2[s];
		}
		move(SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		SIMD g0 = zeros(), g1 = zeros(), g2 = zeros();
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t n = 0; n < nodes; ++n) {
				g0 = g0 + element.dND[gp][n][0] * element.temp[n];
				g1 = g1 + element.dND[gp][n][1] * element.temp[n];
				g2 = g2 + element.dND[gp][n][2] * element.temp[n];
			}
		}
		SIMD sscale; sscale.fill(scale);
		g0 = g0 * sscale;
		g1 = g1 * sscale;
		g2 = g2 * sscale;
		for (size_t s = 0; s < size; ++s) {
			gradient[3 * s + 0] = g0[s];
			gradient[3 * s + 1] = g1[s];
			gradient[3 * s + 2] = g2[s];
		}
		move(size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_GRADIENT_H_ */
