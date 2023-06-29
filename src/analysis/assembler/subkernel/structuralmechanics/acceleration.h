
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ACCELERATION_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ACCELERATION_H_

#include "analysis/assembler/subkernel/boundarycondition.h"
#include "analysis/assembler/subkernel/structuralmechanics/subkernel.h"

namespace espreso {

template <size_t nodes, size_t gps, enum Behaviour behaviour, class Physics> struct AccelerationKernel;

template <size_t nodes, size_t gps, class Physics>
struct AccelerationKernel<nodes, gps, Behaviour::PLANE, Physics>: BoundaryCondition, Physics {
	AccelerationKernel(const BoundaryCondition &base): BoundaryCondition(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]);
				fx = fx + scale * element.ecf.acceleration[gp][0];
				fy = fy + scale * element.ecf.acceleration[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		rhs += SIMD::size * 2 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct AccelerationKernel<nodes, gps, Behaviour::AXISYMMETRIC, Physics>: BoundaryCondition, Physics {
	AccelerationKernel(const BoundaryCondition &base): BoundaryCondition(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
//			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0] * element.ecf.density[gp] * load1(element.N[gp][n]);
//				fx = fx + scale * element.ecf.acceleration[gp][0];
				fy = fy + scale * element.ecf.acceleration[gp][1];
			}
//			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		rhs += SIMD::size * 2 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct AccelerationKernel<nodes, gps, Behaviour::VOLUME, Physics>: BoundaryCondition, Physics {
	AccelerationKernel(const BoundaryCondition &base): BoundaryCondition(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			SIMD fz = load(out + (2 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]);
				fx = fx + scale * element.ecf.acceleration[gp][0];
				fy = fy + scale * element.ecf.acceleration[gp][1];
				fz = fz + scale * element.ecf.acceleration[gp][2];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
			store(out + (2 * nodes + n) * SIMD::size, fz);
		}
		rhs += SIMD::size * 3 * nodes;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ACCELERATION_H_ */
