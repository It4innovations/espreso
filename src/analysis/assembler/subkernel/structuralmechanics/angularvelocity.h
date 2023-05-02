
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ANGULARVELOCITY_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ANGULARVELOCITY_H_

#include "analysis/assembler/subkernel/boundarycondition.h"

namespace espreso {

template <size_t nodes, size_t gps, enum Behaviour behaviour, class Physics> struct AngularVelocityKernel;

template <size_t nodes, size_t gps, class Physics>
struct AngularVelocityKernel<nodes, gps, Behaviour::PLANE, Physics>: BoundaryCondition, Physics {
	AngularVelocityKernel(const BoundaryCondition &base): BoundaryCondition(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]) * element.ecf.angularVelocity[gp][0] * element.ecf.angularVelocity[gp][0];
				fx = fx + scale * element.gpcoords[gp][0];
				fy = fy + scale * element.gpcoords[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		rhs += SIMD::size * 2 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct AngularVelocityKernel<nodes, gps, Behaviour::AXISYMMETRIC, Physics>: BoundaryCondition, Physics {
	AngularVelocityKernel(const BoundaryCondition &base): BoundaryCondition(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0] * element.ecf.density[gp] * load1(element.N[gp][n]) * element.ecf.angularVelocity[gp][0] * element.ecf.angularVelocity[gp][0];
				fx = fx + scale * element.gpcoords[gp][0];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
		}
		rhs += SIMD::size * 2 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct AngularVelocityKernel<nodes, gps, Behaviour::VOLUME, Physics>: BoundaryCondition, Physics {
	AngularVelocityKernel(const BoundaryCondition &base): BoundaryCondition(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			SIMD fz = load(out + (2 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]);
				SIMD vx = element.ecf.angularVelocity[gp][0] * element.ecf.angularVelocity[gp][0];
				SIMD vy = element.ecf.angularVelocity[gp][1] * element.ecf.angularVelocity[gp][1];
				SIMD vz = element.ecf.angularVelocity[gp][2] * element.ecf.angularVelocity[gp][2];
				fx = fx + scale * element.gpcoords[gp][0] * (vy + vz);
				fy = fy + scale * element.gpcoords[gp][1] * (vx + vz);
				fz = fz + scale * element.gpcoords[gp][2] * (vx + vy);
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
			store(out + (2 * nodes + n) * SIMD::size, fz);
		}
		rhs += SIMD::size * 3 * nodes;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ANGULARVELOCITY_H_ */
