
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_NORMALPRESSURE_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_NORMALPRESSURE_H_

#include "rhs.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behavour, class Physics> struct NormalPressureKernel;

template <size_t nodes, size_t gps, class Physics>
struct NormalPressureKernel<nodes, gps, 2, 1, Behaviour::PLANE, Physics>: StructuralMechanicsRHS, Physics {
	NormalPressureKernel(const StructuralMechanicsRHS &base): StructuralMechanicsRHS(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD pressure = element.det[gp] * load1(element.w[gp]) * load1(element.N[gp][n]) * element.ecf.normalPressure[gp];
				fx = fx + pressure * element.normal[gp][0];
				fy = fy + pressure * element.normal[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		rhs += SIMD::size * 2 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct NormalPressureKernel<nodes, gps, 2, 1, Behaviour::AXISYMMETRIC, Physics>: StructuralMechanicsRHS, Physics {
	NormalPressureKernel(const StructuralMechanicsRHS &base): StructuralMechanicsRHS(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD pressure = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0] * load1(element.N[gp][n]) * element.ecf.normalPressure[gp];
				fx = fx + pressure * element.normal[gp][0];
				fy = fy + pressure * element.normal[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		rhs += SIMD::size * 2 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct NormalPressureKernel<nodes, gps, 3, 2, Behaviour::VOLUME, Physics>: StructuralMechanicsRHS, Physics {
	NormalPressureKernel(const StructuralMechanicsRHS &base): StructuralMechanicsRHS(base) {}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			SIMD fz = load(out + (2 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD pressure = element.det[gp] * load1(element.w[gp]) * load1(element.N[gp][n]) * element.ecf.normalPressure[gp];
				fx = fx + pressure * element.normal[gp][0];
				fy = fy + pressure * element.normal[gp][1];
				fz = fz + pressure * element.normal[gp][2];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
			store(out + (2 * nodes + n) * SIMD::size, fz);
		}
		rhs += SIMD::size * 3 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct NormalPressureKernel<nodes, gps, 3, 1, Behaviour::VOLUME, Physics>: StructuralMechanicsRHS, Physics {
	NormalPressureKernel(const StructuralMechanicsRHS &base): StructuralMechanicsRHS(base) {}

	void simd(typename Physics::Element &element)
	{

	}
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_NORMALPRESSURE_H_ */
