
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"
#include "analysis/assembler/module/heattransfer.element.h"

namespace espreso {

struct HeatTransferStiffnessBase: public ActionOperator {
	OutputParameterIterator stiffness;

	HeatTransferStiffnessBase(size_t interval, ParameterData &stiffness)
	: stiffness(stiffness, interval)
	{
		isconst = false;
	}

	void move(int n)
	{
		stiffness += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct HeatTransferStiffness;

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double scale = element.ecf.thickness[gp] * element.det[gp] * element.w[gp] * element.conductivity[gp];
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t m = 0; m < nodes; ++m) {
					stiffness[n * nodes + m] += scale * (element.dND[gp][n][0] * element.dND[gp][m][0] + element.dND[gp][n][1] * element.dND[gp][m][1]);
				}
			}
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.ecf.thickness[gp] * element.det[gp] * element.w[gp] * element.conductivity[gp];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD res = load(stiffness.data + (n * nodes + m) * SIMD::size);
					res = res + scale * (nx * mx + ny * my);
					store(stiffness.data + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double scale = element.det[gp] * element.w[gp] * element.conductivity[gp];
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t m = 0; m < nodes; ++m) {
					stiffness[n * nodes + m] += scale * (
							element.dND[gp][n][0] * element.dND[gp][m][0] +
							element.dND[gp][n][1] * element.dND[gp][m][1] +
							element.dND[gp][n][2] * element.dND[gp][m][2]);
				}
			}
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = stiffness.data;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * element.w[gp] * element.conductivity[gp];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD nz = element.dND[gp][n][2];
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD res = load(out + (n * nodes + m) * SIMD::size);

					res = res + scale * (nx * mx + ny * my + nz * mz);
					store(out + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double scale = element.ecf.thickness[gp] * element.det[gp] * element.w[gp];
			for (size_t n = 0; n < nodes; ++n) {
				double a = element.dND[gp][n][0] * element.conductivity[gp][0] + element.dND[gp][n][1] * element.conductivity[gp][2];
				double b = element.dND[gp][n][0] * element.conductivity[gp][1] + element.dND[gp][n][1] * element.conductivity[gp][3];
				for (size_t m = 0; m < nodes; ++m) {
					stiffness[n * nodes + m] += scale * (a * element.dND[gp][m][0] + b * element.dND[gp][m][1]);
				}
			}
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = stiffness.data;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.conductivity[gp][0], c01 = element.conductivity[gp][1];
			SIMD c10 = element.conductivity[gp][2], c11 = element.conductivity[gp][3];
			double scale = element.ecf.thickness[gp] * element.det[gp] * element.w[gp];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD a = nx * c00 + ny * c01;
				SIMD b = nx * c10 + ny * c11;
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD res = load(stiffness.data + (n * nodes + m) * SIMD::size);
					res = res + scale * (a * mx + b * my);
					store(stiffness.data + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double scale = element.det[gp] * element.w[gp];
			for (size_t n = 0; n < nodes; ++n) {
				double a = element.dND[gp][n][0] * element.conductivity[gp][0] + element.dND[gp][n][1] * element.conductivity[gp][3] + element.dND[gp][n][2] * element.conductivity[gp][6];
				double b = element.dND[gp][n][0] * element.conductivity[gp][1] + element.dND[gp][n][1] * element.conductivity[gp][4] + element.dND[gp][n][2] * element.conductivity[gp][7];
				double c = element.dND[gp][n][0] * element.conductivity[gp][2] + element.dND[gp][n][1] * element.conductivity[gp][5] + element.dND[gp][n][2] * element.conductivity[gp][8];
				for (size_t m = 0; m < nodes; ++m) {
					stiffness[n * nodes + m] += scale * (a * element.dND[gp][m][0] + b * element.dND[gp][m][1] + c * element.dND[gp][m][2]);
				}
			}
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = stiffness.data;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * element.w[gp];
			SIMD c00 = element.conductivity[gp][0], c01 = element.conductivity[gp][1], c02 = element.conductivity[gp][2];
			SIMD c10 = element.conductivity[gp][3], c11 = element.conductivity[gp][4], c12 = element.conductivity[gp][5];
			SIMD c20 = element.conductivity[gp][6], c21 = element.conductivity[gp][7], c22 = element.conductivity[gp][8];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD nz = element.dND[gp][n][2];
				SIMD a = nx * c00 + ny * c01 + nz * c02;
				SIMD b = nx * c10 + ny * c11 + nz * c12;
				SIMD c = nx * c20 + ny * c21 + nz * c22;
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD res = load(stiffness.data + (n * nodes + m) * SIMD::size);
					res = res + scale * (a * mx + b * my + c * mz);
					store(stiffness.data + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_ */
