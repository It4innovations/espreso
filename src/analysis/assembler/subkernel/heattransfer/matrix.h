
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_MATRIX_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_MATRIX_H_

#include "subkernels.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t etype, class Physics> struct HeatTransferStiffnessKernel;

template <size_t nodes, size_t gps, class Physics>
struct HeatTransferStiffnessKernel<nodes, gps, 2, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: HeatTransferMatrixKernel, Physics {
	HeatTransferStiffnessKernel(const HeatTransferMatrixKernel &base): HeatTransferMatrixKernel(base) {}

	void move(int n) { K += n * nodes * nodes; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * element.conductivity[gp];
			for (size_t n = 0, i = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = n; m < nodes; ++m, ++i) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];

					SIMD res = load(out + i * SIMD::size);
					res = res + scale * (nx * mx + ny * my);
					store(out + i * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct HeatTransferStiffnessKernel<nodes, gps, 3, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: HeatTransferMatrixKernel, Physics {
	HeatTransferStiffnessKernel(const HeatTransferMatrixKernel &base): HeatTransferMatrixKernel(base) {}

	void move(int n) { K += n * nodes * nodes; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * load1(element.w[gp]) * element.conductivity[gp];
			for (size_t n = 0, i = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD nz = element.dND[gp][n][2];
				for (size_t m = n; m < nodes; ++m, ++i) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];

					SIMD res = load(out + i * SIMD::size);
					res = res + scale * (nx * mx + ny * my + nz * mz);
					store(out + i * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};


template <size_t nodes, size_t gps, class Physics>
struct HeatTransferStiffnessKernel<nodes, gps, 2, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferMatrixKernel, Physics {
	HeatTransferStiffnessKernel(const HeatTransferMatrixKernel &base): HeatTransferMatrixKernel(base) {}

	void move(int n) { K += n * nodes * nodes; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.conductivity[gp][0];
			SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][2];
			SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
			for (size_t n = 0, i = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD a = nx * c00 + ny * c10;
				SIMD b = nx * c10 + ny * c11;
				for (size_t m = n; m < nodes; ++m, ++i) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];

					SIMD res = load(out + i * SIMD::size);
					res = res + scale * (a * mx + b * my);
					store(out + i * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct HeatTransferStiffnessKernel<nodes, gps, 3, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferMatrixKernel, Physics {
	HeatTransferStiffnessKernel(const HeatTransferMatrixKernel &base): HeatTransferMatrixKernel(base) {}

	void move(int n) { K += n * nodes * nodes; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * load1(element.w[gp]);
			SIMD c00 = element.conductivity[gp][0];
			SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][3];
			SIMD c20 = element.conductivity[gp][2], c21 = element.conductivity[gp][4], c22 = element.conductivity[gp][5];
			for (size_t n = 0, i = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD nz = element.dND[gp][n][2];
				SIMD a = nx * c00 + ny * c10 + nz * c20;
				SIMD b = nx * c10 + ny * c11 + nz * c21;
				SIMD c = nx * c20 + ny * c21 + nz * c22;
				for (size_t m = n; m < nodes; ++m, ++i) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];

					SIMD res = load(out + i * SIMD::size);
					res = res + scale * (a * mx + b * my + c * mz);
					store(out + i * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct HeatTransferStiffnessKernel<nodes, gps, 2, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Physics>: HeatTransferMatrixKernel, Physics {
	HeatTransferStiffnessKernel(const HeatTransferMatrixKernel &base): HeatTransferMatrixKernel(base) {}

	void move(int n) { K += n * nodes * nodes; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * element.conductivity[gp];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];

					SIMD res = load(out + (n * nodes + m) * SIMD::size);
					res = res + scale * (nx * mx + ny * my);
					store(out + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct HeatTransferStiffnessKernel<nodes, gps, 3, HeatTransferElementType::ASYMMETRIC_ISOTROPIC, Physics>: HeatTransferMatrixKernel, Physics {
	HeatTransferStiffnessKernel(const HeatTransferMatrixKernel &base): HeatTransferMatrixKernel(base) {}

	void move(int n) { K += n * nodes * nodes; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * load1(element.w[gp]) * element.conductivity[gp];
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

template <size_t nodes, size_t gps, class Physics>
struct HeatTransferStiffnessKernel<nodes, gps, 2, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: HeatTransferMatrixKernel, Physics {
	HeatTransferStiffnessKernel(const HeatTransferMatrixKernel &base): HeatTransferMatrixKernel(base) {}

	void move(int n) { K += n * nodes * nodes; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.conductivity[gp][0], c01 = element.conductivity[gp][2];
			SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][3];
			SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD a = nx * c00 + ny * c01;
				SIMD b = nx * c10 + ny * c11;
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD res = load(out + (n * nodes + m) * SIMD::size);
					res = res + scale * (a * mx + b * my);
					store(out + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct HeatTransferStiffnessKernel<nodes, gps, 3, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: HeatTransferMatrixKernel, Physics {
	HeatTransferStiffnessKernel(const HeatTransferMatrixKernel &base): HeatTransferMatrixKernel(base) {}

	void move(int n) { K += n * nodes * nodes; }

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * load1(element.w[gp]);
			SIMD c00 = element.conductivity[gp][0], c01 = element.conductivity[gp][3], c02 = element.conductivity[gp][6];
			SIMD c10 = element.conductivity[gp][1], c11 = element.conductivity[gp][4], c12 = element.conductivity[gp][7];
			SIMD c20 = element.conductivity[gp][2], c21 = element.conductivity[gp][5], c22 = element.conductivity[gp][8];
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
					SIMD res = load(out + (n * nodes + m) * SIMD::size);
					res = res + scale * (a * mx + b * my + c * mz);
					store(out + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_MATRIX_H_ */
