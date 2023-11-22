
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_MATRIX_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_MATRIX_H_

#include "subkernel.h"

namespace espreso {

struct StructuralMechanicsMatrix: public SubKernel {
	const char* name() const { return "HeatTransferMatrix"; }

	double *K;
	Matrix_Shape shape;

	StructuralMechanicsMatrix()
	: K(nullptr), shape(Matrix_Shape::LOWER)
	{
		isconst = false;
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::ITERATION;
	}

	void activate(double *K)
	{
		this->K = K;
		this->isactive = 1;
	}

	inline void toFull(const size_t &size)
	{
		double * __restrict__ out = K;
		if (shape == Matrix_Shape::FULL) {
			for (size_t n = 0; n < size; ++n) {
				for (size_t m = n + 1; m < size; ++m) {
					SIMD res = load(out + (n * size + m) * SIMD::size);
					store(out + (m * size + n) * SIMD::size, res);
				}
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, enum Behaviour behaviour, enum ElasticityModel model, class Physics> struct StructuralMechanicsStiffness: StructuralMechanicsMatrix, Physics {
	StructuralMechanicsStiffness(const StructuralMechanicsMatrix &base): StructuralMechanicsMatrix(base) {}

	void simd(typename Physics::Element &element) {}
};

template <size_t nodes, size_t gps, enum ElasticityModel model, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 2, Behaviour::PLANE, model, Physics>: StructuralMechanicsMatrix, Physics {
	StructuralMechanicsStiffness(const StructuralMechanicsMatrix &base): StructuralMechanicsMatrix(base) {}

	// C
	// 0 1 _
	//   0 _
	//     2

	// B * C * Bt
	//
	// C = 3x3
	// B = dX  0 dY
	//      0 dY dX
	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		constexpr size_t size = 2 * nodes;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
			SIMD c0 = element.elasticity[gp][0];
			SIMD c1 = element.elasticity[gp][1];
			SIMD c2 = element.elasticity[gp][2];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = nx * c0;
					SIMD c = ny * c2;

					SIMD xx = load(out + (n * size + m) * SIMD::size);
					xx = xx + scale * (a * mx + c * my);
					store(out + (n * size + m) * SIMD::size, xx);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = nx * c1;
					SIMD c = ny * c2;

					SIMD xy = load(out + (n * size + m + nodes) * SIMD::size);
					xy = xy + scale * (b * my + c * mx);
					store(out + (n * size + m + nodes) * SIMD::size, xy);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = ny * c0;
					SIMD c = nx * c2;

					SIMD yy = load(out + ((n + nodes) * size + m + nodes) * SIMD::size);
					yy = yy + scale * (b * my + c * mx);
					store(out + ((n + nodes) * size + m + nodes) * SIMD::size, yy);
				}
			}
		}
		toFull(size);
		K += SIMD::size * 2 * nodes * 2 * nodes;
	}
};

template <size_t nodes, size_t gps, enum ElasticityModel model, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 2, Behaviour::AXISYMMETRIC, model, Physics>: StructuralMechanicsMatrix, Physics {
	StructuralMechanicsStiffness(const StructuralMechanicsMatrix &base): StructuralMechanicsMatrix(base) {}

	// C
	// 0 1 1 _
	//   0 1 _
	//     0 _
	//       2

	// B * C * Bt
	//
	// C = 4x4
	// B = dX  0  C dY
	//      0 dY  0 dX
	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		constexpr size_t size = 2 * nodes;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD coo[nodes];
			for (size_t n = 0; n < nodes; ++n) {
				coo[n] = load1(element.N[gp][n]) / element.gpcoords[gp][0];
			}
			SIMD scale = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0];
			SIMD c0 = element.elasticity[gp][0];
			SIMD c1 = element.elasticity[gp][1];
			SIMD c2 = element.elasticity[gp][2];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = nx * c0 + coo[n] * c1;
					SIMD c = nx * c1 + coo[n] * c0;
					SIMD d = ny * c2;

					SIMD xx = load(out + (n * size + m) * SIMD::size);
					xx = xx + scale * (a * mx + c * coo[m] + d * my);
					store(out + (n * size + m) * SIMD::size, xx);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = nx * c1 + coo[n] * c1;
					SIMD d = ny * c2;

					SIMD xy = load(out + (n * size + m + nodes) * SIMD::size);
					xy = xy + scale * (b * my + d * mx);
					store(out + (n * size + m + nodes) * SIMD::size, xy);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = ny * c0;
					SIMD d = nx * c2;

					SIMD yy = load(out + ((n + nodes) * size + m + nodes) * SIMD::size);
					yy = yy + scale * (b * my + d * mx);
					store(out + ((n + nodes) * size + m + nodes) * SIMD::size, yy);
				}
			}
		}
		toFull(2 * nodes);
		K += SIMD::size * 2 * nodes * 2 * nodes;
	}
};

template <size_t nodes, size_t gps, enum ElasticityModel model, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 3, Behaviour::VOLUME, model, Physics>: StructuralMechanicsMatrix, Physics {
	StructuralMechanicsStiffness(const StructuralMechanicsMatrix &base): StructuralMechanicsMatrix(base) {}

	// C
	// 0 1 1 _ _ _
	//   0 1 _ _ _
	//     0 _ _ _
	//       2 _ _
	//         2 _
	//           2

	// B * C * Bt
	//
	// C = 6x6
	// B = dX  0  0 dY  0 dZ
	//      0 dY  0 dX dZ  0
	//      0  0 dZ  0 dY dX
	//     - - - - - - - - -
	//      0  6 12 18 24 30
	//      a  b  c  d  e  f
	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		constexpr size_t size = 3 * nodes;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * load1(element.w[gp]);
			SIMD c0 = element.elasticity[gp][0];
			SIMD c1 = element.elasticity[gp][1];
			SIMD c2 = element.elasticity[gp][2];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD nz = element.dND[gp][n][2];
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD a = nx * c0;
					SIMD d = ny * c2;
					SIMD f = nz * c2;

					SIMD xx = load(out + (n * size + m) * SIMD::size);
					xx = xx + scale * (a * mx + d * my + f * mz);
					store(out + (n * size + m) * SIMD::size, xx);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD b = nx * c1;
					SIMD c = nx * c1;
					SIMD d = ny * c2;
					SIMD f = nz * c2;

					SIMD xy = load(out + (n * size + m + nodes) * SIMD::size);
					SIMD xz = load(out + (n * size + m + nodes * 2) * SIMD::size);
					xy = xy + scale * (b * my + d * mx);
					xz = xz + scale * (c * mz + f * mx);
					store(out + (n * size + m + nodes) * SIMD::size, xy);
					store(out + (n * size + m + nodes * 2) * SIMD::size, xz);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD b = ny * c0;
					SIMD d = nx * c2;
					SIMD e = nz * c2;

					SIMD yy = load(out + ((n + nodes) * size + m + nodes) * SIMD::size);
					yy = yy + scale * (b * my + d * mx + e * mz);
					store(out + ((n + nodes) * size + m + nodes) * SIMD::size, yy);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD c = ny * c1;
					SIMD e = nz * c2;

					SIMD yz = load(out + ((n + nodes) * size + m + nodes * 2) * SIMD::size);
					yz = yz + scale * (c * mz + e * my);
					store(out + ((n + nodes) * size + m + nodes * 2) * SIMD::size, yz);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD c = nz * c0;
					SIMD e = ny * c2;
					SIMD f = nx * c2;

					SIMD zz = load(out + ((n + nodes * 2) * size + m + nodes * 2) * SIMD::size);
					zz = zz + scale * (c * mz + e * my + f * mx);
					store(out + ((n + nodes * 2) * size + m + nodes * 2) * SIMD::size, zz);
				}
			}
		}
		toFull(3 * nodes);
		K += SIMD::size * 3 * nodes * 3 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 3, Behaviour::VOLUME, ElasticityModel::ORTHOTROPIC, Physics>: StructuralMechanicsMatrix, Physics {
	StructuralMechanicsStiffness(const StructuralMechanicsMatrix &base): StructuralMechanicsMatrix(base) {}

	// C
	// 0 1 2 _ _ _
	//   3 4 _ _ _
	//     5 _ _ _
	//       6 _ _
	//         7 _
	//           8

	// B * C * Bt
	//
	// C = 6x6
	// B = dX  0  0 dY  0 dZ
	//      0 dY  0 dX dZ  0
	//      0  0 dZ  0 dY dX
	//     - - - - - - - - -
	//      0  6 12 18 24 30
	//      a  b  c  d  e  f
	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		constexpr size_t size = 3 * nodes;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * load1(element.w[gp]);
			SIMD c00 = element.elasticity[gp][0];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][3];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][4], c22 = element.elasticity[gp][5];
			SIMD c33 = element.elasticity[gp][6];
			SIMD c44 = element.elasticity[gp][7];
			SIMD c55 = element.elasticity[gp][8];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD nz = element.dND[gp][n][2];
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD a = nx * c00;
					SIMD d = ny * c33;
					SIMD f = nz * c55;

					SIMD xx = load(out + (n * size + m) * SIMD::size);
					xx = xx + scale * (a * mx + d * my + f * mz);
					store(out + (n * size + m) * SIMD::size, xx);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD b = nx * c01;
					SIMD c = nx * c02;
					SIMD d = ny * c33;
					SIMD f = nz * c55;

					SIMD xy = load(out + (n * size + m + nodes) * SIMD::size);
					SIMD xz = load(out + (n * size + m + nodes * 2) * SIMD::size);
					xy = xy + scale * (b * my + d * mx);
					xz = xz + scale * (c * mz + f * mx);
					store(out + (n * size + m + nodes) * SIMD::size, xy);
					store(out + (n * size + m + nodes * 2) * SIMD::size, xz);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD b = ny * c11;
					SIMD d = nx * c33;
					SIMD e = nz * c44;

					SIMD yy = load(out + ((n + nodes) * size + m + nodes) * SIMD::size);
					yy = yy + scale * (b * my + d * mx + e * mz);
					store(out + ((n + nodes) * size + m + nodes) * SIMD::size, yy);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD c = ny * c12;
					SIMD e = nz * c44;

					SIMD yz = load(out + ((n + nodes) * size + m + nodes * 2) * SIMD::size);
					yz = yz + scale * (c * mz + e * my);
					store(out + ((n + nodes) * size + m + nodes * 2) * SIMD::size, yz);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD c = nz * c22;
					SIMD e = ny * c44;
					SIMD f = nx * c55;

					SIMD zz = load(out + ((n + nodes * 2) * size + m + nodes * 2) * SIMD::size);
					zz = zz + scale * (c * mz + e * my + f * mx);
					store(out + ((n + nodes * 2) * size + m + nodes * 2) * SIMD::size, zz);
				}
			}
		}
		toFull(3 * nodes);
		K += SIMD::size * 3 * nodes * 3 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 3, Behaviour::VOLUME, ElasticityModel::SYMMETRIC, Physics>: StructuralMechanicsMatrix, Physics {
	StructuralMechanicsStiffness(const StructuralMechanicsMatrix &base): StructuralMechanicsMatrix(base) {}

	// B * C * Bt
	//
	// C = 6x6
	// B = dX  0  0 dY  0 dZ
	//      0 dY  0 dX dZ  0
	//      0  0 dZ  0 dY dX
	//     - - - - - - - - -
	//      0  6 12 18 24 30
	//      a  b  c  d  e  f
	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = K;
		constexpr size_t size = 3 * nodes;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.det[gp] * load1(element.w[gp]);
			SIMD c00 = element.elasticity[gp][0];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][ 6];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][ 7], c22 = element.elasticity[gp][11];
			SIMD c03 = element.elasticity[gp][3], c13 = element.elasticity[gp][ 8], c23 = element.elasticity[gp][12], c33 = element.elasticity[gp][15];
			SIMD c04 = element.elasticity[gp][4], c14 = element.elasticity[gp][ 9], c24 = element.elasticity[gp][13], c34 = element.elasticity[gp][16], c44 = element.elasticity[gp][18];
			SIMD c05 = element.elasticity[gp][5], c15 = element.elasticity[gp][10], c25 = element.elasticity[gp][14], c35 = element.elasticity[gp][17], c45 = element.elasticity[gp][19], c55 = element.elasticity[gp][20];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD nz = element.dND[gp][n][2];
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD a = nx * c00 + ny * c03 + nz * c05;
					SIMD d = nx * c03 + ny * c33 + nz * c35;
					SIMD f = nx * c05 + ny * c35 + nz * c55;

					SIMD xx = load(out + (n * size + m) * SIMD::size);
					xx = xx + scale * (a * mx + d * my + f * mz);
					store(out + (n * size + m) * SIMD::size, xx);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD b = nx * c01 + ny * c13 + nz * c15;
					SIMD c = nx * c02 + ny * c23 + nz * c25;
					SIMD d = nx * c03 + ny * c33 + nz * c35;
					SIMD e = nx * c04 + ny * c34 + nz * c45;
					SIMD f = nx * c05 + ny * c35 + nz * c55;

					SIMD xy = load(out + (n * size + m + nodes) * SIMD::size);
					SIMD xz = load(out + (n * size + m + nodes * 2) * SIMD::size);
					xy = xy + scale * (b * my + d * mx + e * mz);
					xz = xz + scale * (c * mz + e * my + f * mx);
					store(out + (n * size + m + nodes) * SIMD::size, xy);
					store(out + (n * size + m + nodes * 2) * SIMD::size, xz);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD b = ny * c11 + nx * c13 + nz * c14;
					SIMD d = ny * c13 + nx * c33 + nz * c34;
					SIMD e = ny * c14 + nx * c34 + nz * c44;

					SIMD yy = load(out + ((n + nodes) * size + m + nodes) * SIMD::size);
					yy = yy + scale * (b * my + d * mx + e * mz);
					store(out + ((n + nodes) * size + m + nodes) * SIMD::size, yy);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD c = ny * c12 + nx * c23 + nz * c24;
					SIMD e = ny * c14 + nx * c34 + nz * c44;
					SIMD f = ny * c15 + nx * c35 + nz * c45;

					SIMD yz = load(out + ((n + nodes) * size + m + nodes * 2) * SIMD::size);
					yz = yz + scale * (c * mz + e * my + f * mx);
					store(out + ((n + nodes) * size + m + nodes * 2) * SIMD::size, yz);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD c = nz * c22 + ny * c24 + nx * c25;
					SIMD e = nz * c24 + ny * c44 + nx * c45;
					SIMD f = nz * c25 + ny * c45 + nx * c55;

					SIMD zz = load(out + ((n + nodes * 2) * size + m + nodes * 2) * SIMD::size);
					zz = zz + scale * (c * mz + e * my + f * mx);
					store(out + ((n + nodes * 2) * size + m + nodes * 2) * SIMD::size, zz);
				}
			}
		}
		toFull(3 * nodes);
		K += SIMD::size * 3 * nodes * 3 * nodes;
	}
};

template <size_t nodes, size_t gps, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 3, Behaviour::VOLUME, ElasticityModel::ANISOTROPIC, Physics>: StructuralMechanicsMatrix, Physics {
	StructuralMechanicsStiffness(const StructuralMechanicsMatrix &base): StructuralMechanicsMatrix(base) {}

	// B * C * Bt
	//
	// C = 6x6
	// B = dX  0  0 dY  0 dZ
	//      0 dY  0 dX dZ  0
	//      0  0 dZ  0 dY dX
	//     - - - - - - - - -
	//      0  6 12 18 24 30
	//      a  b  c  d  e  f
	void simd(typename Physics::Element &element)
	{
//		double * __restrict__ out = K;
//		constexpr size_t size = 3 * nodes;
//		for (size_t gp = 0; gp < gps; ++gp) {
//			SIMD scale = element.det[gp] * load1(element.w[gp]);
//			SIMD c00 = element.elasticity[gp][0], c10 = element.elasticity[gp][ 6], c20 = element.elasticity[gp][12], c30 = element.elasticity[gp][18], c40 = element.elasticity[gp][24], c50 = element.elasticity[gp][30];
//			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][ 7], c21 = element.elasticity[gp][13], c31 = element.elasticity[gp][19], c41 = element.elasticity[gp][25], c51 = element.elasticity[gp][31];
//			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][ 8], c22 = element.elasticity[gp][14], c32 = element.elasticity[gp][20], c42 = element.elasticity[gp][26], c52 = element.elasticity[gp][32];
//			SIMD c03 = element.elasticity[gp][3], c13 = element.elasticity[gp][ 9], c23 = element.elasticity[gp][15], c33 = element.elasticity[gp][21], c43 = element.elasticity[gp][27], c53 = element.elasticity[gp][33];
//			SIMD c04 = element.elasticity[gp][4], c14 = element.elasticity[gp][10], c24 = element.elasticity[gp][16], c34 = element.elasticity[gp][22], c44 = element.elasticity[gp][28], c54 = element.elasticity[gp][34];
//			SIMD c05 = element.elasticity[gp][5], c15 = element.elasticity[gp][11], c25 = element.elasticity[gp][17], c35 = element.elasticity[gp][23], c45 = element.elasticity[gp][29], c55 = element.elasticity[gp][35];
//			for (size_t n = 0; n < nodes; ++n) {
//				SIMD nx = element.dND[gp][n][0];
//				SIMD ny = element.dND[gp][n][1];
//				SIMD nz = element.dND[gp][n][2];
//				for (size_t m = 0; m < nodes; ++m) {
//					SIMD mx = element.dND[gp][m][0];
//					SIMD my = element.dND[gp][m][1];
//					SIMD mz = element.dND[gp][m][2];
//					SIMD a = nx * c00 + ny * c30 + nz * c50;
//					SIMD d = nx * c03 + ny * c33 + nz * c53;
//					SIMD f = nx * c05 + ny * c35 + nz * c55;
//
//					SIMD xx = load(out + (n * size + m) * SIMD::size);
//					xx = xx + scale * (a * mx + d * my + f * mz);
//					store(out + (n * size + m) * SIMD::size, xx);
//				}
//				for (size_t m = 0; m < nodes; ++m) {
//					SIMD mx = element.dND[gp][m][0];
//					SIMD my = element.dND[gp][m][1];
//					SIMD mz = element.dND[gp][m][2];
//					SIMD b = nx * c01 + ny * c31 + nz * c51;
//					SIMD c = nx * c02 + ny * c32 + nz * c52;
//					SIMD d = nx * c03 + ny * c33 + nz * c53;
//					SIMD e = nx * c04 + ny * c34 + nz * c54;
//					SIMD f = nx * c05 + ny * c35 + nz * c55;
//
//					SIMD xy = load(out + (n * size + m + nodes) * SIMD::size);
//					SIMD xz = load(out + (n * size + m + nodes * 2) * SIMD::size);
//					xy = xy + scale * (b * my + d * mx);
//					xz = xz + scale * (c * mz + f * mx);
//					store(out + (n * size + m + nodes) * SIMD::size, xy);
//					store(out + (n * size + m + nodes * 2) * SIMD::size, xz);
//				}
//				for (size_t m = n; m < nodes; ++m) {
//					SIMD mx = element.dND[gp][m][0];
//					SIMD my = element.dND[gp][m][1];
//					SIMD mz = element.dND[gp][m][2];
//					SIMD b = ny * c11 + nx * c31 + nz * c41;
//					SIMD d = ny * c13 + nx * c33 + nz * c43;
//					SIMD e = ny * c14 + nx * c34 + nz * c44;
//
//					SIMD yy = load(out + ((n + nodes) * size + m + nodes) * SIMD::size);
//					yy = yy + scale * (b * my + d * mx + e * mz);
//					store(out + ((n + nodes) * size + m + nodes) * SIMD::size, yy);
//				}
//				for (size_t m = 0; m < nodes; ++m) {
//					SIMD mx = element.dND[gp][m][0];
//					SIMD my = element.dND[gp][m][1];
//					SIMD mz = element.dND[gp][m][2];
//					SIMD c = ny * c12 + nx * c32 + nz * c42;
//					SIMD e = ny * c14 + nx * c34 + nz * c44;
//					SIMD f = ny * c15 + nx * c35 + nz * c45;
//
//					SIMD yz = load(out + ((n + nodes) * size + m + nodes * 2) * SIMD::size);
//					yz = yz + scale * (c * mz + e * my);
//					store(out + ((n + nodes) * size + m + nodes * 2) * SIMD::size, yz);
//				}
//				for (size_t m = n; m < nodes; ++m) {
//					SIMD mx = element.dND[gp][m][0];
//					SIMD my = element.dND[gp][m][1];
//					SIMD mz = element.dND[gp][m][2];
//					SIMD a = nz * c20 + ny * c40 + nx * c50;
//					SIMD b = nz * c21 + ny * c41 + nx * c51;
//					SIMD c = nz * c22 + ny * c42 + nx * c52;
//					SIMD d = nz * c23 + ny * c43 + nx * c53;
//					SIMD e = nz * c24 + ny * c44 + nx * c54;
//					SIMD f = nz * c25 + ny * c45 + nx * c55;
//
//					SIMD zz = load(out + ((n + nodes * 2) * size + m + nodes * 2) * SIMD::size);
//					zz = zz + scale * (c * mz + e * my + f * mx);
//					store(out + ((n + nodes * 2) * size + m + nodes * 2) * SIMD::size, zz);
//				}
//			}
//		}
//		K += SIMD::size * 3 * nodes * 3 * nodes;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_MATRIX_H_ */
