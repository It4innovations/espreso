
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_K_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_K_H_

#include <analysis/assembler/subkernel/operator.h>

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct HeatTransferStiffness;

struct StructuralMechanicsStiffnessBase: public ActionOperator {
	const char* name() const { return "StructuralMechanicsStiffnessBase"; }

	OutputParameterIterator stiffness;

	StructuralMechanicsStiffnessBase(size_t interval, ParameterData &stiffness)
	: stiffness(stiffness, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void move(int n)
	{
		stiffness += n;
	}

	inline size_t index(const size_t &r, const size_t &c, const size_t &size)
	{
		return r * size + c - ((r + 1) * r / 2);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct StructuralMechanicsStiffness;

template <size_t nodes, size_t gps, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 2, 2, StructuralMechanicsElementType::SYMMETRIC_PLANE, Physics>: StructuralMechanicsStiffnessBase, Physics {
	using StructuralMechanicsStiffnessBase::StructuralMechanicsStiffnessBase;

	// B * C * Bt
	//
	// C = 3x3
	// B = dX  0 dY
	//      0 dY dX
	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = stiffness.data;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]);
			SIMD c00 = element.elasticity[gp][0];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][3];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][4], c22 = element.elasticity[gp][5];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = nx * c00 + ny * c02;
					SIMD c = nx * c02 + ny * c22;

					size_t i = index(n, m, 2 * nodes);
					SIMD xx = load(out + i * SIMD::size);
					xx = xx + scale * (a * mx + c * my);
					store(out + i * SIMD::size, xx);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = nx * c01 + ny * c12;
					SIMD c = nx * c02 + ny * c22;

					size_t i = index(n, m + nodes, 2 * nodes);
					SIMD xy = load(out + i * SIMD::size);
					xy = xy + scale * (b * my + c * mx);
					store(out + i * SIMD::size, xy);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = ny * c11 + nx * c12;
					SIMD c = ny * c12 + nx * c22;

					size_t i = index(n + nodes, m + nodes, 2 * nodes);
					SIMD yy = load(out + i * SIMD::size);
					yy = yy + scale * (b * my + c * mx);
					store(out + i * SIMD::size, yy);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 2, 2, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC, Physics>: StructuralMechanicsStiffnessBase, Physics {
	using StructuralMechanicsStiffnessBase::StructuralMechanicsStiffnessBase;

	// B * C * Bt
	//
	// C = 4x4
	// B = dX  0  C dY
	//      0 dY  0 dX
	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = stiffness.data;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD coo[nodes];
			for (size_t n = 0; n < nodes; ++n) {
				coo[n] = load1(element.N[gp][n]) / element.gpcoords[gp][0];
			}
			SIMD scale = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0];
			SIMD c00 = element.elasticity[gp][0];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][4];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][5], c22 = element.elasticity[gp][7];
			SIMD c03 = element.elasticity[gp][3], c13 = element.elasticity[gp][6], c23 = element.elasticity[gp][8], c33 = element.elasticity[gp][9];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = nx * c00 + coo[n] * c02 + ny * c03;
					SIMD c = nx * c02 + coo[n] * c22 + ny * c23;
					SIMD d = nx * c03 + coo[n] * c23 + ny * c33;

					size_t i = index(n, m, 2 * nodes);
					SIMD xx = load(out + i * SIMD::size);
					xx = xx + scale * (a * mx + c * coo[m] + d * my);
					store(out + i * SIMD::size, xx);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = nx * c01 + coo[n] * c12 + ny * c13;
					SIMD d = nx * c03 + coo[n] * c23 + ny * c33;

					size_t i = index(n, m + nodes, 2 * nodes);
					SIMD xy = load(out + i * SIMD::size);
					xy = xy + scale * (b * my + d * mx);
					store(out + i * SIMD::size, xy);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD b = ny * c11 + nx * c13;
					SIMD d = ny * c13 + nx * c33;

					size_t i = index(n + nodes, m + nodes, 2 * nodes);
					SIMD yy = load(out + i * SIMD::size);
					yy = yy + scale * (b * my + d * mx);
					store(out + i * SIMD::size, yy);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct StructuralMechanicsStiffness<nodes, gps, 3, 3, StructuralMechanicsElementType::SYMMETRIC_VOLUME, Physics>: StructuralMechanicsStiffnessBase, Physics {
	using StructuralMechanicsStiffnessBase::StructuralMechanicsStiffnessBase;

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
		double * __restrict__ out = stiffness.data;
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

					size_t i = index(n, m, 3 * nodes);
					SIMD xx = load(out + i * SIMD::size);
					xx = xx + scale * (a * mx + d * my + f * mz);
					store(out + i * SIMD::size, xx);
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

					size_t iy = index(n, m + nodes, 3 * nodes);
					size_t iz = index(n, m + nodes * 2, 3 * nodes);
					SIMD xy = load(out + iy * SIMD::size);
					SIMD xz = load(out + iz * SIMD::size);
					xy = xy + scale * (b * my + d * mx + e * mz);
					xz = xz + scale * (c * mz + e * my + f * mx);
					store(out + iy * SIMD::size, xy);
					store(out + iz * SIMD::size, xz);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD b = ny * c11 + nx * c13 + nz * c14;
					SIMD d = ny * c13 + nx * c33 + nz * c34;
					SIMD e = ny * c14 + nx * c34 + nz * c44;

					size_t i = index(n + nodes, m + nodes, 3 * nodes);
					SIMD yy = load(out + i * SIMD::size);
					yy = yy + scale * (b * my + d * mx + e * mz);
					store(out + i * SIMD::size, yy);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD c = ny * c12 + nx * c23 + nz * c24;
					SIMD e = ny * c14 + nx * c34 + nz * c44;
					SIMD f = ny * c15 + nx * c35 + nz * c45;

					size_t i = index(n + nodes, m + nodes * 2, 3 * nodes);
					SIMD yz = load(out + i * SIMD::size);
					yz = yz + scale * (c * mz + e * my + f * mx);
					store(out + i * SIMD::size, yz);
				}
				for (size_t m = n; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD a = nz * c02 + ny * c04 + nx * c05;
					SIMD b = nz * c12 + ny * c14 + nx * c15;
					SIMD c = nz * c22 + ny * c24 + nx * c25;
					SIMD d = nz * c23 + ny * c34 + nx * c35;
					SIMD e = nz * c24 + ny * c44 + nx * c45;
					SIMD f = nz * c25 + ny * c45 + nx * c55;

					size_t i = index(n + nodes * 2, m + nodes * 2, 3 * nodes);
					SIMD zz = load(out + i * SIMD::size);
					zz = zz + scale * (c * mz + e * my + f * mx);
					store(out + i * SIMD::size, zz);
				}
			}
		}
		move(SIMD::size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_K_H_ */
