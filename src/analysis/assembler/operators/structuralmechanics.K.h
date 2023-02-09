
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_K_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_K_H_

#include "analysis/assembler/operator.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct HeatTransferStiffness;

struct StructuralMechanicsStiffnessBase: public ActionOperator {
	OutputParameterIterator stiffness;

	StructuralMechanicsStiffnessBase(size_t interval, ParameterData &stiffness)
	: stiffness(stiffness, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		stiffness += n;
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
			SIMD c00 = element.elasticity[gp][0], c10 = element.elasticity[gp][3], c20 = element.elasticity[gp][6];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][4], c21 = element.elasticity[gp][7];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][5], c22 = element.elasticity[gp][8];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = nx * c00 + ny * c20;
					SIMD b = nx * c01 + ny * c21;
					SIMD c = nx * c02 + ny * c22;

					SIMD xx = load(out + (2 * n * nodes + m        ) * SIMD::size);
					SIMD xy = load(out + (2 * n * nodes + m + nodes) * SIMD::size);
					xx = xx + scale * (a * mx + c * my);
					xy = xy + scale * (b * my + c * mx);
					store(out + (2 * n * nodes + m        ) * SIMD::size, xx);
					store(out + (2 * n * nodes + m + nodes) * SIMD::size, xy);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = ny * c10 + nx * c20;
					SIMD b = ny * c11 + nx * c21;
					SIMD c = ny * c12 + nx * c22;

					SIMD yx = load(out + (2 * nodes * nodes + 2 * n * nodes + m        ) * SIMD::size);
					SIMD yy = load(out + (2 * nodes * nodes + 2 * n * nodes + m + nodes) * SIMD::size);
					yx = yx + scale * (a * mx + c * my);
					yy = yy + scale * (b * my + c * mx);
					store(out + (2 * nodes * nodes + 2 * n * nodes + m        ) * SIMD::size, yx);
					store(out + (2 * nodes * nodes + 2 * n * nodes + m + nodes) * SIMD::size, yy);
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
			SIMD c00 = element.elasticity[gp][0], c10 = element.elasticity[gp][4], c20 = element.elasticity[gp][ 8], c30 = element.elasticity[gp][12];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][5], c21 = element.elasticity[gp][ 9], c31 = element.elasticity[gp][13];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][6], c22 = element.elasticity[gp][10], c32 = element.elasticity[gp][14];
			SIMD c03 = element.elasticity[gp][3], c13 = element.elasticity[gp][7], c23 = element.elasticity[gp][11], c33 = element.elasticity[gp][15];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = nx * c00 + coo[n] * c20 + ny * c30;
					SIMD b = nx * c01 + coo[n] * c21 + ny * c31;
					SIMD c = nx * c02 + coo[n] * c22 + ny * c32;
					SIMD d = nx * c03 + coo[n] * c23 + ny * c33;

					SIMD xx = load(out + (2 * n * nodes + m        ) * SIMD::size);
					SIMD xy = load(out + (2 * n * nodes + m + nodes) * SIMD::size);
					xx = xx + scale * (a * mx + c * coo[m] + d * my);
					xy = xy + scale * (b * my              + d * mx);
					store(out + (2 * n * nodes + m        ) * SIMD::size, xx);
					store(out + (2 * n * nodes + m + nodes) * SIMD::size, xy);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD a = ny * c10 + nx * c30;
					SIMD b = ny * c11 + nx * c31;
					SIMD c = ny * c12 + nx * c32;
					SIMD d = ny * c13 + nx * c33;

					SIMD yx = load(out + (2 * nodes * nodes + 2 * n * nodes + m        ) * SIMD::size);
					SIMD yy = load(out + (2 * nodes * nodes + 2 * n * nodes + m + nodes) * SIMD::size);
					yx = yx + scale * (a * mx + c * coo[m] + d * my);
					yy = yy + scale * (b * my              + d * mx);
					store(out + (2 * nodes * nodes + 2 * n * nodes + m        ) * SIMD::size, yx);
					store(out + (2 * nodes * nodes + 2 * n * nodes + m + nodes) * SIMD::size, yy);
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
			SIMD c00 = element.elasticity[gp][0], c10 = element.elasticity[gp][ 6], c20 = element.elasticity[gp][12]; SIMD c30 = element.elasticity[gp][18], c40 = element.elasticity[gp][24], c50 = element.elasticity[gp][30];
			SIMD c01 = element.elasticity[gp][1], c11 = element.elasticity[gp][ 7], c21 = element.elasticity[gp][13]; SIMD c31 = element.elasticity[gp][19], c41 = element.elasticity[gp][25], c51 = element.elasticity[gp][31];
			SIMD c02 = element.elasticity[gp][2], c12 = element.elasticity[gp][ 8], c22 = element.elasticity[gp][14]; SIMD c32 = element.elasticity[gp][20], c42 = element.elasticity[gp][26], c52 = element.elasticity[gp][32];
			SIMD c03 = element.elasticity[gp][3], c13 = element.elasticity[gp][ 9], c23 = element.elasticity[gp][15]; SIMD c33 = element.elasticity[gp][21], c43 = element.elasticity[gp][27], c53 = element.elasticity[gp][33];
			SIMD c04 = element.elasticity[gp][4], c14 = element.elasticity[gp][10], c24 = element.elasticity[gp][16]; SIMD c34 = element.elasticity[gp][22], c44 = element.elasticity[gp][28], c54 = element.elasticity[gp][34];
			SIMD c05 = element.elasticity[gp][5], c15 = element.elasticity[gp][11], c25 = element.elasticity[gp][17]; SIMD c35 = element.elasticity[gp][23], c45 = element.elasticity[gp][29], c55 = element.elasticity[gp][35];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gp][n][0];
				SIMD ny = element.dND[gp][n][1];
				SIMD nz = element.dND[gp][n][2];
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD a = nx * c00 + ny * c30 + nz * c50;
					SIMD b = nx * c01 + ny * c31 + nz * c51;
					SIMD c = nx * c02 + ny * c32 + nz * c52;
					SIMD d = nx * c03 + ny * c33 + nz * c53;
					SIMD e = nx * c04 + ny * c34 + nz * c54;
					SIMD f = nx * c05 + ny * c35 + nz * c55;

					SIMD xx = load(out + (3 * n * nodes + 0 * nodes + m) * SIMD::size);
					SIMD xy = load(out + (3 * n * nodes + 1 * nodes + m) * SIMD::size);
					SIMD xz = load(out + (3 * n * nodes + 2 * nodes + m) * SIMD::size);
					xx = xx + scale * (a * mx + d * my + f * mz);
					xy = xy + scale * (b * my + d * mx + e * mz);
					xz = xz + scale * (c * mz + e * my + f * mx);
					store(out + (3 * n * nodes + 0 * nodes + m) * SIMD::size, xx);
					store(out + (3 * n * nodes + 1 * nodes + m) * SIMD::size, xy);
					store(out + (3 * n * nodes + 2 * nodes + m) * SIMD::size, xz);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD a = ny * c10 + nx * c30 + nz * c40;
					SIMD b = ny * c11 + nx * c31 + nz * c41;
					SIMD c = ny * c12 + nx * c32 + nz * c42;
					SIMD d = ny * c13 + nx * c33 + nz * c43;
					SIMD e = ny * c14 + nx * c34 + nz * c44;
					SIMD f = ny * c15 + nx * c35 + nz * c45;

					SIMD yx = load(out + (3 * nodes * nodes + 3 * n * nodes + 0 * nodes + m) * SIMD::size);
					SIMD yy = load(out + (3 * nodes * nodes + 3 * n * nodes + 1 * nodes + m) * SIMD::size);
					SIMD yz = load(out + (3 * nodes * nodes + 3 * n * nodes + 2 * nodes + m) * SIMD::size);
					yx = yx + scale * (a * mx + d * my + f * mz);
					yy = yy + scale * (b * my + d * mx + e * mz);
					yz = yz + scale * (c * mz + e * my + f * mx);
					store(out + (3 * nodes * nodes + 3 * n * nodes + 0 * nodes + m) * SIMD::size, yx);
					store(out + (3 * nodes * nodes + 3 * n * nodes + 1 * nodes + m) * SIMD::size, yy);
					store(out + (3 * nodes * nodes + 3 * n * nodes + 2 * nodes + m) * SIMD::size, yz);
				}
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gp][m][0];
					SIMD my = element.dND[gp][m][1];
					SIMD mz = element.dND[gp][m][2];
					SIMD a = nz * c20 + ny * c40 + nx * c50;
					SIMD b = nz * c21 + ny * c41 + nx * c51;
					SIMD c = nz * c22 + ny * c42 + nx * c52;
					SIMD d = nz * c23 + ny * c43 + nx * c53;
					SIMD e = nz * c24 + ny * c44 + nx * c54;
					SIMD f = nz * c25 + ny * c45 + nx * c55;

					SIMD zx = load(out + (6 * nodes * nodes + 3 * n * nodes + 0 * nodes + m) * SIMD::size);
					SIMD zy = load(out + (6 * nodes * nodes + 3 * n * nodes + 1 * nodes + m) * SIMD::size);
					SIMD zz = load(out + (6 * nodes * nodes + 3 * n * nodes + 2 * nodes + m) * SIMD::size);
					zx = zx + scale * (a * mx + d * my + f * mz);
					zy = zy + scale * (b * my + d * mx + e * mz);
					zz = zz + scale * (c * mz + e * my + f * mx);
					store(out + (6 * nodes * nodes + 3 * n * nodes + 0 * nodes + m) * SIMD::size, zx);
					store(out + (6 * nodes * nodes + 3 * n * nodes + 1 * nodes + m) * SIMD::size, zy);
					store(out + (6 * nodes * nodes + 3 * n * nodes + 2 * nodes + m) * SIMD::size, zz);
				}
			}
		}
		move(SIMD::size);
	}
};

//struct StructuralMechanicsStiffness: public ActionOperator {
//	StructuralMechanicsStiffness(
//			int interval,
//			const ParameterData &N,
//			const ParameterData &dND,
//			const ParameterData &weight,
//			const ParameterData &determinant,
//			const ParameterData &coordinates,
//			const ParameterData &elasticity,
//			const ParameterData &thickness,
//			ParameterData &stiffness)
//	: N(N, interval),
//	  dND(dND, interval),
//	  weight(weight, interval, 0),
//	  determinant(determinant, interval),
//	  coordinates(coordinates, interval),
//	  elasticity(elasticity, interval),
//	  thickness(thickness, interval),
//	  stiffness(stiffness, interval)
//	{
//
//	}
//
//	InputParameterIterator N, dND, weight, determinant, coordinates, elasticity, thickness;
//	OutputParameterIterator stiffness;
//
//	void operator++()
//	{
//		++dND; ++determinant; ++coordinates; ++elasticity; ++thickness;
//		++stiffness;
//	}
//
//	void move(int n)
//	{
//		dND += n; determinant += n; coordinates += n; elasticity += n; thickness += n;
//		stiffness += n;
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct StiffnessPlane: public StructuralMechanicsStiffness {
//	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;
//
//	void operator()()
//	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDDXDY33DXDYN<nodes>(determinant[gpindex] * weight[gpindex], elasticity.data + 9 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
//		}
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct StiffnessPlaneWithThickness: public StructuralMechanicsStiffness {
//	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;
//
//	void operator()()
//	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDDXDY33DXDYN<nodes>(determinant[gpindex] * weight[gpindex] * thickness[gpindex], elasticity.data + 9 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
//		}
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct StiffnessAxisymmetric: public StructuralMechanicsStiffness {
//	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;
//
//	void operator()()
//	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			double coo[nodes];
//			for (int n = 0; n < nodes; ++n) {
//				coo[n] = N[gpindex * nodes + n] / coordinates[2 * gpindex];
//			}
//			ADDDXDYCOO44DXDYN<nodes>(determinant[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex], elasticity.data + 16 * gpindex, dND.data + 2 * nodes * gpindex, coo, stiffness.data);
//		}
//	}
//};
//
//
//template<size_t nodes, size_t gps>
//struct Stiffness3DElasticity: public StructuralMechanicsStiffness {
//	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;
//
//	void operator()()
//	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDDXDYDZ66DXDYDZN<nodes>(determinant[gpindex] * weight[gpindex], elasticity.data + 36 * gpindex, dND.data + 3 * nodes * gpindex, stiffness.data);
//		}
//	}
//};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_K_H_ */
