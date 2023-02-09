
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_

#include "analysis/assembler/operator.h"
#include "math/simd/simd.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneStrain: ActionOperator, Physics {

	ElasticityIsotropicPlaneStrain(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD k = ex * (C1 - mi) / ((C1 + mi) * (C1 - C2 * mi));
			element.elasticity[gp][0] = k;
			element.elasticity[gp][1] = k * (mi / (C1 - mi));
			element.elasticity[gp][2] = zeros();
			element.elasticity[gp][3] = k * (mi / (C1 - mi));
			element.elasticity[gp][4] = k;
			element.elasticity[gp][5] = zeros();
			element.elasticity[gp][6] = zeros();
			element.elasticity[gp][7] = zeros();
			element.elasticity[gp][8] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneStress: ActionOperator, Physics {

	ElasticityIsotropicPlaneStress(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C1 = load1(1.);
			SIMD C2 = load1(.5);
			SIMD k = ex / (C1 - mi * mi);
			element.elasticity[gp][0] = k;
			element.elasticity[gp][1] = k * mi;
			element.elasticity[gp][2] = zeros();
			element.elasticity[gp][3] = k * mi;
			element.elasticity[gp][4] = k;
			element.elasticity[gp][5] = zeros();
			element.elasticity[gp][6] = zeros();
			element.elasticity[gp][7] = zeros();
			element.elasticity[gp][8] = k * ((C1 -  mi) * C2);
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicPlaneAxisymmetric: ActionOperator, Physics {

	ElasticityIsotropicPlaneAxisymmetric(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD k = ex * (C1 - mi) / ((C1 + mi) * (C1 - C2 * mi));
			element.elasticity[gp][ 0] = k;
			element.elasticity[gp][ 1] = k * (mi / (C1 - mi));
			element.elasticity[gp][ 2] = zeros();
			element.elasticity[gp][ 3] = k * (mi / (C1 - mi));
			element.elasticity[gp][ 4] = k * (mi / (C1 - mi));
			element.elasticity[gp][ 5] = k;
			element.elasticity[gp][ 6] = zeros();
			element.elasticity[gp][ 7] = k * (mi / (C1 - mi));
			element.elasticity[gp][ 8] = zeros();
			element.elasticity[gp][ 9] = zeros();
			element.elasticity[gp][10] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
			element.elasticity[gp][11] = zeros();
			element.elasticity[gp][12] = k * (mi / (C1 - mi));
			element.elasticity[gp][13] = k * (mi / (C1 - mi));
			element.elasticity[gp][14] = zeros();
			element.elasticity[gp][15] = k;
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct ElasticityIsotropicVolume: ActionOperator, Physics {

	ElasticityIsotropicVolume(size_t interval)
	{
		action = Action::ASSEMBLE;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C05 = load1(.5);
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD ee = ex / ((C1 + mi) * (C1 - C2 * mi));
			element.elasticity[gp][ 0] = ee * (C1 - mi);
			element.elasticity[gp][ 1] = ee * mi;
			element.elasticity[gp][ 2] = ee * mi;
			element.elasticity[gp][ 3] = zeros();
			element.elasticity[gp][ 4] = zeros();
			element.elasticity[gp][ 5] = zeros();
			element.elasticity[gp][ 6] = ee * mi;
			element.elasticity[gp][ 7] = ee * (C1 - mi);
			element.elasticity[gp][ 8] = ee * mi;
			element.elasticity[gp][ 9] = zeros();
			element.elasticity[gp][10] = zeros();
			element.elasticity[gp][11] = zeros();
			element.elasticity[gp][12] = ee * mi;
			element.elasticity[gp][13] = ee * mi;
			element.elasticity[gp][14] = ee * (C1 - mi);
			element.elasticity[gp][15] = zeros();
			element.elasticity[gp][16] = zeros();
			element.elasticity[gp][17] = zeros();
			element.elasticity[gp][18] = zeros();
			element.elasticity[gp][19] = zeros();
			element.elasticity[gp][20] = zeros();
			element.elasticity[gp][21] = ee * (C05 - mi);
			element.elasticity[gp][22] = zeros();
			element.elasticity[gp][23] = zeros();
			element.elasticity[gp][24] = zeros();
			element.elasticity[gp][25] = zeros();
			element.elasticity[gp][26] = zeros();
			element.elasticity[gp][27] = zeros();
			element.elasticity[gp][28] = ee * (C05 - mi);
			element.elasticity[gp][29] = zeros();
			element.elasticity[gp][30] = zeros();
			element.elasticity[gp][31] = zeros();
			element.elasticity[gp][32] = zeros();
			element.elasticity[gp][33] = zeros();
			element.elasticity[gp][34] = zeros();
			element.elasticity[gp][35] = ee * (C05 - mi);
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_ */
