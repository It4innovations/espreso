
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SIGMA_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SIGMA_H_

#include "element.h"
#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct Sigma: SubKernel {
	const char* name() const { return "Sigma"; }

	Sigma()
	: behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN),
	  model(LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC),
	  rotated(false)
	{
		isconst = false;
		action = SubKernel::SOLUTION;
	}

	void activate(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour, LinearElasticPropertiesConfiguration::MODEL model, bool rotated)
	{
		this->behaviour = behaviour;
		this->model = model;
		this->rotated = rotated;
		this->isactive = 1;
	}

	StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
	LinearElasticPropertiesConfiguration::MODEL model;
	bool rotated;
};

template <size_t nodes, size_t ndim> struct SigmaKernel;

template <size_t nodes>
struct SigmaKernel<nodes, 2>: Sigma {
	SigmaKernel(const Sigma &base): Sigma(base) {}

	template <typename Element>
	void reset(Element &element)
	{

	}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{

	}
};

template <size_t nodes>
struct SigmaKernel<nodes, 3>: Sigma {
	SigmaKernel(const Sigma &base): Sigma(base) {}

	template <typename Element>
	void reset(Element &element)
	{
		element.sigma[0] = element.sigma[1] = element.sigma[2] = element.sigma[3] = element.sigma[4] = element.sigma[5] = zeros();
	}

	template <typename Element>
	void simd(Element &element, size_t gp)
	{
		if (rotated && model != LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC) {
			SIMD c00 = element.elasticity[0];
			SIMD c01 = element.elasticity[1], c11 = element.elasticity[ 6];
			SIMD c02 = element.elasticity[2], c12 = element.elasticity[ 7], c22 = element.elasticity[11];
			SIMD c03 = element.elasticity[3], c13 = element.elasticity[ 8], c23 = element.elasticity[12], c33 = element.elasticity[15];
			SIMD c04 = element.elasticity[4], c14 = element.elasticity[ 9], c24 = element.elasticity[13], c34 = element.elasticity[16], c44 = element.elasticity[18];
			SIMD c05 = element.elasticity[5], c15 = element.elasticity[10], c25 = element.elasticity[14], c35 = element.elasticity[17], c45 = element.elasticity[19], c55 = element.elasticity[20];

			SIMD uB0 = element.smallStrainTensor[0];
			SIMD uB1 = element.smallStrainTensor[1];
			SIMD uB2 = element.smallStrainTensor[2];
			SIMD uB3 = element.smallStrainTensor[3];
			SIMD uB4 = element.smallStrainTensor[4];
			SIMD uB5 = element.smallStrainTensor[5];

			element.sigma[0] = element.sigma[0] + uB0 * c00 + uB1 * c01 + uB2 * c02 + uB3 * c03 + uB4 * c04 + uB5 * c05;
			element.sigma[1] = element.sigma[1] + uB0 * c01 + uB1 * c11 + uB2 * c12 + uB3 * c13 + uB4 * c14 + uB5 * c15;
			element.sigma[2] = element.sigma[2] + uB0 * c02 + uB1 * c12 + uB2 * c22 + uB3 * c23 + uB4 * c24 + uB5 * c25;
			element.sigma[3] = element.sigma[3] + uB0 * c03 + uB1 * c13 + uB2 * c23 + uB3 * c33 + uB4 * c34 + uB5 * c35;
			element.sigma[4] = element.sigma[4] + uB0 * c04 + uB1 * c14 + uB2 * c24 + uB3 * c34 + uB4 * c44 + uB5 * c45;
			element.sigma[5] = element.sigma[5] + uB0 * c05 + uB1 * c15 + uB2 * c25 + uB3 * c35 + uB4 * c45 + uB5 * c55;
			return;
		}
		switch (model) {
		case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
		{
			// C
			// 0 1 1 _ _ _
			//   0 1 _ _ _
			//     0 _ _ _
			//       2 _ _
			//         2 _
			//           2
			SIMD c0 = element.elasticity[0];
			SIMD c1 = element.elasticity[1];
			SIMD c2 = element.elasticity[21];

			element.sigma[0] = element.sigma[0] + element.smallStrainTensor[0] * c0 + element.smallStrainTensor[1] * c1 + element.smallStrainTensor[2] * c1;
			element.sigma[1] = element.sigma[1] + element.smallStrainTensor[0] * c1 + element.smallStrainTensor[1] * c0 + element.smallStrainTensor[2] * c1;
			element.sigma[2] = element.sigma[2] + element.smallStrainTensor[0] * c1 + element.smallStrainTensor[1] * c1 + element.smallStrainTensor[2] * c0;
			element.sigma[3] = element.sigma[3] + element.smallStrainTensor[3] * c2;
			element.sigma[4] = element.sigma[4] + element.smallStrainTensor[4] * c2;
			element.sigma[5] = element.sigma[5] + element.smallStrainTensor[5] * c2;
		} break;
		case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
		{
			// C
			// 0 1 2 _ _ _
			//   3 4 _ _ _
			//     5 _ _ _
			//       6 _ _
			//         7 _
			//           8
			SIMD c00 = element.elasticity[0];
			SIMD c01 = element.elasticity[1], c11 = element.elasticity[3];
			SIMD c02 = element.elasticity[2], c12 = element.elasticity[4], c22 = element.elasticity[5];
			SIMD c33 = element.elasticity[6];
			SIMD c44 = element.elasticity[7];
			SIMD c55 = element.elasticity[8];

			element.sigma[0] = element.sigma[0] + element.smallStrainTensor[0] * c00 + element.smallStrainTensor[1] * c01 + element.smallStrainTensor[2] * c02;
			element.sigma[1] = element.sigma[1] + element.smallStrainTensor[0] * c01 + element.smallStrainTensor[1] * c11 + element.smallStrainTensor[2] * c12;
			element.sigma[2] = element.sigma[2] + element.smallStrainTensor[0] * c02 + element.smallStrainTensor[1] * c12 + element.smallStrainTensor[2] * c22;
			element.sigma[3] = element.sigma[3] + element.smallStrainTensor[3] * c33;
			element.sigma[4] = element.sigma[4] + element.smallStrainTensor[4] * c44;
			element.sigma[5] = element.sigma[5] + element.smallStrainTensor[5] * c55;
		} break;
		case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
		{
			SIMD c00 = element.elasticity[0], c10 = element.elasticity[ 6], c20 = element.elasticity[12], c30 = element.elasticity[18], c40 = element.elasticity[24], c50 = element.elasticity[30];
			SIMD c01 = element.elasticity[1], c11 = element.elasticity[ 7], c21 = element.elasticity[13], c31 = element.elasticity[19], c41 = element.elasticity[25], c51 = element.elasticity[31];
			SIMD c02 = element.elasticity[2], c12 = element.elasticity[ 8], c22 = element.elasticity[14], c32 = element.elasticity[20], c42 = element.elasticity[26], c52 = element.elasticity[32];
			SIMD c03 = element.elasticity[3], c13 = element.elasticity[ 9], c23 = element.elasticity[15], c33 = element.elasticity[21], c43 = element.elasticity[27], c53 = element.elasticity[33];
			SIMD c04 = element.elasticity[4], c14 = element.elasticity[10], c24 = element.elasticity[16], c34 = element.elasticity[22], c44 = element.elasticity[28], c54 = element.elasticity[34];
			SIMD c05 = element.elasticity[5], c15 = element.elasticity[11], c25 = element.elasticity[17], c35 = element.elasticity[23], c45 = element.elasticity[29], c55 = element.elasticity[35];

			SIMD uB0 = element.smallStrainTensor[0];
			SIMD uB1 = element.smallStrainTensor[1];
			SIMD uB2 = element.smallStrainTensor[2];
			SIMD uB3 = element.smallStrainTensor[3];
			SIMD uB4 = element.smallStrainTensor[4];
			SIMD uB5 = element.smallStrainTensor[5];

			element.sigma[0] = element.sigma[0] + uB0 * c00 + uB1 * c10 + uB2 * c20 + uB3 * c30 + uB4 * c40 + uB5 * c50;
			element.sigma[1] = element.sigma[1] + uB0 * c01 + uB1 * c11 + uB2 * c21 + uB3 * c31 + uB4 * c41 + uB5 * c51;
			element.sigma[2] = element.sigma[2] + uB0 * c02 + uB1 * c12 + uB2 * c22 + uB3 * c32 + uB4 * c42 + uB5 * c52;
			element.sigma[3] = element.sigma[3] + uB0 * c03 + uB1 * c13 + uB2 * c23 + uB3 * c33 + uB4 * c43 + uB5 * c53;
			element.sigma[4] = element.sigma[4] + uB0 * c04 + uB1 * c14 + uB2 * c24 + uB3 * c34 + uB4 * c44 + uB5 * c54;
			element.sigma[5] = element.sigma[5] + uB0 * c05 + uB1 * c15 + uB2 * c25 + uB3 * c35 + uB4 * c45 + uB5 * c55;
		} break;
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_SIGMA_H_ */
