
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ELASTICITY_H_

#include "subkernel.h"

namespace espreso {

struct Elasticity: SubKernel {
	const char* name() const { return "ConductivityKernel"; }

	Elasticity()
	: behaviour(StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN), configuration(nullptr), indirect(false)
	{
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::ITERATION;
		isactive = true;
	}

	void activate(StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour, const LinearElasticPropertiesConfiguration *configuration, bool indirect)
	{
		this->behaviour = behaviour;
		this->configuration = configuration;
		this->indirect = indirect;
		this->isconst = !(this->configuration->needCoordinates() || this->configuration->needTemperature());
	}

	StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour;
	const LinearElasticPropertiesConfiguration *configuration;
	bool indirect;
};

template <size_t gps, size_t ndim, enum ElasticityModel model, class Physics> struct ElasticityKernel;

template <size_t gps, class Physics> struct ElasticityKernel<gps, 2, ElasticityModel::ISOTROPIC, Physics>: Elasticity, Physics {
	ElasticityKernel(const Elasticity &base): Elasticity(base) {}

	void simd(typename Physics::Element &element)
	{
		switch (behaviour) {
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD ex = element.ecf.youngModulus[gp];
				SIMD mi = element.ecf.poissonRatio[gp];
				SIMD C1 = load1(1.);
				SIMD C2 = load1(2.);
				SIMD k = ex * (C1 - mi) / ((C1 + mi) * (C1 - C2 * mi));
				if (indirect) {
					element.ecf.elasticity[gp][0] = k;
					element.ecf.elasticity[gp][1] = k * (mi / (C1 - mi));
					element.ecf.elasticity[gp][2] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
				} else {
					element.elasticity[gp][0] = k;
					element.elasticity[gp][1] = k * (mi / (C1 - mi));
					element.elasticity[gp][2] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
				}
			} break;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD ex = element.ecf.youngModulus[gp];
				SIMD mi = element.ecf.poissonRatio[gp];
				SIMD C1 = load1(1.);
				SIMD C2 = load1(.5);
				SIMD k = ex / (C1 - mi * mi);
				if (indirect) {
					element.ecf.elasticity[gp][0] = k;
					element.ecf.elasticity[gp][1] = k * mi;
					element.ecf.elasticity[gp][2] = k * ((C1 -  mi) * C2);
				} else {
					element.elasticity[gp][0] = k;
					element.elasticity[gp][1] = k * mi;
					element.elasticity[gp][2] = k * ((C1 -  mi) * C2);
				}
			} break;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD ex = element.ecf.youngModulus[gp];
				SIMD mi = element.ecf.poissonRatio[gp];
				SIMD C05 = load1(.5);
				SIMD C1 = load1(1.);
				SIMD C2 = load1(2.);
				SIMD k = ex / ((C1 + mi) * (C1 - C2 * mi));
				if (indirect) {
					element.ecf.elasticity[gp][0] = k * (C1 - mi);
					element.ecf.elasticity[gp][1] = k * mi;
					element.ecf.elasticity[gp][2] = k * (C1 - C2 * mi) * C05;
				} else {
					element.elasticity[gp][0] = k * (C1 - mi);
					element.elasticity[gp][1] = k * mi;
					element.elasticity[gp][2] = k * (C1 - C2 * mi) * C05;
				}
			} break;
		}
	}
};

template <size_t gps, class Physics> struct ElasticityKernel<gps, 3, ElasticityModel::ISOTROPIC, Physics>: Elasticity, Physics {
	ElasticityKernel(const Elasticity &base): Elasticity(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp];
			SIMD mi = element.ecf.poissonRatio[gp];
			SIMD C05 = load1(.5);
			SIMD C1 = load1(1.);
			SIMD C2 = load1(2.);
			SIMD ee = ex / ((C1 + mi) * (C1 - C2 * mi));
			if (indirect) {
				element.ecf.elasticity[gp][0] = ee * (C1 - mi);
				element.ecf.elasticity[gp][1] = ee * mi;
				element.ecf.elasticity[gp][2] = ee * (C05 - mi);
			} else {
				element.elasticity[gp][0] = ee * (C1 - mi);
				element.elasticity[gp][1] = ee * mi;
				element.elasticity[gp][2] = ee * (C05 - mi);
			}
		}

	}
};

template <size_t gps, class Physics> struct ElasticityKernel<gps, 2, ElasticityModel::ORTHOTROPIC, Physics>: Elasticity, Physics {
	ElasticityKernel(const Elasticity &base): Elasticity(base) {}

	void simd(typename Physics::Element &element)
	{
		switch (behaviour) {
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			for (size_t gp = 0; gp < gps; ++gp) {

			} break;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
			for (size_t gp = 0; gp < gps; ++gp) {

			} break;
		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
			for (size_t gp = 0; gp < gps; ++gp) {

			} break;
		}
	}
};

template <size_t gps, class Physics> struct ElasticityKernel<gps, 3, ElasticityModel::ORTHOTROPIC, Physics>: Elasticity, Physics {
	ElasticityKernel(const Elasticity &base): Elasticity(base) {}

	// C
	// 0 1 2 _ _ _
	//   3 4 _ _ _
	//     5 _ _ _
	//       6 _ _
	//         7 _
	//           8
	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD ex = element.ecf.youngModulus[gp][0];
			SIMD ey = element.ecf.youngModulus[gp][1];
			SIMD ez = element.ecf.youngModulus[gp][2];
			SIMD miXY = element.ecf.poissonRatio[gp][0];
			SIMD miXZ = element.ecf.poissonRatio[gp][1];
			SIMD miYZ = element.ecf.poissonRatio[gp][2];

			SIMD miYX = miXY * ey / ex;
			SIMD miZY = miYZ * ez / ey;
			SIMD miZX = miXZ * ex / ez;

			SIMD C1 = load1(1);
			SIMD ksi = C1 / (C1 - (miXY * miYX + miYZ * miZY + miXZ * miZX) - (miXY * miYZ * miZX + miYX * miZY * miXZ));

			if (indirect) {
				element.ecf.elasticity[gp][0] = ksi * ex * (C1   - miYZ * miZY);
				element.ecf.elasticity[gp][1] = ksi * ey * (miXY + miXZ * miZY);
				element.ecf.elasticity[gp][2] = ksi * ez * (miXZ + miYZ * miXY);
				element.ecf.elasticity[gp][3] = ksi * ey * (C1   - miXZ * miZX);
				element.ecf.elasticity[gp][4] = ksi * ez * (miYZ + miYX * miXZ);
				element.ecf.elasticity[gp][5] = ksi * ez * (C1   - miYX * miXY);
				element.ecf.elasticity[gp][6] = element.ecf.shearModulus[gp][0];
				element.ecf.elasticity[gp][7] = element.ecf.shearModulus[gp][1];
				element.ecf.elasticity[gp][8] = element.ecf.shearModulus[gp][2];
			} else {
				element.elasticity[gp][0] = ksi * ex * (C1   - miYZ * miZY);
				element.elasticity[gp][1] = ksi * ey * (miXY + miXZ * miZY);
				element.elasticity[gp][2] = ksi * ez * (miXZ + miYZ * miXY);
				element.elasticity[gp][3] = ksi * ey * (C1   - miXZ * miZX);
				element.elasticity[gp][4] = ksi * ez * (miYZ + miYX * miXZ);
				element.elasticity[gp][5] = ksi * ez * (C1   - miYX * miXY);
				element.elasticity[gp][6] = element.ecf.shearModulus[gp][0];
				element.elasticity[gp][7] = element.ecf.shearModulus[gp][1];
				element.elasticity[gp][8] = element.ecf.shearModulus[gp][2];
			}
		}
	}
};

template <size_t gps, class Physics> struct ElasticityKernel<gps, 2, ElasticityModel::ANISOTROPIC, Physics>: Elasticity, Physics {
	ElasticityKernel(const Elasticity &base): Elasticity(base) {}

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t gps, class Physics> struct ElasticityKernel<gps, 3, ElasticityModel::ANISOTROPIC, Physics>: Elasticity, Physics {
	ElasticityKernel(const Elasticity &base): Elasticity(base) {}

	void simd(typename Physics::Element &element)
	{

	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_ELASTICITY_H_ */
