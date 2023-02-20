
#include "elasticity.h"
#include "elasticity.coordinatesystem.h"

#include "analysis/assembler/module/structuralmechanics.h"
#include "analysis/assembler/module/structuralmechanics.generator.h"

#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"

namespace espreso {

void StructuralMechanics::generateElasticity()
{
	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];

		bool isconst = true, rotate = mat->linear_elastic_properties.model != LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC;
		if (mat->linear_elastic_properties.model != LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC) {
			if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
				if (info::mesh->dimension == 2) {
					rotate &= mat->coordinate_system.rotation.z.isset;
				}
				if (info::mesh->dimension == 3) {
					rotate &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
				}
			}
		}

		elementOps[interval].push_back(generateExpression<ExternalGPsExpression>(interval, etype[interval], mat->density.evaluator,
					[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.density[gp][s] = value; }));
		elementOps[interval].push_back(generateExpression<ExternalGPsExpression>(interval, etype[interval], mat->heat_capacity.evaluator,
					[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.heatCapacity[gp][s] = value; }));

		const LinearElasticPropertiesConfiguration &le = mat->linear_elastic_properties;
		switch (info::mesh->dimension) {
		case 2:
			switch (mat->linear_elastic_properties.model) {
			case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
				elementOps[interval].push_back(generateExpression2D<ExternalGPsExpression>(interval, etype[interval], le.young_modulus.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.youngModulus[gp][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression2D<ExternalGPsExpression>(interval, etype[interval], le.poisson_ratio.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.poissonRatio[gp][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				break;
			case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
			case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
				break;
			}

			switch (settings.element_behaviour) {
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
				elementOps[interval].push_back(generateElementOperator2D<ElasticityIsotropicPlaneStrain>(interval, etype[interval]));
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
				elementOps[interval].push_back(generateElementOperator2D<ElasticityIsotropicPlaneStress>(interval, etype[interval]));
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
				elementOps[interval].push_back(generateElementOperator2D<ElasticityIsotropicPlaneStress>(interval, etype[interval]));
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
				elementOps[interval].push_back(generateElementOperator2D<ElasticityIsotropicPlaneAxisymmetric>(interval, etype[interval]));
				break;
			}
			elementOps[interval].back()->isconst &= isconst;

			if (rotate) {
				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
//					elementOps[interval].push_back(generateExpression2D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.rotation.z.evaluator,
//								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
//					isconst &= elementOps[interval].back()->isconst;
//					elementOps[interval].push_back(generateElementOperator2D<CoordinateSystemCartesian>(interval, etype[interval]));
//					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
//					elementOps[interval].push_back(generateExpression2D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.center.x.evaluator,
//								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
//					isconst &= elementOps[interval].back()->isconst;
//					elementOps[interval].push_back(generateExpression2D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.center.y.evaluator,
//								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][1][s] = value; }));
//					isconst &= elementOps[interval].back()->isconst;
//					elementOps[interval].push_back(generateElementOperator2D<CoordinateSystemCylindric>(interval, etype[interval]));
//					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::SPHERICAL:
					break;
				}
//				elementOps[interval].push_back(generateElementOperator2D<CoordinateSystemApply>(interval, etype[interval]));
//				elementOps[interval].back()->isconst &= isconst;
			} else {
				elementOps[interval].push_back(generateElementOperator2D<ElasticityCoordinateSystemCopy>(interval, etype[interval]));
				elementOps[interval].back()->isconst &= isconst;
			}
			break;
		case 3:
			switch (mat->linear_elastic_properties.model) {
			case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.young_modulus.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.youngModulus[gp][0][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.poisson_ratio.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.poissonRatio[gp][0][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;

				elementOps[interval].push_back(generateElementOperator3D<ElasticityIsotropicVolume>(interval, etype[interval]));
				elementOps[interval].back()->isconst &= isconst;
				break;
			case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.young_modulus.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.youngModulus[gp][0][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.young_modulus.get(1, 1).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.youngModulus[gp][1][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.young_modulus.get(2, 2).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.youngModulus[gp][2][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.poisson_ratio.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.poissonRatio[gp][0][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.poisson_ratio.get(1, 1).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.poissonRatio[gp][1][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.poisson_ratio.get(2, 2).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.poissonRatio[gp][2][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.shear_modulus.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.shearModulus[gp][0][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.shear_modulus.get(1, 1).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.shearModulus[gp][1][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.shear_modulus.get(2, 2).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.shearModulus[gp][2][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;

				elementOps[interval].push_back(generateElementOperator3D<ElasticityOrthotropicVolume>(interval, etype[interval]));
				elementOps[interval].back()->isconst &= isconst;
				break;
			case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
				break;
			}

			if (rotate) {
				switch (mat->coordinate_system.type) {
				case CoordinateSystemConfiguration::TYPE::CARTESIAN:
					elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.rotation.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.rotation.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.rotation.z.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementOperator3D<ElasticityCoordinateSystemCartesian>(interval, etype[interval]));
					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
					elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.center.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.center.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementOperator3D<ElasticityCoordinateSystemCylindric>(interval, etype[interval]));
					elementOps[interval].back()->isconst &= isconst;
					break;
				case CoordinateSystemConfiguration::TYPE::SPHERICAL:
					elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.center.x.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][0][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.center.y.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][1][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], mat->coordinate_system.center.z.evaluator,
								[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.center[gp][2][s] = value; }));
					isconst &= elementOps[interval].back()->isconst;
					elementOps[interval].push_back(generateElementOperator3D<ElasticityCoordinateSystemSpherical>(interval, etype[interval]));
					elementOps[interval].back()->isconst &= isconst;
					break;
				}
				elementOps[interval].push_back(generateElementOperator3D<ElasticityCoordinateSystemApply>(interval, etype[interval]));
				elementOps[interval].back()->isconst &= isconst;
			} else {
				elementOps[interval].push_back(generateElementOperator3D<ElasticityCoordinateSystemCopy>(interval, etype[interval]));
				elementOps[interval].back()->isconst &= isconst;
			}
			break;
		}

	}
}

}
