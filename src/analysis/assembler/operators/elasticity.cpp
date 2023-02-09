
#include "elasticity.h"

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

		bool isconst = true, rotate = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;

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
			break;
		case 3:
			switch (mat->linear_elastic_properties.model) {
			case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.young_modulus.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.youngModulus[gp][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				elementOps[interval].push_back(generateExpression3D<ExternalGPsExpression>(interval, etype[interval], le.poisson_ratio.get(0, 0).evaluator,
							[] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.poissonRatio[gp][s] = value; }));
				isconst &= elementOps[interval].back()->isconst;
				break;
			case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
			case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
				break;
			}
			elementOps[interval].push_back(generateElementOperator3D<ElasticityIsotropicVolume>(interval, etype[interval]));
			elementOps[interval].back()->isconst &= isconst;
			break;
		}
	}
}

}
