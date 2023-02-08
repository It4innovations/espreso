
#include "elasticity.h"
#include "gausspoints.h"
#include "copy.h"

#include "analysis/assembler/module/structuralmechanics.h"

#include "config/ecf/material/material.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"

namespace espreso {

void elasticity(StructuralMechanics &module)
{
//	bool cooToGP = false;
//	for(size_t interval = 0; !cooToGP && interval < info::mesh->elements->eintervals.size(); ++interval) {
//		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//		if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//			cooToGP = true;
//		}
//	}
//
//	if (module.coords.gp.data == nullptr && cooToGP) {
//		module.controller.addInput(module.coords.gp, module.coords.node);
//		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
//			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//			if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//				module.coords.gp.isconst[interval] = true;
//				module.coords.gp.update[interval] = -1;
//			}
//		}
//	}
//	if (module.coords.gp.data == nullptr) {
//		module.controller.prepare(module.coords.gp);
//		for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
//			const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//			if (mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN) {
//				if (info::mesh->dimension == 2) {
//					module.elementOps[interval].emplace_back(instantiate<StructuralMechanics::NGP, 2, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
//				}
//				if (info::mesh->dimension == 3) {
//					module.elementOps[interval].emplace_back(instantiate<StructuralMechanics::NGP, 3, FromNodesToGaussPoints>(interval, module.controller, module.integration.N, module.coords.node, module.coords.gp));
//				}
//			}
//		}
//	}
//
//	if (cooToGP) {
//		if (info::mesh->dimension == 2) {
//			module.controller.addInput(module.cooSystem.angle2D, module.coords.gp);
//		}
//		if (info::mesh->dimension == 3) {
//			module.controller.addInput(module.cooSystem.angle3D, module.coords.gp);
//		}
//	}
//
//	if (info::mesh->dimension == 2) {
////		module.controller.addInput(module.cooSystem.angle2D, module.cooSystem.cartesian2D, module.cooSystem.cylindric);
//		module.controller.prepare(module.cooSystem.angle2D);
//		module.controller.addInput(module.material.elasticityPlane, module.cooSystem.angle2D);
//		module.controller.addInput(module.material.elasticityAxisymm, module.cooSystem.angle2D);
//	}
//	if (info::mesh->dimension == 3) {
////		module.controller.addInput(module.cooSystem.angle3D, module.cooSystem.cartesian3D, module.cooSystem.cylindric, module.cooSystem.spherical);
//		module.controller.prepare(module.cooSystem.angle3D);
//		module.controller.addInput(module.material.elasticity3D, module.cooSystem.angle3D);
//	}
//
//	module.controller.addInput(module.material.elasticityPlane, module.material.model.isoYoungModulus, module.material.model.isoPoissonRatio, module.material.model.youngModulus, module.material.model.poissonRatio);
//	module.controller.addInput(module.material.elasticityAxisymm, module.material.model.isoYoungModulus, module.material.model.isoPoissonRatio, module.material.model.youngModulus, module.material.model.poissonRatio);
//	module.controller.addInput(module.material.elasticity3D, module.material.model.isoYoungModulus, module.material.model.isoPoissonRatio, module.material.model.youngModulus, module.material.model.poissonRatio, module.material.model.shearModulus, module.material.model.anisotropic3D);
//
//	if (info::mesh->dimension == 2) {
//		switch (module.settings.element_behaviour) {
//		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
//		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
//		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
//			module.controller.prepare(module.material.elasticityPlane);
//			break;
//		case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
//			module.controller.prepare(module.material.elasticityAxisymm);
//			break;
//		}
//	}
//	if (info::mesh->dimension == 3) {
//		module.controller.prepare(module.material.elasticity3D);
//	}
//
//	for(size_t interval = 0; interval < info::mesh->elements->eintervals.size(); ++interval) {
//		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
//		if (info::mesh->dimension == 2) {
//			switch (mat->linear_elastic_properties.model) {
//			case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
//				switch (module.settings.element_behaviour) {
//				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
//					module.elementOps[interval].emplace_back(instantiate<StructuralMechanics::NGP, ElasticityIsotropicPlaneStrain>(interval, module.controller,
//							module.material.model.isoYoungModulus, module.material.model.isoPoissonRatio, module.material.model.shearModulus, module.material.elasticityPlane));
//					break;
//				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
//					module.elementOps[interval].emplace_back(instantiate<StructuralMechanics::NGP, ElasticityIsotropicPlaneStress>(interval, module.controller,
//							module.material.model.isoYoungModulus, module.material.model.isoPoissonRatio, module.material.model.shearModulus, module.material.elasticityPlane));
//					break;
//				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
//					module.elementOps[interval].emplace_back(instantiate<StructuralMechanics::NGP, ElasticityIsotropicPlaneStress>(interval, module.controller,
//							module.material.model.isoYoungModulus, module.material.model.isoPoissonRatio, module.material.model.shearModulus, module.material.elasticityPlane));
//					break;
//				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
//					module.elementOps[interval].emplace_back(instantiate<StructuralMechanics::NGP, ElasticityIsotropicAxisymmetric>(interval, module.controller,
//							module.material.model.isoYoungModulus, module.material.model.isoPoissonRatio, module.material.model.shearModulus, module.material.elasticityAxisymm));
//					break;
//				default:
//					break;
//				}
//				break;
//			case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
//				break;
//			case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
//				break;
//			}
//		}
//		if (info::mesh->dimension == 3) {
//			switch (mat->linear_elastic_properties.model) {
//			case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
//				module.elementOps[interval].emplace_back(instantiate<StructuralMechanics::NGP, ElasticityIsotropic3D>(interval, module.controller,
//						module.material.model.isoYoungModulus, module.material.model.isoPoissonRatio, module.material.model.shearModulus, module.material.elasticity3D));
//				break;
//			case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
//				break;
//			case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
//				break;
//			}
//		}
//	}
}

}
