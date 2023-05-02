
#include "kernel.h"

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour>
void runPlaneElasticity(StructuralMechanics::SubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.elasticity.configuration->model) {
	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
		if (subkernels.coosystem.rotated || subkernels.plasticity.isactive) {
			compute<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ISOTROPIC, ElasticityModel::SYMMETRIC>(subkernels, action);
		} else {
			compute<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ISOTROPIC, ElasticityModel::ISOTROPIC>(subkernels, action);
		}
		break;
	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
		if (subkernels.coosystem.rotated) {
			compute<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ORTHOTROPIC, ElasticityModel::SYMMETRIC>(subkernels, action);
		} else {
			compute<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ORTHOTROPIC, ElasticityModel::ORTHOTROPIC>(subkernels, action);
		}
		break;
	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
		compute<code, nodes, gps, ndim, edim, behaviour, ElasticityModel::ANISOTROPIC, ElasticityModel::ANISOTROPIC>(subkernels, action);
		break;
	}
}

void StructuralMechanics::runPlane(Action action, size_t interval)
{
	switch (subkernels[interval].code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): runPlaneElasticity<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, Behaviour::PLANE>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): runPlaneElasticity<Element::CODE::TRIANGLE6, 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, Behaviour::PLANE>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): runPlaneElasticity<Element::CODE::SQUARE4  , 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, Behaviour::PLANE>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::SQUARE8  ): runPlaneElasticity<Element::CODE::SQUARE8  , 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, Behaviour::PLANE>(subkernels[interval], action); break;
	}
}

}
