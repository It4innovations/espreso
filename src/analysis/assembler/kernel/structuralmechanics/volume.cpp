
#include "kernel.h"

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour>
void runVolumeElasticity(StructuralMechanics::SubKernels &subkernels, Assembler::Action action)
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

void StructuralMechanics::runVolume(Action action, size_t interval)
{
	switch (subkernels[interval].code) {
	case static_cast<size_t>(Element::CODE::TETRA4   ): runVolumeElasticity<Element::CODE::TETRA4   ,  4, StructuralMechanicsGPC::TETRA4    , 3, 3, Behaviour::VOLUME>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::TETRA10  ): runVolumeElasticity<Element::CODE::TETRA10  , 10, StructuralMechanicsGPC::TETRA10   , 3, 3, Behaviour::VOLUME>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5 ): runVolumeElasticity<Element::CODE::PYRAMID5 ,  5, StructuralMechanicsGPC::PYRAMID5  , 3, 3, Behaviour::VOLUME>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): runVolumeElasticity<Element::CODE::PYRAMID13, 13, StructuralMechanicsGPC::PYRAMID13 , 3, 3, Behaviour::VOLUME>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PRISMA6  ): runVolumeElasticity<Element::CODE::PRISMA6  ,  6, StructuralMechanicsGPC::PRISMA6   , 3, 3, Behaviour::VOLUME>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PRISMA15 ): runVolumeElasticity<Element::CODE::PRISMA15 , 15, StructuralMechanicsGPC::PRISMA15  , 3, 3, Behaviour::VOLUME>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::HEXA8    ): runVolumeElasticity<Element::CODE::HEXA8    ,  8, StructuralMechanicsGPC::HEXA8     , 3, 3, Behaviour::VOLUME>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::HEXA20   ): runVolumeElasticity<Element::CODE::HEXA20   , 20, StructuralMechanicsGPC::HEXA20    , 3, 3, Behaviour::VOLUME>(subkernels[interval], action); break;
	}
}

}

