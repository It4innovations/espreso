
#include "kernel.h"

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runVolumeConductivity(HeatTransfer::SubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.conductivity.conductivity->model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		if (subkernels.advection.isactive) {
			compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::ISOTROPIC, ThermalConductivityConfiguration::MODEL::DIAGONAL>(subkernels, action);
		} else {
			compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::ISOTROPIC, ThermalConductivityConfiguration::MODEL::ISOTROPIC>(subkernels, action);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		if (subkernels.coosystem.rotated) {
			compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::SYMMETRIC>(subkernels, action);
		} else {
			compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::DIAGONAL>(subkernels, action);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::SYMMETRIC, ThermalConductivityConfiguration::MODEL::SYMMETRIC>(subkernels, action);
		break;
	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
		compute<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, ThermalConductivityConfiguration::MODEL::ANISOTROPIC>(subkernels, action);
		break;
	}
}

void HeatTransfer::runVolume(Action action, size_t interval)
{
	switch (subkernels[interval].code) {
	case static_cast<size_t>(Element::CODE::TETRA4   ): runVolumeConductivity<Element::CODE::TETRA4   ,  4, HeatTransferGPC::TETRA4    , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::TETRA10  ): runVolumeConductivity<Element::CODE::TETRA10  , 10, HeatTransferGPC::TETRA10   , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5 ): runVolumeConductivity<Element::CODE::PYRAMID5 ,  5, HeatTransferGPC::PYRAMID5  , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): runVolumeConductivity<Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13 , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PRISMA6  ): runVolumeConductivity<Element::CODE::PRISMA6  ,  6, HeatTransferGPC::PRISMA6   , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PRISMA15 ): runVolumeConductivity<Element::CODE::PRISMA15 , 15, HeatTransferGPC::PRISMA15  , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::HEXA8    ): runVolumeConductivity<Element::CODE::HEXA8    ,  8, HeatTransferGPC::HEXA8     , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::HEXA20   ): runVolumeConductivity<Element::CODE::HEXA20   , 20, HeatTransferGPC::HEXA20    , 3, 3>(subkernels[interval], action); break;
	}
}

}
