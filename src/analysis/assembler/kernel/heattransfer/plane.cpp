
#include "kernel.h"

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runPlaneConductivity(HeatTransfer::SubKernels &subkernels, Assembler::Action action)
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

void HeatTransfer::runPlane(Action action, size_t interval)
{
	switch (subkernels[interval].code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): runPlaneConductivity<Element::CODE::TRIANGLE3, 3, HeatTransferGPC::TRIANGLE3, 2, 2>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): runPlaneConductivity<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): runPlaneConductivity<Element::CODE::SQUARE4  , 4, HeatTransferGPC::SQUARE4  , 2, 2>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::SQUARE8  ): runPlaneConductivity<Element::CODE::SQUARE8  , 8, HeatTransferGPC::SQUARE8  , 2, 2>(subkernels[interval], action); break;
	}
}

}
