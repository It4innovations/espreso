
#include "kernel.h"

namespace espreso {

void HeatTransfer::runVolume(Action action, size_t interval)
{
	switch (subkernels[interval].code) {
	case static_cast<size_t>(Element::CODE::TETRA4   ): runConductivity<Element::CODE::TETRA4   ,  4, HeatTransferGPC::TETRA4    , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::TETRA10  ): runConductivity<Element::CODE::TETRA10  , 10, HeatTransferGPC::TETRA10   , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5 ): runConductivity<Element::CODE::PYRAMID5 ,  5, HeatTransferGPC::PYRAMID5  , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): runConductivity<Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13 , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PRISMA6  ): runConductivity<Element::CODE::PRISMA6  ,  6, HeatTransferGPC::PRISMA6   , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::PRISMA15 ): runConductivity<Element::CODE::PRISMA15 , 15, HeatTransferGPC::PRISMA15  , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::HEXA8    ): runConductivity<Element::CODE::HEXA8    ,  8, HeatTransferGPC::HEXA8     , 3, 3>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::HEXA20   ): runConductivity<Element::CODE::HEXA20   , 20, HeatTransferGPC::HEXA20    , 3, 3>(subkernels[interval], action); break;
	}
}

}
