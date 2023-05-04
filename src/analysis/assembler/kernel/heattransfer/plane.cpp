
#include "kernel.h"

namespace espreso {

void HeatTransfer::runPlane(Action action, size_t interval)
{
	switch (subkernels[interval].code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): runConductivity<Element::CODE::TRIANGLE3, 3, HeatTransferGPC::TRIANGLE3, 2, 2>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): runConductivity<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): runConductivity<Element::CODE::SQUARE4  , 4, HeatTransferGPC::SQUARE4  , 2, 2>(subkernels[interval], action); break;
	case static_cast<size_t>(Element::CODE::SQUARE8  ): runConductivity<Element::CODE::SQUARE8  , 8, HeatTransferGPC::SQUARE8  , 2, 2>(subkernels[interval], action); break;
	}
}

}
