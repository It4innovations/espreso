
#include "kernel.h"
#include "analysis/assembler/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runElement<Element::CODE::SQUARE8>(SubKernel::Action action, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setElementKernel<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 2, 2>(subkernels[interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runElementKernel<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 2, 2>(subkernels[interval], action); break;
	default: break;
	}
}

template <>
void HeatTransfer::runBoundary<Element::CODE::SQUARE8>(SubKernel::Action action, size_t region, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setBoundaryKernel<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 3, 2>(boundary[region][interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runBoundaryKernel<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 3, 2>(boundary[region][interval], action); break;
	default: break;
	}
}

}
