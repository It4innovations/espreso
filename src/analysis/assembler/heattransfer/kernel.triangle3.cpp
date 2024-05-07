
#include "kernel.h"
#include "analysis/assembler/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runElement<Element::CODE::TRIANGLE6>(const step::Step &step, SubKernel::Action action, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setElementKernel<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2>(subkernels[interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runElementKernel<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2>(step, subkernels[interval], action); break;
	default: break;
	}
}

template <>
void HeatTransfer::runBoundary<Element::CODE::TRIANGLE6>(const step::Step &step, SubKernel::Action action, size_t region, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setBoundaryKernel<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 3, 2>(boundary[region][interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runBoundaryKernel<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 3, 2>(boundary[region][interval], action); break;
	default: break;
	}
}

}
