
#include "kernel.h"
#include "analysis/assembler/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runElement<Element::CODE::PRISMA15>(const step::Step &step, SubKernel::Action action, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setElementKernel<Element::CODE::PRISMA15, 15, HeatTransferGPC::PRISMA15, 3, 3>(subkernels[interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runElementKernel<Element::CODE::PRISMA15, 15, HeatTransferGPC::PRISMA15, 3, 3>(step, subkernels[interval], action); break;
	default: break;
	}
}

}
