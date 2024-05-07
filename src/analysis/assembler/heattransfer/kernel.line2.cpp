
#include "kernel.h"
#include "analysis/assembler/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runBoundary<Element::CODE::LINE2>(const step::Step &step, SubKernel::Action action, size_t region, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (action) {
		case SubKernel::Action::PREPROCESS:
			setBoundaryKernel<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 2, 1>(boundary[region][interval], action); break;
		case SubKernel::Action::ASSEMBLE:
		case SubKernel::Action::REASSEMBLE:
		case SubKernel::Action::ITERATION:
		case SubKernel::Action::SOLUTION:
			runBoundaryKernel<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 2, 1>(boundary[region][interval], action); break;
		default: break;
		} break;
	case 3:
		switch (action) {
		case SubKernel::Action::PREPROCESS:
			setBoundaryKernel<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 3, 1>(boundary[region][interval], action); break;
		case SubKernel::Action::ASSEMBLE:
		case SubKernel::Action::REASSEMBLE:
		case SubKernel::Action::ITERATION:
		case SubKernel::Action::SOLUTION:
			runBoundaryKernel<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 3, 1>(boundary[region][interval], action); break;
		default: break;
		} break;
	}
}

}
