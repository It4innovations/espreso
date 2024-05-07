
#include "kernel.h"
#include "analysis/assembler/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runBoundary<Element::CODE::POINT1>(const step::Step &step, SubKernel::Action action, size_t region, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (action) {
		case SubKernel::Action::PREPROCESS:
			setNodeKernel<2>(boundary[region][interval], action); break;
		case SubKernel::Action::ASSEMBLE:
		case SubKernel::Action::REASSEMBLE:
		case SubKernel::Action::ITERATION:
		case SubKernel::Action::SOLUTION:
			runNodeKernel<2>(boundary[region][interval], action); break;
		default: break;
		} break;
	case 3:
		switch (action) {
		case SubKernel::Action::PREPROCESS:
			setNodeKernel<3>(boundary[region][interval], action); break;
		case SubKernel::Action::ASSEMBLE:
		case SubKernel::Action::REASSEMBLE:
		case SubKernel::Action::ITERATION:
		case SubKernel::Action::SOLUTION:
			runNodeKernel<3>(boundary[region][interval], action); break;
		default: break;
		} break;
	}
}

}
