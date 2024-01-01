
#include "kernel.h"
#include "analysis/assembler/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runElement<Element::CODE::SQUARE4>(SubKernel::Action action, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setElementKernel<Element::CODE::SQUARE4, 4, StructuralMechanicsGPC::SQUARE4, 2, 2>(subkernels[interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runElementKernel<Element::CODE::SQUARE4, 4, StructuralMechanicsGPC::SQUARE4, 2, 2>(subkernels[interval], action); break;
	default: break;
	}
}

template <>
void StructuralMechanics::runBoundary<Element::CODE::SQUARE4>(SubKernel::Action action, size_t region, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setBoundaryKernel<Element::CODE::SQUARE4, 4, StructuralMechanicsGPC::SQUARE4, 3, 2>(boundary[region][interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runBoundaryKernel<Element::CODE::SQUARE4, 4, StructuralMechanicsGPC::SQUARE4, 3, 2>(boundary[region][interval], action); break;
	default: break;
	}
}

}
