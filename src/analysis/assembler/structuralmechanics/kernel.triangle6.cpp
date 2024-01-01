
#include "kernel.h"
#include "analysis/assembler/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runElement<Element::CODE::TRIANGLE3>(SubKernel::Action action, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setElementKernel<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2>(subkernels[interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runElementKernel<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2>(subkernels[interval], action); break;
	default: break;
	}
}

template <>
void StructuralMechanics::runBoundary<Element::CODE::TRIANGLE3>(SubKernel::Action action, size_t region, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setBoundaryKernel<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2>(boundary[region][interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runBoundaryKernel<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2>(boundary[region][interval], action); break;
	default: break;
	}
}

}
