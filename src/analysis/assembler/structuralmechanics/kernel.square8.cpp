
#include "kernel.h"
#include "analysis/assembler/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runElement<Element::CODE::SQUARE8>(const step::Step &step, SubKernel::Action action, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setElementKernel<Element::CODE::SQUARE8, 8, StructuralMechanicsGPC::SQUARE8, 2, 2>(subkernels[interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runElementKernel<Element::CODE::SQUARE8, 8, StructuralMechanicsGPC::SQUARE8, 2, 2>(step, subkernels[interval], action); break;
	default: break;
	}
}

template <>
void StructuralMechanics::runBoundary<Element::CODE::SQUARE8>(const step::Step &step, SubKernel::Action action, size_t region, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setBoundaryKernel<Element::CODE::SQUARE8, 8, StructuralMechanicsGPC::SQUARE8, 3, 2>(boundary[region][interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runBoundaryKernel<Element::CODE::SQUARE8, 8, StructuralMechanicsGPC::SQUARE8, 3, 2>(boundary[region][interval], action); break;
	default: break;
	}
}

}
