
#include "kernel.h"
#include "analysis/assembler/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runElement<Element::CODE::HEXA8>(const step::Step &step, SubKernel::Action action, size_t interval)
{
	switch (action) {
	case SubKernel::Action::PREPROCESS:
		setElementKernel<Element::CODE::HEXA8, 8, StructuralMechanicsGPC::HEXA8, 3, 3>(subkernels[interval], action); break;
	case SubKernel::Action::ASSEMBLE:
	case SubKernel::Action::REASSEMBLE:
	case SubKernel::Action::ITERATION:
	case SubKernel::Action::SOLUTION:
		runElementKernel<Element::CODE::HEXA8, 8, StructuralMechanicsGPC::HEXA8, 3, 3>(step, subkernels[interval], action); break;
	default: break;
	}
}

}
