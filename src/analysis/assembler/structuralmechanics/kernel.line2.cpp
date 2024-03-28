
#include "kernel.h"
#include "analysis/assembler/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runBoundary<Element::CODE::LINE2>(SubKernel::Action action, size_t region, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (action) {
		case SubKernel::Action::PREPROCESS:
			setBoundaryKernel<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 2, 1>(boundary[region][interval], action); break;
		case SubKernel::Action::ASSEMBLE:
		case SubKernel::Action::REASSEMBLE:
		case SubKernel::Action::ITERATION:
		case SubKernel::Action::SOLUTION:
			runBoundaryKernel<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 2, 1>(boundary[region][interval], action); break;
		default: break;
		} break;
	case 3:
		switch (action) {
		case SubKernel::Action::PREPROCESS:
			setBoundaryKernel<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 3, 1>(boundary[region][interval], action); break;
		case SubKernel::Action::ASSEMBLE:
		case SubKernel::Action::REASSEMBLE:
		case SubKernel::Action::ITERATION:
		case SubKernel::Action::SOLUTION:
			runBoundaryKernel<Element::CODE::LINE2, 2, StructuralMechanicsGPC::LINE2, 3, 1>(boundary[region][interval], action); break;
		default: break;
		} break;
	}
}

}