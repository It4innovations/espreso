
#include "kernel.h"
#include "analysis/assembler/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runElement<Element::CODE::TETRA4>(const step::Step &step, SubKernel::Action action, size_t interval)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setElementKernel<Element::CODE::TETRA4, 4, StructuralMechanicsGPC::TETRA4, 3, 3>(subkernels[interval], action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runElementKernel<Element::CODE::TETRA4, 4, StructuralMechanicsGPC::TETRA4, 3, 3>(step, subkernels[interval], action); break;
    default: break;
    }
}

}
