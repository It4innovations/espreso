
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runElement<Element::CODE::HEXA8>(const step::Step &step, StructuralMechanicsElementOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setElementKernel<Element::CODE::HEXA8, 8, StructuralMechanicsGPC::HEXA8, 3, 3>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runElementKernel<Element::CODE::HEXA8, 8, StructuralMechanicsGPC::HEXA8, 3, 3>(step, operators, action); break;
    default: break;
    }
}

}
