
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runElement<Element::CODE::TRIANGLE3>(const step::Step &step, StructuralMechanicsElementOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setElementKernel<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runElementKernel<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2>(step, operators, action); break;
    default: break;
    }
}

template <>
void runBoundary<Element::CODE::TRIANGLE3>(const step::Step &step, StructuralMechanicsBoundaryOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setBoundaryKernel<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runBoundaryKernel<Element::CODE::TRIANGLE3, 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2>(operators, action); break;
    default: break;
    }
}

}
