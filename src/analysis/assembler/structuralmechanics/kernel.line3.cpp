
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runBoundary<Element::CODE::LINE3>(const step::Step &step, StructuralMechanicsFaceOperators &operators, SubKernel::Action action)
{
    switch (info::mesh->dimension) {
    case 2:
        switch (action) {
        case SubKernel::Action::PREPROCESS:
            setBoundaryKernel<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 2, 1>(operators, action); break;
        case SubKernel::Action::ASSEMBLE:
        case SubKernel::Action::REASSEMBLE:
        case SubKernel::Action::ITERATION:
        case SubKernel::Action::SOLUTION:
            runBoundaryKernel<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 2, 1>(operators, action); break;
        default: break;
        } break;
    case 3:
        switch (action) {
        case SubKernel::Action::PREPROCESS:
            setBoundaryKernel<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 3, 1>(operators, action); break;
        case SubKernel::Action::ASSEMBLE:
        case SubKernel::Action::REASSEMBLE:
        case SubKernel::Action::ITERATION:
        case SubKernel::Action::SOLUTION:
            runBoundaryKernel<Element::CODE::LINE3, 3, StructuralMechanicsGPC::LINE3, 3, 1>(operators, action); break;
        default: break;
        } break;
    }
}

}
