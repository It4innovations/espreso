
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runElement<Element::CODE::SQUARE8>(const step::Step &step, HeatTransferElementOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setElementKernel<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 2, 2>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runElementKernel<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 2, 2>(step, operators, action); break;
    default: break;
    }
}

template <>
void runBoundary<Element::CODE::SQUARE8>(const step::Step &step, HeatTransferBoundaryOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setBoundaryKernel<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 3, 2>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runBoundaryKernel<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 3, 2>(operators, action); break;
    default: break;
    }
}

}
