
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runElement<Element::CODE::TRIANGLE6>(const step::Step &step, HeatTransferElementOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setElementKernel<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runElementKernel<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2>(step, operators, action); break;
    default: break;
    }
}

template <>
void runBoundary<Element::CODE::TRIANGLE6>(const step::Step &step, HeatTransferBoundaryOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setBoundaryKernel<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 3, 2>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runBoundaryKernel<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 3, 2>(operators, action); break;
    default: break;
    }
}

}
