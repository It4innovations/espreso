
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runBoundary<Element::CODE::LINE2>(const step::Step &step, HeatTransferBoundaryOperators &operators, SubKernel::Action action)
{
    switch (info::mesh->dimension) {
    case 2:
        switch (action) {
        case SubKernel::Action::PREPROCESS:
            setBoundaryKernel<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 2, 1>(operators, action); break;
        case SubKernel::Action::ASSEMBLE:
        case SubKernel::Action::REASSEMBLE:
        case SubKernel::Action::ITERATION:
        case SubKernel::Action::SOLUTION:
            runBoundaryKernel<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 2, 1>(operators, action); break;
        default: break;
        } break;
    case 3:
        switch (action) {
        case SubKernel::Action::PREPROCESS:
            setBoundaryKernel<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 3, 1>(operators, action); break;
        case SubKernel::Action::ASSEMBLE:
        case SubKernel::Action::REASSEMBLE:
        case SubKernel::Action::ITERATION:
        case SubKernel::Action::SOLUTION:
            runBoundaryKernel<Element::CODE::LINE2, 2, HeatTransferGPC::LINE2, 3, 1>(operators, action); break;
        default: break;
        } break;
    }
}

}
