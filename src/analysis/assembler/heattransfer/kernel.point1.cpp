
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runBoundary<Element::CODE::POINT1>(const step::Step &step, const step::Time &time, HeatTransferBoundaryOperators &operators, SubKernel::Action action)
{
    switch (info::mesh->dimension) {
    case 2:
        switch (action) {
        case SubKernel::Action::PREPROCESS:
            setNodeKernel<2>(operators, action); break;
        case SubKernel::Action::ASSEMBLE:
        case SubKernel::Action::REASSEMBLE:
        case SubKernel::Action::ITERATION:
        case SubKernel::Action::SOLUTION:
            runNodeKernel<2>(operators, action); break;
        default: break;
        } break;
    case 3:
        switch (action) {
        case SubKernel::Action::PREPROCESS:
            setNodeKernel<3>(operators, action); break;
        case SubKernel::Action::ASSEMBLE:
        case SubKernel::Action::REASSEMBLE:
        case SubKernel::Action::ITERATION:
        case SubKernel::Action::SOLUTION:
            runNodeKernel<3>(operators, action); break;
        default: break;
        } break;
    }
}

}
