
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runElement<Element::CODE::HEXA20>(const step::Step &step, HeatTransferElementOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setElementKernel<Element::CODE::HEXA20, 20, HeatTransferGPC::HEXA20, 3, 3>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runElementKernel<Element::CODE::HEXA20, 20, HeatTransferGPC::HEXA20, 3, 3>(step, operators, action); break;
    default: break;
    }
}

}
