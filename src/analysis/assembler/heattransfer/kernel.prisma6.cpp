
#include "kernel.h"
#include "element.h"

namespace espreso {

template <>
void runElement<Element::CODE::PRISMA6>(const step::Step &step, const step::Time &time, HeatTransferElementOperators &operators, SubKernel::Action action)
{
    switch (action) {
    case SubKernel::Action::PREPROCESS:
        setElementKernel<Element::CODE::PRISMA6, 6, HeatTransferGPC::PRISMA6, 3, 3>(operators, action); break;
    case SubKernel::Action::ASSEMBLE:
    case SubKernel::Action::REASSEMBLE:
    case SubKernel::Action::ITERATION:
    case SubKernel::Action::SOLUTION:
        runElementKernel<Element::CODE::PRISMA6, 6, HeatTransferGPC::PRISMA6, 3, 3>(step, time, operators, action); break;
    default: break;
    }
}

}
