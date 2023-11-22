
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::TRIANGLE3>(Action action, size_t interval)
{
	runConductivity<Element::CODE::TRIANGLE3, 3, HeatTransferGPC::TRIANGLE3, 2, 2>(subkernels[interval], action);
}

}
