
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::PYRAMID5>(Action action, size_t interval)
{
	runConductivity<Element::CODE::PYRAMID5, 5, HeatTransferGPC::PYRAMID5, 3, 3>(subkernels[interval], action);
}

}
