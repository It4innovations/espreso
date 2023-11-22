
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::PYRAMID13>(Action action, size_t interval)
{
	runConductivity<Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13, 3, 3>(subkernels[interval], action);
}

}
