
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::HEXA8>(Action action, size_t interval)
{
	runConductivity<Element::CODE::HEXA8, 8, HeatTransferGPC::HEXA8, 3, 3>(subkernels[interval], action);
}

}
