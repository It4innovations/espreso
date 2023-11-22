
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::TETRA10>(Action action, size_t interval)
{
	runConductivity<Element::CODE::TETRA10, 10, HeatTransferGPC::TETRA10, 3, 3>(subkernels[interval], action);
}

}
