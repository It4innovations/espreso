
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::TETRA4>(Action action, size_t interval)
{
	runConductivity<Element::CODE::TETRA4, 4, HeatTransferGPC::TETRA4, 3, 3>(subkernels[interval], action);
}

}
