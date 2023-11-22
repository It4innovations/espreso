
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::SQUARE4>(Action action, size_t interval)
{
	runConductivity<Element::CODE::SQUARE4, 4, HeatTransferGPC::SQUARE4, 2, 2>(subkernels[interval], action);

}

}
