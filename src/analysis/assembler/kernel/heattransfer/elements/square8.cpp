
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::SQUARE8>(Action action, size_t interval)
{
	runConductivity<Element::CODE::SQUARE8, 8, HeatTransferGPC::SQUARE8, 2, 2>(subkernels[interval], action);
}

}
