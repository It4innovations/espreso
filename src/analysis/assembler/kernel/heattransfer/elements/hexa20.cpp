
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::HEXA20>(Action action, size_t interval)
{
	runConductivity<Element::CODE::HEXA20, 20, HeatTransferGPC::HEXA20, 3, 3>(subkernels[interval], action);
}

}
