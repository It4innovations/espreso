
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::TRIANGLE6>(Action action, size_t interval)
{
	runConductivity<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2>(subkernels[interval], action);
}

}
