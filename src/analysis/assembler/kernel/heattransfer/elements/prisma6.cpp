
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::PRISMA6>(Action action, size_t interval)
{
	runConductivity<Element::CODE::PRISMA6, 6, HeatTransferGPC::PRISMA6, 3, 3>(subkernels[interval], action);
}

}
