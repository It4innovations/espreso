
#include "analysis/assembler/kernel/heattransfer/kernel.h"
#include "analysis/assembler/module/heattransfer.h"

namespace espreso {

template <>
void HeatTransfer::runGroup<Element::CODE::PRISMA15>(Action action, size_t interval)
{
	runConductivity<Element::CODE::PRISMA15, 15, HeatTransferGPC::PRISMA15, 3, 3>(subkernels[interval], action);
}

}
