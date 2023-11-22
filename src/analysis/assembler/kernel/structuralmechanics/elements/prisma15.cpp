
#include "analysis/assembler/kernel/structuralmechanics/kernel.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runGroup<Element::CODE::PRISMA15>(Action action, size_t interval, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
{
	runElasticity<Element::CODE::PRISMA15, 15, StructuralMechanicsGPC::PRISMA15, 3, 3, Behaviour::VOLUME>(subkernels[interval], action);
}

}
