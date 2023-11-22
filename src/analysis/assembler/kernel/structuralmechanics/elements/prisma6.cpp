
#include "analysis/assembler/kernel/structuralmechanics/kernel.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runGroup<Element::CODE::PRISMA6>(Action action, size_t interval, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
{
	runElasticity<Element::CODE::PRISMA6, 6, StructuralMechanicsGPC::PRISMA6, 3, 3, Behaviour::VOLUME>(subkernels[interval], action);
}

}
