
#include "analysis/assembler/kernel/structuralmechanics/kernel.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runGroup<Element::CODE::PYRAMID5>(Action action, size_t interval, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
{
	runElasticity<Element::CODE::PYRAMID5, 5, StructuralMechanicsGPC::PYRAMID5, 3, 3, Behaviour::VOLUME>(subkernels[interval], action);
}

}
