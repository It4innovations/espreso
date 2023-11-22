
#include "analysis/assembler/kernel/structuralmechanics/kernel.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runGroup<Element::CODE::HEXA8>(Action action, size_t interval, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
{
	runElasticity<Element::CODE::HEXA8, 8, StructuralMechanicsGPC::HEXA8, 3, 3, Behaviour::VOLUME>(subkernels[interval], action);
}

}
