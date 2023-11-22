
#include "analysis/assembler/kernel/structuralmechanics/kernel.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runGroup<Element::CODE::HEXA20>(Action action, size_t interval, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
{
	runElasticity<Element::CODE::HEXA20, 20, StructuralMechanicsGPC::HEXA20, 3, 3, Behaviour::VOLUME>(subkernels[interval], action);
}

}
