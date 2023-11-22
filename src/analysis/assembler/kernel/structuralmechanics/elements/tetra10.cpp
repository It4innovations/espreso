
#include "analysis/assembler/kernel/structuralmechanics/kernel.h"
#include "analysis/assembler/module/structuralmechanics.h"

namespace espreso {

template <>
void StructuralMechanics::runGroup<Element::CODE::TETRA10>(Action action, size_t interval, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
{
	runElasticity<Element::CODE::TETRA10, 10, StructuralMechanicsGPC::TETRA10, 3, 3, Behaviour::VOLUME>(subkernels[interval], action);
}

}
