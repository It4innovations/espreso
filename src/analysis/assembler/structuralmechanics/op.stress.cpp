
#include "op.stress.h"
#include "analysis/assembler/structuralmechanics.h"

using namespace espreso;

Stress::Stress()
: enodes(info::mesh->elements->nodes->cbegin()),
  end(info::mesh->elements->nodes->cend()),
  multiplicity(nullptr),
  principalStress(nullptr),
  principalStrain(nullptr),
  componentStress(nullptr),
  componentStrain(nullptr),
  vonMisesStress(nullptr),
  vonMisesStrain(nullptr),
  principalStressAvg(nullptr),
  principalStrainAvg(nullptr),
  componentStressAvg(nullptr),
  componentStrainAvg(nullptr),
  vonMisesStressAvg(nullptr),
  vonMisesStrainAvg(nullptr)
{
    isconst = false;
    action = SubKernel::SOLUTION;
}


void Stress::activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, size_t interval, double *multiplicity)
{
    this->enodes = enodes;
    this->end = end;
    this->multiplicity = multiplicity;

    this->principalStress = StructuralMechanics::Results::principalStress->data.data() + info::mesh->elements->eintervals[interval].begin * info::mesh->dimension;
    this->principalStrain = StructuralMechanics::Results::principalStrain->data.data() + info::mesh->elements->eintervals[interval].begin * info::mesh->dimension;
    this->componentStress = StructuralMechanics::Results::componentStress->data.data() + info::mesh->elements->eintervals[interval].begin * 3 * (info::mesh->dimension - 1);
    this->componentStrain = StructuralMechanics::Results::componentStrain->data.data() + info::mesh->elements->eintervals[interval].begin * 3 * (info::mesh->dimension - 1);
    this->vonMisesStress  = StructuralMechanics::Results::vonMisesStress ->data.data() + info::mesh->elements->eintervals[interval].begin;
    this->vonMisesStrain  = StructuralMechanics::Results::vonMisesStrain ->data.data() + info::mesh->elements->eintervals[interval].begin;

    this->principalStressAvg = StructuralMechanics::Results::principalStressAvg->data.data();
    this->principalStrainAvg = StructuralMechanics::Results::principalStrainAvg->data.data();
    this->componentStressAvg = StructuralMechanics::Results::componentStressAvg->data.data();
    this->componentStrainAvg = StructuralMechanics::Results::componentStrainAvg->data.data();
    this->vonMisesStressAvg  = StructuralMechanics::Results::vonMisesStressAvg ->data.data();
    this->vonMisesStrainAvg  = StructuralMechanics::Results::vonMisesStrainAvg ->data.data();
    isactive = 1;
}
