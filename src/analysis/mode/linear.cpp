
#include "linear.h"

#include "esinfo/ecfinfo.h"

using namespace espreso;

void AX_Linear::init(AX_LinearSystem<double> *system)
{
	this->system = system;
}

bool AX_Linear::solve(AX_Scheme &scheme, AX_HeatTransfer &assembler)
{
	bool updateA = false, updateRHS = false, updateDirichlet = false;
	scheme.reassemble(assembler, updateA, updateRHS);
	updateDirichlet |= assembler.fillDirichlet(system->dirichlet);
	system->update(assembler, updateA, updateRHS, updateDirichlet);

	bool solved = system->solve();
	if (solved) {
		assembler.updateSolution(system->solver.x);
	}
	return solved;
}
