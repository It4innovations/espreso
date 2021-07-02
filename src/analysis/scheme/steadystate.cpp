
#include "steadystate.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "output/output.h"

using namespace espreso;

void AX_SteadyState::init(AX_LinearSystem<double> *solver)
{
	K = solver->assembler.A->copyPattern();
	f = solver->assembler.b->copyPattern();
	mappingK = solver->mapping(K);
	mappingF = solver->mapping(f);
	this->solver = solver;
}

void AX_SteadyState::reassemble(AX_HeatTransfer &assembler, bool &A, bool &b)
{
	if (A |= assembler.fillK(K, mappingK)) {
		solver->solver.A->fillData(K);
	}
	if (b |= assembler.fillRHS(f, mappingF)) {
		solver->solver.b->fillData(f);
	}

	if (info::ecf->output.print_matrices) {
		K->store(utils::filename(utils::debugDirectory() + "/scheme", "K").c_str());
		f->store(utils::filename(utils::debugDirectory() + "/scheme", "f").c_str());
	}
}

void AX_SteadyState::solved()
{
	info::mesh->output->updateSolution();
}

