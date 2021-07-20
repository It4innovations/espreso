
#include "steadystate.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "output/output.h"

using namespace espreso;

void AX_SteadyState::init(AX_LinearSystem<double> *system, AX_HeatTransfer &assembler)
{
	system->setMapping(K = system->assembler.A->copyPattern());
	system->setMapping(f = system->assembler.b->copyPattern());
	assembler.setK(K);
	assembler.setRHS(f);
	assembler.init();
	this->system = system;
}

void AX_SteadyState::reassemble(AX_HeatTransfer &assembler, bool &updatedA, bool &updatedB)
{
	bool updatedK, updatedM, updatedRHS;
	assembler.next(updatedK, updatedM, updatedRHS);
	system->solver.A->fillData(K);
	system->solver.b->fillData(f);
	updatedA |= updatedK;
	updatedB |= updatedRHS;

	if (info::ecf->output.print_matrices) {
		K->store(utils::filename(utils::debugDirectory() + "/scheme", "K").c_str());
		f->store(utils::filename(utils::debugDirectory() + "/scheme", "f").c_str());
	}
}

void AX_SteadyState::solved()
{
	info::mesh->output->updateSolution();
}

