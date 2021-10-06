
#include "scheme.h"
#include "steadystate.h"
#include "analysis/linearsystem/linearsystem.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/physicssolver/harmonic.h"

#include <cmath>

using namespace espreso;

AX_SteadyState::AX_SteadyState()
: K{}, f{}, x{}, dirichlet{}
{

}

AX_SteadyState::~AX_SteadyState()
{
	clear(K, f, x, dirichlet);
}

void AX_SteadyState::setTime(step::Time &time, double current)
{
	time.shift = current;
	time.start = 0;
	time.current = current;
	time.final = current;
}

void AX_SteadyState::init(AX_LinearSystem<double> *system)
{
	system->setMapping(K = system->assembler.A->copyPattern());
	system->setMapping(f = system->assembler.b->copyPattern());
	system->setMapping(x = system->assembler.x->copyPattern());
	system->setDirichletMapping(dirichlet = system->assembler.dirichlet->copyPattern());
}

void AX_SteadyState::composeSystem(step::Step &step, AX_LinearSystem<double> *system)
{
	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{K, f}\n");
		K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
		f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f").c_str());
		dirichlet->store(utils::filename(utils::debugDirectory(step) + "/scheme", "dirichlet").c_str());
	}

	system->solver.A->touched = true;
	system->solver.A->copy(K);
	system->solver.b->touched = true;
	system->solver.b->copy(f);
	system->solver.dirichlet->touched = true;
	system->solver.dirichlet->copy(dirichlet);
}

void AX_SteadyState::extractSolution(step::Step &step, AX_LinearSystem<double> *system)
{
	x->copy(system->solver.x);

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{x}\n");
		x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
	}
}
