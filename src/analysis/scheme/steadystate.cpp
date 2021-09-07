
#include "analysis/scheme/steadystate.h"
#include "analysis/linearsystem/linearsystem.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/physicssolver/harmonic.h"

#include <cmath>

using namespace espreso;

AX_SteadyState::AX_SteadyState()
: K{}, f{}, x{}
{

}

AX_SteadyState::~AX_SteadyState()
{
	if (K) {
		delete K;
	}
	if (f) {
		delete f;
	}
	if (x) {
		delete x;
	}
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
	system->setMapping(x = system->assembler.b->copyPattern());
}

void AX_SteadyState::composeSystem(step::Step &step, AX_LinearSystem<double> *system)
{
	system->solver.A->touched = true;
	system->solver.A->fillData(K);
	system->solver.b->touched = true;
	system->solver.b->fillData(f);

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{K, f}\n");
		K->store(utils::filename(utils::debugDirectory(step) + "/scheme", "K").c_str());
		f->store(utils::filename(utils::debugDirectory(step) + "/scheme", "f").c_str());
	}
}

void AX_SteadyState::extractSolution(step::Step &step, AX_LinearSystem<double> *system)
{
	x->fillData(system->solver.x);

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{x}\n");
		x->store(utils::filename(utils::debugDirectory(step) + "/scheme", "x").c_str());
	}
}
