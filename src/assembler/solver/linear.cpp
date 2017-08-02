
#include "linear.h"
#include "../step.h"
#include "../instance.h"
#include "../physics/physics.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

Linear::Linear(
		Mesh *mesh,
		Physics* physics,
		FETISolver* linearSolver,
		output::Store* store,
		double duration,
		Matrices restriction)
: Solver("STEADY STATE", mesh, physics, linearSolver, store, duration, restriction)
{

}

void Linear::run(Step &step)
{
	ESINFO(PROGRESS1) << "Run " << _name << " solver for " << physics->name();

	_restriction &= ~Matrices::M;
	init(step);
	preprocess(step);
	solve(step);
	postprocess(step);
	finalize(step);

	instance->clear();
}

void Linear::init(Step &step)
{
	preprocessData(step);
	updateMatrices(step, Matrices::M | Matrices::K | Matrices::f);
	composeGluing(step, Matrices::B1);
	regularizeMatrices(step, Matrices::K);
	composeGluing(step, Matrices::B0);
}

void Linear::preprocess(Step &step)
{
	initLinearSolver(step);
}

void Linear::solve(Step &step)
{
	runLinearSolver(step);
}

void Linear::postprocess(Step &step)
{
	processSolution(step);
	storeSolution(step);
}

void Linear::finalize(Step &step)
{
	finalizeLinearSolver(step);
}
