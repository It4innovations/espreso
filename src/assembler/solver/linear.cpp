
#include "linear.h"
#include "../step.h"
#include "../instance.h"
#include "../physics/physics.h"
#include "../../basis/logging/logging.h"

using namespace espreso;

Linear::Linear(
		Mesh *mesh,
		Physics* physics,
		LinearSolver* linearSolver,
		output::Store* store,
		Matrices restriction)
: Solver("LINEAR", mesh, physics, linearSolver, store, restriction)
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
}

void Linear::init(Step &step)
{
	assembleMatrices(step, Matrices::M | Matrices::K | Matrices::f);
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
