
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
		store::ResultStore* store)
: Solver("LINEAR", mesh, physics, linearSolver, store)
{

}

void Linear::run(Step &step)
{
	ESINFO(PROGRESS1) << "Run " << _name << " solver for " << physics->name();

	assembleMatrices(step, Matrices::K | Matrices::f);
	composeGluing(step, Matrices::B1);
	regularizeMatrices(step, Matrices::K);
	composeGluing(step, Matrices::B0);

	initLinearSolver();
	startLinearSolver();
	processSolution(step);
	storeSolution(step);

	finalizeLinearSolver();
}
