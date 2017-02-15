
#include "linear.h"
#include "../step.h"

using namespace espreso;

Linear::Linear(
		Mesh *mesh,
		std::vector<Physics*> &physics,
		std::vector<Instance*> &instances,
		std::vector<LinearSolver*> &linearSolvers,
		store::ResultStore* store)
: Solver(mesh, physics, instances, linearSolvers, store)
{

}

void Linear::run(Step &step)
{
	assembleStiffnessMatrices(step);
	assembleB1(step);
	makeStiffnessMatricesRegular(step);
	assembleB0(step);

	initLinearSolver();
	startLinearSolver();
	processSolution(step);
	storeSolution(step);

	finalizeLinearSolver();
}
