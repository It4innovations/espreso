
#include "linear.h"
#include "../step.h"

using namespace espreso;

Linear::Linear(
		Mesh *mesh,
		std::vector<NewPhysics*> &physics,
		std::vector<Instance*> &instances,
		std::vector<LinearSolver*> &linearSolvers,
		store::ResultStore* store)
: Solver(mesh, physics, instances, linearSolvers, store)
{

}

void Linear::run(const Step &step)
{
	assembleStiffnessMatrices(step);
	assembleB1(step);
	makeStiffnessMatricesRegular();
	assembleB0(step);

	initLinearSolver();
	startLinearSolver(step);
	finalizeLinearSolver();
}
