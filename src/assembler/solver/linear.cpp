
#include "linear.h"
#include "../step.h"

using namespace espreso;

Linear::Linear(
		Mesh *mesh,
		std::vector<NewPhysics*> &physics,
		std::vector<NewInstance*> &instances,
		std::vector<LinearSolver*> &linearSolvers,
		store::ResultStore* store)
: Solver(mesh, physics, instances, linearSolvers, store)
{

}

void Linear::init()
{
	Step step;
	meshPreprocessing();
	assembleStiffnessMatrices(step);
	assembleB1(step);
	makeStiffnessMatricesRegular();
	assembleB0(step);

	initLinearSolver();
}

void Linear::solve(std::vector<std::vector<double> > &solution)
{
	Step step;
	startLinearSolver(step, solution);
}

void Linear::finalize()
{
	finalizeLinearSolver();
}

