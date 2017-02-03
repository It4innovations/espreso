
#include "linear.h"

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
	meshPreprocessing();
	assembleStiffnessMatrices();
	assembleB1();
	makeStiffnessMatricesRegular();
	assembleB0();

	initLinearSolver();
}

void Linear::solve(std::vector<std::vector<double> > &solution)
{
	startLinearSolver(solution);
}

void Linear::finalize()
{
	finalizeLinearSolver();
}

