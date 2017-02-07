
#include "newtonrhapson.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"
#include "../physics/physics.h"

#include "../../basis/utilities/utils.h"
#include "../../mesh/settings/property.h"
#include "../../solver/generic/SparseMatrix.h"
#include "../../solver/generic/LinearSolver.h"

using namespace espreso;

NewtonRhapson::NewtonRhapson(
		Mesh *mesh,
		std::vector<Physics*> &physics,
		std::vector<Instance*> &instances,
		std::vector<LinearSolver*> &linearSolvers,
		store::ResultStore* store)
: Solver(mesh, physics, instances, linearSolvers, store)
{

}

void NewtonRhapson::run(Step &step)
{
	assembleStiffnessMatrices(step);
	assembleB1(step);
	makeStiffnessMatricesRegular(step);
	assembleB0(step);

	initLinearSolver();
	startLinearSolver();
	storeSolution(step);

	double precision = 1e-3;
	double normDeltaT = 2, normT = 1;

	while (normDeltaT / normT > 1e-3) {
		std::vector<std::vector<double> > T = instances.front()->primalSolution;
		step.solver++;

		instances.front() = new Instance(instances.front()->domains);
		instances.front()->DOFs = physics.front()->_instance->DOFs;
		instances.front()->primalSolution = physics.front()->_instance->primalSolution;
		instances.front()->solutions = physics.front()->_instance->solutions;
		physics.front()->_instance->solutions.resize(0);
		linearSolvers.front() = new LinearSolver(linearSolvers.front()->configuration, linearSolvers.front()->physics, linearSolvers.front()->constraints);
		physics.front()->_instance = instances.front();


		assembleStiffnessMatrices(step);
		subtractResidualForces(step);
		assembleB1(step);

		subtractSolutionFromB1c(step);

		makeStiffnessMatricesRegular(step);

		initLinearSolver();
		startLinearSolver();

		normDeltaT = norm(instances.front()->primalSolution);
		addToPrimar(0, T);
		normT = norm(instances.front()->primalSolution);
		storeSolution(step);
	}

	finalizeLinearSolver();
}



